import os
import re
import ast
import sys
import time
import math
import shlex
import logging
import threading
import numpy as np
import mappy as mp
import pyabpoa as pa
from pyfaidx import Fasta
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional
from subprocess import PIPE, Popen,run,CalledProcessError
from collections import defaultdict, namedtuple, Counter

min_block_len = 30
min_supporting_reads = 2
max_gap_length = 200
min_cseq_length = 30
bp_cluster_distance = 30
MAX_GAP_TWO_ALIGNMENTS = 15
MAX_GAP_ONE_ALIGNMENT = 30
skip_chrM = True
only_keep_pcgene = True

MIN_GAP = -15
GENE_EXTENSION = 100
MIN_SOFT_CLIP_LENGTH = 30


ReadInfo = namedtuple("ReadInfo", ["read_id", "read_length", "read_start_offset", "read_end_offset", "strand", "seq","path",
                                   "ref_len", "ref_start_offset", "ref_end_offset","nr_matches", "block_len", "map_quality",
                                   "s1_score", "divergence","nm_errors", "s2", "cg", "tp"])

SamAlignmentInfo = namedtuple("SamAlignmentInfo", [
    "query_name", "flag", "ref_name", "start", "end", "mapq", "cigar", "tags","is_mapped",
    "is_primary", "is_secondary", "is_supplementary", "is_reverse",
    "q_len_cigar",       # Total query bases consumed by CIGAR ops in SEQ (S, M, I, =, X)
    "r_len_cigar",       # Total reference bases spanned by CIGAR (M, D, N, =, X)
    "q_aln_start",       # 0-based inclusive start of core alignment (M,I,=,X) on query segment (SEQ)
    "q_aln_end",         # 0-based exclusive end of core alignment (M,I,=,X) on query segment (SEQ)
    "s_start_len",       # Length of leading 'S' (soft-clip) operation on SEQ
    "s_end_len",         # Length of trailing 'S' (soft-clip) operation on SEQ
    "q_core_len",        # Sum of lengths of M, I, =, X operations (core alignment on SEQ)
    "h_start_len",       # Length of leading 'H' (hard-clip) operation
    "h_end_len"          # Length of trailing 'H' (hard-clip) operation
])

def parse_readinfo(line):
    parts = line.split("\t")
    read_id = parts[0]
    read_length = int(parts[1])
    read_start_offset = int(parts[2])
    read_end_offset = int(parts[3])
    strand = parts[4]
    seq = parts[5] 
    path = seq.lstrip("<").split("<") if seq[0] == '<' else seq.lstrip(">").split(">")
    ref_len = int(parts[6])
    ref_start_offset = int(parts[7])
    ref_end_offset = int(parts[8])
    nr_matches = int(parts[9])
    block_len = int(parts[10])
    map_quality = int(parts[11])
    tags = parts[12:]
    tag_dict = {tag.split(':')[0]: tag.split(':')[2] for tag in tags if ':' in tag}
    s1_score = int(tag_dict.get("s1", 0))
    divergence = float(tag_dict.get("dv", 1.0))
    nm_errors = int(tag_dict.get("NM", 0))
    s2 = int(tag_dict.get("s2", 0))
    cg = tag_dict.get("cg", " ")
    tp = tag_dict.get("tp", 0)
    #cm = tag_dict.get("cm", 0)
    return ReadInfo(read_id, read_length, read_start_offset, read_end_offset, strand, seq, path,
                    ref_len, ref_start_offset, ref_end_offset, nr_matches, block_len, map_quality,
                    s1_score, divergence, nm_errors, s2, cg, tp)

def parse_sam_line(line):
    if line.startswith('@'):
        return None
    fields = line.strip().split('\t')
    if len(fields) < 11:
        return None
    query_name = fields[0]
    flag = int(fields[1])
    ref_name = fields[2]
    start = int(fields[3]) if fields[3] != '*' else -1 # 1-based
    mapq = int(fields[4])
    cigar = fields[5]
    # Parse CIGAR details
    q_len_cigar, r_len_cigar, q_aln_start, q_aln_end, \
    s_start_len, s_end_len, q_core_len, \
    h_start_len, h_end_len = parse_cigar_for_query_details(cigar)
    ref_end_calculated = (start + r_len_cigar - 1) if (start != -1 and r_len_cigar > 0) else -1
    if cigar == '*':
        ref_end_calculated = -1
    tags = {}
    for tag_str in fields[11:]:
        try:
            tag_name, tag_type, tag_value = tag_str.split(':', 2)
            if tag_type == 'i':
                tags[tag_name] = int(tag_value)
            elif tag_type == 'f':
                tags[tag_name] = float(tag_value)
            else:
                tags[tag_name] = tag_value
        except ValueError:
            pass
    is_mapped = not bool(flag & 0x4)
    is_primary = not bool(flag & 0x100) and not bool(flag & 0x800) # Standard primary
    is_secondary = bool(flag & 0x100)
    is_supplementary = bool(flag & 0x800)
    is_reverse = bool(flag & 0x10)

    return SamAlignmentInfo(
        query_name=query_name, flag=flag, ref_name=ref_name, start=start, end=ref_end_calculated,
        mapq=mapq, cigar=cigar, tags=tags, is_mapped=is_mapped,is_primary=is_primary,
        is_secondary=is_secondary, is_supplementary=is_supplementary, is_reverse=is_reverse,
        q_len_cigar=q_len_cigar, r_len_cigar=r_len_cigar, q_aln_start=q_aln_start, q_aln_end=q_aln_end,
        s_start_len=s_start_len, s_end_len=s_end_len, q_core_len=q_core_len,
        h_start_len=h_start_len, h_end_len=h_end_len
    )

def setup_logging(log_path='fusion_detect_middle.log'):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s', 
        datefmt='%Y-%m-%d %H:%M:%S',  
        filename=log_path, 
        filemode='w'  
    )

def parse_fusion_file(file_path):
    fusion_dict = {}
    with open(file_path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or "--" not in line:
                continue
            fields = line.split('\t')
            if len(fields) != 3:
                continue
            gene_pair = fields[0]
            if fields[1] != '.':
                tags = ast.literal_eval(fields[1])
            fusion_dict[gene_pair] = tags
    return fusion_dict

def get_file_type(file_path: Path) -> str:
    extension = file_path.suffix.lower()
    if extension in ['.fa', '.fasta', '.fna']:
        return 'fasta'
    elif extension in ['.fq', '.fastq']:
        return 'fastq'
    else:
        print(f"Warning: Unrecognized file extension '{extension}' for '{file_path.name}'. Assuming FASTA format.")
        return 'fasta'

def convert_fastq_to_fasta_if_needed(input_file_path_str: str) -> Optional[str]:
    input_file_path = Path(input_file_path_str)
    if not input_file_path.is_file():
        sys.exit(f"Error: Input file not found: '{input_file_path_str}'")
    file_type = get_file_type(input_file_path)
    if file_type == 'fasta':
        return str(input_file_path)
    elif file_type == 'fastq':
        output_file_path = input_file_path.parent / (input_file_path.name + ".fasta")
        if output_file_path.exists():
            print(f"Output file '{output_file_path.name}' already exists. Skipping conversion for '{input_file_path.name}'.")
            return str(output_file_path)
        print(f"Detected FASTQ file: '{input_file_path.name}'. Converting to FASTA...")
        print(f"Output file will be: '{output_file_path.name}'")
        command = ['seqtk', 'seq', '-a', str(input_file_path)]
        try:
            with open(output_file_path, 'w') as outfile:
                run(
                    command,
                    stdout=outfile,
                    stderr=PIPE,
                    check=True,
                    text=True,
                    encoding='utf-8'
                )
            print(f"Successfully converted '{input_file_path.name}' to '{output_file_path.name}'")
            return str(output_file_path)

        except FileNotFoundError:
            print(f"Error: 'seqtk' command not found.", file=sys.stderr)
            print("Please ensure 'seqtk' is installed and available in your system's PATH.", file=sys.stderr)
            return None
        except CalledProcessError as e:
            print(f"Error: 'seqtk' failed during conversion of '{input_file_path.name}'.", file=sys.stderr)
            print(f"Command: {' '.join(e.cmd)}", file=sys.stderr)
            print(f"Return code: {e.returncode}", file=sys.stderr)
            print(f"Stderr:\n{e.stderr}", file=sys.stderr)
            if output_file_path.exists():
                try:
                    output_file_path.unlink()
                    print(f"Removed incomplete output file: '{output_file_path.name}'")
                except OSError as del_err:
                    print(f"Warning: Could not remove incomplete output file '{output_file_path.name}': {del_err}", file=sys.stderr)
            return None
        except Exception as e:
            print(f"An unexpected error occurred during conversion: {e}", file=sys.stderr)
            if output_file_path.exists() and output_file_path.stat().st_size == 0:
                try:
                   output_file_path.unlink()
                   print(f"Removed potentially empty output file: '{output_file_path.name}'")
                except OSError as del_err:
                    print(f"Warning: Could not remove potentially empty output file '{output_file_path.name}': {del_err}", file=sys.stderr)
            return None
    return None

def check_fusion_pair(geneA, geneB, fusion_annot_dict):
    special_normal_tags = [
        "ConjoinG", "Babiceanu_Normal", "Greger_Normal", "GTEx_recurrent_StarF2019", "DGD_PARALOGS", "BodyMap"
    ]
    pair_keys = [f"{geneA}--{geneB}", f"{geneB}--{geneA}"]
    for key in pair_keys:
        if key in fusion_annot_dict:
            tags = fusion_annot_dict[key]
            if any(tag in tags for tag in special_normal_tags):
                return "NORMAL_FUSION"
            else:
                return "PASS"
    return None

def load_gene_positions(file_path):
    gene_positions = {}
    try:
        with open(file_path, 'r') as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    gene_key = parts[0]
                    gene_name = gene_key.split("_")[0]
                    chrom = gene_key.split("_")[-1]
                    strand = parts[1]
                    type = parts[2]
                    start = int(parts[3])
                    end = int(parts[4])
                    Fexons = list(parts[5].split(';'))
                    Bexons = list(parts[6].split(';'))
                    gene_positions[gene_key] = {'gene_name': gene_name,'chrom': chrom, 'start': start, 'strand': strand, 'type': type, 'end': end, 'Fexons': Fexons, 'Bexons': Bexons}
    except Exception as e:
        logging.error(f"Error reading gene position file: {e}")
    return gene_positions

def is_valid_fusion(gene1_genome_strand, gene1_graph_strand, gene2_genome_strand, gene2_graph_strand):
    transcript_dir1 = (gene1_genome_strand == "+") == (gene1_graph_strand == ">")
    transcript_dir2 = (gene2_genome_strand == "+") == (gene2_graph_strand == ">")
    return transcript_dir1 == transcript_dir2

def validate_fusion(gene_positions,upper_node,lower_node,only_keep_pcgene = True):
    pos_key_up = f"{upper_node.path[0].split(':')[1]}_{upper_node.path[0].split(':')[0]}"
    pos_key_lo = f"{lower_node.path[0].split(':')[1]}_{lower_node.path[0].split(':')[0]}"
    gene1_genome_strand = gene_positions[pos_key_up]["strand"]
    gene2_genome_strand = gene_positions[pos_key_lo]["strand"]
    gene1_type = gene_positions[pos_key_up]["type"]
    gene2_type = gene_positions[pos_key_lo]["type"]
    gene1_graph_strand = upper_node.seq[0]
    gene2_graph_strand = lower_node.seq[0]
    if only_keep_pcgene:
        if gene1_type == 'protein_coding' and gene2_type == 'protein_coding':
            return is_valid_fusion(gene1_genome_strand, gene1_graph_strand, gene2_genome_strand, gene2_graph_strand)
        else:
            return False
    else:
        return is_valid_fusion(gene1_genome_strand, gene1_graph_strand, gene2_genome_strand, gene2_graph_strand)

def get_all_genes_from_path(path_str):
    parts = path_str.split(":")
    if len(parts) < 2:
        return set()
    gene_part = parts[1]
    if gene_part.replace("-", "").isdigit():
        return set()
    if "--" in gene_part:
        return set(gene_part.split("--"))
    else:
        return {gene_part}

def merge_fusion_records(sorted_cluster):
    if not sorted_cluster:
        return []
    merged_records = []
    current_record = sorted_cluster[0]
    current_gene = current_record.path[0].split(':')[1]
    current_chrom = current_record.path[0].split(":")[0]

    for next_record in sorted_cluster[1:]:
        next_gene = next_record.path[0].split(':')[1]
        next_chrom = next_record.path[0].split(":")[0]
        gap = next_record.read_start_offset - current_record.read_end_offset
        current_path_set = set(current_record.path)
        next_path_set = set(next_record.path)
        has_intersection = not current_path_set.isdisjoint(next_path_set)
        if current_gene == next_gene and current_chrom == next_chrom and gap >= -5:
            if current_record.seq[0] == next_record.seq[0] and not has_intersection:
                merged_seq = current_record.seq + next_record.seq
                merged_path = current_record.path + next_record.path
                merged_s1 = current_record.s1_score + next_record.s1_score
                merged_nm = current_record.nm_errors + next_record.nm_errors
                merged_alignment_length = current_record.block_len + next_record.block_len
                merged_record = ReadInfo(
                    read_id=current_record.read_id,
                    read_length=current_record.read_length,
                    read_start_offset=current_record.read_start_offset,
                    read_end_offset=next_record.read_end_offset,
                    strand=current_record.strand,
                    seq=merged_seq,
                    path=merged_path,
                    ref_len=current_record.ref_len,
                    ref_start_offset=current_record.ref_start_offset,
                    ref_end_offset=next_record.ref_end_offset,
                    nr_matches=current_record.nr_matches + next_record.nr_matches,
                    block_len=merged_alignment_length,
                    map_quality=max(current_record.map_quality, next_record.map_quality),
                    s1_score=merged_s1,
                    divergence=(current_record.divergence + next_record.divergence) / 2,
                    nm_errors=merged_nm,
                    s2=current_record.s2 + next_record.s2,
                    cg=current_record.cg + f'{gap}D' + next_record.cg,
                    tp=current_record.tp
                )

                current_record = merged_record
            else:
                merged_records.append(current_record)
                current_record = next_record
                current_gene = next_gene
                current_chrom = next_chrom
                continue
        else:
            merged_records.append(current_record)
            current_record = next_record
            current_gene = next_gene
            current_chrom = next_chrom
    merged_records.append(current_record)
    return merged_records

def check_gene_similarity(gene1, gene2):
    if (gene1 and gene1.startswith("ENSG")) or (gene2 and gene2.startswith("ENSG")):
        return False
    if not gene1 or not gene2:
        return False
    if gene1 == gene2:
        return True
    min_len = min(len(gene1), len(gene2))
    if (gene1 in gene2 or gene2 in gene1) and min_len >= 3:
        return True
    return False

def check_genes_similarity(sorted_cluster):
    num_records = len(sorted_cluster)
    if num_records == 0 or num_records == 1:
        return False

    all_genes_sets = []
    maingene_list = []
    for record in sorted_cluster:
        genes_in_record = set()
        if hasattr(record, 'path') and record.path is not None:
            for path_item in record.path:
                genes_in_record.update(get_all_genes_from_path(path_item))
        main_gene = None
        if hasattr(record, 'seq'):
            main_gene = record.path[0].split(':')[1]
        all_genes_sets.append(genes_in_record)
        maingene_list.append(main_gene)
    common_genes_intersection = all_genes_sets[0].copy()
    for i in range(1, num_records):
        common_genes_intersection.intersection_update(all_genes_sets[i])
    if common_genes_intersection: 
        return True
    all_main_gene_pairs_similar = True
    if num_records >= 2: 
        for i in range(num_records):
            for j in range(i + 1, num_records):
                gene1 = maingene_list[i]
                gene2 = maingene_list[j]
                if not check_gene_similarity(gene1, gene2):
                    all_main_gene_pairs_similar = False
                    break
            if not all_main_gene_pairs_similar:
                break
        if all_main_gene_pairs_similar:
            return True
    return False

def truncate_to_first_two_genes(merged_cluster):
    all_genes = set()
    for record in merged_cluster:
        gene = record.path[0].split(':')[1]
        all_genes.add(gene)
    if len(all_genes) > 2:
        return merged_cluster, None
    found_genes = set()
    result = []
    for i, record1 in enumerate(merged_cluster):
        gene = record1.path[0].split(':')[1]
        result.append(record1)
        found_genes.add(gene)
        if i == 1:
            break
    if len(found_genes) == 1:
        return merged_cluster, None

    truncation_point = (merged_cluster[1].read_end_offset + merged_cluster[2].read_start_offset) // 2

    return result, truncation_point

def check_overlaps(sorted_cluster,overlap = 1):
    def check_overlap(record1, record2, overlap):
        overlap_start = max(record1.read_start_offset, record2.read_start_offset)
        overlap_end = min(record1.read_end_offset, record2.read_end_offset)
        overlap_length = overlap_end - overlap_start
        return overlap_length >= overlap
    for i in range(len(sorted_cluster) - 1):
        for j in range(i + 1, len(sorted_cluster)):
            if check_overlap(sorted_cluster[i], sorted_cluster[j], overlap):
                return True
    return False

def gene2_dict():
    return {
        0: [],
        1: []
    }

def compute_breakpoints_with_position(node,position):
    if position == 'upper':
        if node.seq[0] == '>':
            bp =  node.path[-1].split(':')[-1].split('-')[-1]
        else:
            bp =  node.path[-1].split(':')[-1].split('-')[0]
    elif position == 'lower':
        if node.seq[0] == '>':
            bp =  node.path[0].split(':')[-1].split('-')[0]
        else:
            bp =  node.path[0].split(':')[-1].split('-')[-1]
    else:
        bp =  node.path[-1].split(':')[-1].split('-')[-1]
    return bp

def extract_position_info(path_str):
    try:
        parts = path_str.split(':')
        chrom = parts[0] 
        if len(parts) > 1:
            position_part = parts[-1]
            return f"{chrom}:{position_part}"
        return path_str
    except:
        return path_str

def build_total_reference_fasta(gene_set, gene_positions, reference_genome_path, output_fasta, extension=100):
    region_params = []
    header_map = []
    for gene, chrom in gene_set:
        pos_key = f"{gene}_{chrom}"
        start = max(1, gene_positions[pos_key]['start'] - extension)
        end = gene_positions[pos_key]['end'] + extension
        region = f"{chrom}:{start}-{end}"
        region_params.append(region)
        header_map.append((region, f"{chrom}:{gene}:{start}-{end}"))
    region_file = output_fasta + ".regions.txt"
    with open(region_file, "w") as f:
        for region in region_params:
            f.write(region + "\n")
    cmd = f"samtools faidx {reference_genome_path} -r {region_file}"
    result = run(cmd, stdout=PIPE, shell=True, check=True, encoding='utf-8')
    fasta_lines = result.stdout.splitlines()

    header_map_dict = dict(header_map)
    new_lines = []
    for line in fasta_lines:
        if line.startswith('>'):
            origin = line[1:]
            new = header_map_dict.get(origin, origin)
            new_lines.append(f">{new}")
        else:
            new_lines.append(line)
    with open(output_fasta, 'w') as f:
        f.write('\n'.join(new_lines) + '\n')
        
def validate_fusion_genes_with_mappy(gene2_cluster, total_ref_fasta_path, fasta_sequences, threads=8):
    validate_gene2_cluster = defaultdict(lambda: {0: [], 1: []})
    try:
        aligner = mp.Aligner(fn_idx_in=total_ref_fasta_path, preset="splice", k=16,min_dp_score=20, min_chain_score=10,n_threads=threads)
    except Exception as e:
        return validate_gene2_cluster

    for fusion_key, read_info_dict in gene2_cluster.items():
        gene_names = fusion_key.split(':')
        if len(gene_names) != 2:
            continue
        gene1, gene2 = gene_names
        read_ids = set()
        read_details = {} 

        for i in range(2):
            for node in read_info_dict[i]:
                read_id = node.read_id
                read_ids.add(read_id)
                if read_id not in read_details:
                    read_details[read_id] = {
                        'length': node.read_length,
                        'nodes': [None, None]
                    }
                read_details[read_id]['nodes'][i] = node

        mappy_alignments = defaultdict(list)
        for read_id, details in read_details.items():
            try:
                seq = str(fasta_sequences[read_id])[:details['length']]
                mappy_alignments[read_id] = list(aligner.map(seq))
            except Exception as e:
                continue

        valid_read_ids = set()
        for read_id, alignments in mappy_alignments.items():
            if read_id not in read_details:
                continue
            gaf_node1 = read_details[read_id]['nodes'][0]
            gaf_node2 = read_details[read_id]['nodes'][1]
            read_length = read_details[read_id]['length']
            choose_alignments = [aln for aln in alignments if aln.ctg.split(':')[1] in (gene1, gene2)]
            gene1_alignments = [aln for aln in choose_alignments if aln.ctg.split(':')[1] == gene1]
            gene2_alignments = [aln for aln in choose_alignments if aln.ctg.split(':')[1] == gene2]

            if len(gene1_alignments) >= 1 and len(gene2_alignments) >= 1:
                if gene1_alignments and gene2_alignments:
                    aln1 = max(gene1_alignments, key=lambda x: x.q_en)
                    aln2 = min(gene2_alignments, key=lambda x: x.q_st)
                    gap = aln2.q_st - aln1.q_en
                    if (MIN_GAP <= gap <= MAX_GAP_TWO_ALIGNMENTS) :
                        valid_read_ids.add(read_id)
            elif (len(gene1_alignments) >= 1 and len(gene2_alignments) == 0) or (len(gene1_alignments) == 0 and len(gene2_alignments) >= 1):
                sam_record = max(choose_alignments, key=lambda x: x.blen)
                gaf_of_aligned_gene = gaf_node1 if sam_record.ctg.split(':')[1] == gene1 else gaf_node2
                gaf_of_missing_gene = gaf_node2 if gaf_of_aligned_gene == gaf_node1 else gaf_node1
                if gaf_node1 and gaf_node2:
                    gap = gaf_node2.read_start_offset - gaf_node1.read_end_offset
                    if MIN_GAP <= gap <= MAX_GAP_ONE_ALIGNMENT:
                        map_len = sam_record.q_en - sam_record.q_st
                        if abs(read_length - map_len) <= MIN_SOFT_CLIP_LENGTH:
                            continue
                        missing_gaf_segment_len_on_read = gaf_of_missing_gene.read_end_offset - gaf_of_missing_gene.read_start_offset
                        if gaf_of_missing_gene == gaf_node2:
                            clip_len_to_check = read_length - sam_record.q_en
                        else:
                            clip_len_to_check = sam_record.q_st
                        if clip_len_to_check >= MIN_SOFT_CLIP_LENGTH and \
                           abs(clip_len_to_check - missing_gaf_segment_len_on_read) <= MIN_SOFT_CLIP_LENGTH:
                            valid_read_ids.add(read_id)
            else:
                pass
        for i in range(2):
            for node in read_info_dict[i]:
                if node.read_id in valid_read_ids:
                    validate_gene2_cluster[fusion_key][i].append(node)

    return validate_gene2_cluster

def find_nearest_breakpoint(bp, breakpoints, threshold = 30):
    nearest = None
    min_diff = threshold + 1
    for b in breakpoints:
        try:
            b_int = int(b)
            diff = abs(bp - b_int)
            if diff <= threshold and diff < min_diff:
                min_diff = diff
                nearest = b_int
        except Exception:
            continue
    return nearest

def compute_breakpoints_with_position_correct(node, position):
    td = 20
    if position == 'upper':
        distance = node.ref_len - node.ref_end_offset + 1
        if node.seq[0] == '>':
            if distance <= td:
                bp = int(node.path[-1].split(':')[-1].split('-')[-1])
            else:
                bp = int(node.path[-1].split(':')[-1].split('-')[-1]) - distance
        else:
            if distance <= td:
                bp = int(node.path[-1].split(':')[-1].split('-')[0])
            else:
                bp = int(node.path[-1].split(':')[-1].split('-')[0]) + distance
    elif position == 'lower':
        if node.seq[0] == '>':
            if node.ref_start_offset <= td:
                bp = int(node.path[0].split(':')[-1].split('-')[0])
            else:
                bp = int(node.path[0].split(':')[-1].split('-')[0]) + node.ref_start_offset
        else:
            if node.ref_start_offset <= td:
                bp = int(node.path[0].split(':')[-1].split('-')[-1])
            else:
                bp = int(node.path[0].split(':')[-1].split('-')[-1]) - node.ref_start_offset
    else:
        bp = -1
    return bp

class FusionAggregator:
    def __init__(self):
        self.candidates = []
        self.read_alignments = {}  # {read_id: {'upper': {'start': start1, 'end': end1, 'strand': strand1, 'read_length': len1},
                                   #            'lower': {'start': start2, 'end': end2, 'strand': strand2, 'read_length': len2}}}

    def _calculate_base_weight(self, map_quality):
        if map_quality >= 30:
            return 10 + map_quality * 0.15
        elif 10 <= map_quality < 30:
            return (map_quality // 10) * 3.5
        else:
            return map_quality * 0.2

    def _extract_path_info(self, node, position):
        direction = node.seq[0] if hasattr(node, 'seq') and node.seq else '>'
        simplified_path = []
        if hasattr(node, 'path') and node.path:
            for path_item in node.path:
                simplified_path.append(extract_position_info(path_item))
            chrom = node.path[0].split(':')[0]
        else:
            simplified_path = []
            chrom = 'unknown'

        path = simplified_path if direction == '>' else list(reversed(simplified_path))
        chrom = node.path[0].split(':')[0]
        gene_name = node.path[0].split(':')[1]
        breakpoint = compute_breakpoints_with_position(node, position)

        read_start = getattr(node, 'read_start_offset', 0)
        read_end = getattr(node, 'read_end_offset', 0)
        read_length = getattr(node, 'read_length', 0)
        strand = direction  

        return {
            'direction': direction,
            'path': path,
            'chrom': chrom,
            'gene_name': gene_name,
            'path_length': len(path), 
            'map_quality': getattr(node, 'map_quality', 0),
            'breakpoint': breakpoint,
            'read_start': read_start,
            'read_end': read_end,
            'read_length': read_length,
            'strand': strand
        }

    def process_read_pair(self, upper_node, lower_node):

        upper_info = self._extract_path_info(upper_node, 'upper')
        lower_info = self._extract_path_info(lower_node, 'lower')
        upper_m_bp = compute_breakpoints_with_position_correct(upper_node, 'upper')
        lower_m_bp = compute_breakpoints_with_position_correct(lower_node, 'lower')
        read_id = upper_node.read_id 
        self.read_alignments[read_id] = {
            'upper': {
                'start': upper_info['read_start'],
                'end': upper_info['read_end'],
                'strand': upper_info['strand'],
                'read_length': upper_info['read_length'],
                'm_bp': upper_m_bp
            },
            'lower': {
                'start': lower_info['read_start'],
                'end': lower_info['read_end'],
                'strand': lower_info['strand'],
                'read_length': lower_info['read_length'],
                'm_bp': lower_m_bp
            }
        }
        upper_weight = self._calculate_base_weight(upper_info['map_quality'])
        lower_weight = self._calculate_base_weight(lower_info['map_quality'])
        upper_path_len = upper_info['path_length']
        lower_path_len = lower_info['path_length']
        length_bonus = (1 + math.log1p(upper_path_len + lower_path_len)) * 0.5
        if upper_info['map_quality'] >= 30 and lower_info['map_quality'] >= 30:
            quality_bonus = 1.5
        elif upper_info['map_quality'] >= 10 and lower_info['map_quality'] >= 10:
            quality_bonus = 1.2
        else:
            quality_bonus = 0.8
        combined_score = ((upper_weight + lower_weight) ** 0.5) * length_bonus * quality_bonus
        metadata = {
            'upper_dir': upper_info['direction'],
            'upper_path': upper_info['path'], # list
            'upper_chrom': upper_info['chrom'],
            'upper_bp': upper_info['breakpoint'],
            'upper_quality': upper_info['map_quality'],
            'upper_read_start': upper_info['read_start'],
            'upper_read_end': upper_info['read_end'],
            'lower_dir': lower_info['direction'],
            'lower_path': lower_info['path'], # list
            'lower_chrom': lower_info['chrom'],
            'lower_bp': lower_info['breakpoint'],
            'lower_quality': lower_info['map_quality'],
            'lower_read_start': lower_info['read_start'],
            'lower_read_end': lower_info['read_end'],
            'read_id': read_id,
            'read_length': max(upper_info['read_length'], lower_info['read_length'])
        }
        self.candidates.append((metadata, combined_score))

    def _cluster(self, max_breakpoint_distance=15):
        if not self.candidates:
            return {}
        return self._process_clusters(self.candidates, max_breakpoint_distance)

    def _process_clusters(self, candidates, bp_cluster_distance=0):
        if not candidates:
            return {}
        breakpoint_clusters = []
        for meta, score in candidates:
            upper_bp = meta['upper_bp']
            lower_bp = meta['lower_bp']
            upper_chrom = meta['upper_chrom']
            lower_chrom = meta['lower_chrom']
            if upper_bp is None or lower_bp is None: continue
            try:
                int_upper_bp = int(upper_bp)
                int_lower_bp = int(lower_bp)
            except (ValueError, TypeError): continue

            found_cluster = False
            for cluster in breakpoint_clusters:
                example_meta = cluster[0][0] 
                try:
                    example_upper_bp = int(example_meta['upper_bp'])
                    example_lower_bp = int(example_meta['lower_bp'])
                except (ValueError, TypeError, KeyError):
                    continue
                example_upper_chrom = example_meta['upper_chrom']
                example_lower_chrom = example_meta['lower_chrom']

                if upper_chrom != example_upper_chrom or lower_chrom != example_lower_chrom: continue

                upper_bp_match = abs(int_upper_bp - example_upper_bp) <= bp_cluster_distance
                lower_bp_match = abs(int_lower_bp - example_lower_bp) <= bp_cluster_distance

                if upper_bp_match and lower_bp_match:
                    cluster.append((meta, score))
                    found_cluster = True
                    break

            if not found_cluster:
                breakpoint_clusters.append([(meta, score)])
        final_result = {}
        bp_group_counter = 0 
        for bp_cluster in breakpoint_clusters:
            if len(bp_cluster) > 1:
                path_subclusters = defaultdict(list)
                for meta, score in bp_cluster:
                    path_key = (tuple(meta['upper_path']), tuple(meta['lower_path']))
                    path_subclusters[path_key].append((meta, score))
                path_group_counter = 0
                for path_key, sub_cluster in path_subclusters.items():
                    if not sub_cluster: continue 
                    sub_cluster.sort(key=lambda x: x[1], reverse=True)
                    best_meta, best_score = sub_cluster[0]
                    read_ids = [meta['read_id'] for meta, _ in sub_cluster]
                    unique_read_ids = set(read_ids)
                    total_score = sum(score for _, score in sub_cluster)
                    upper_bps = [int(meta['upper_bp']) for meta, _ in bp_cluster if meta['upper_bp'] is not None]
                    lower_bps = [int(meta['lower_bp']) for meta, _ in bp_cluster if meta['lower_bp'] is not None]
                    if not upper_bps or not lower_bps: continue
                    rep_upper_bp = Counter(upper_bps).most_common(1)[0][0]
                    rep_lower_bp = Counter(lower_bps).most_common(1)[0][0]
                    result_key = f"{best_meta['upper_chrom']}:{rep_upper_bp}-{best_meta['lower_chrom']}:{rep_lower_bp}_bpGroup_{bp_group_counter}_pathGroup_{path_group_counter}"
                    read_alignment_info = {read_id: self.read_alignments[read_id] for read_id in unique_read_ids if read_id in self.read_alignments}
                    upper_data = {'direction': best_meta['upper_dir'], 'path': best_meta['upper_path'], 'chrom': best_meta['upper_chrom'], 'breakpoint': rep_upper_bp, 'score': total_score / 2}
                    lower_data = {'direction': best_meta['lower_dir'], 'path': best_meta['lower_path'], 'chrom': best_meta['lower_chrom'], 'breakpoint': rep_lower_bp, 'score': total_score / 2}
                    final_result[result_key] = {
                        'upper': upper_data, 'lower': lower_data, 'score': total_score,
                        'supporting_reads': len(unique_read_ids), 'read_ids': list(unique_read_ids),
                        'bp_group': bp_group_counter, 
                        'path_group': path_group_counter, 
                        'read_alignment_info': read_alignment_info
                    }
                    path_group_counter += 1
            else:
                if not bp_cluster: continue 
                meta, score = bp_cluster[0]
                upper_bp = meta['upper_bp']
                lower_bp = meta['lower_bp']
                if upper_bp is None or lower_bp is None: continue
                try:
                    rep_upper_bp = int(upper_bp)
                    rep_lower_bp = int(lower_bp)
                except (ValueError, TypeError): continue
                read_id = meta['read_id']
                unique_read_ids = {read_id}
                total_score = score
                result_key = f"{meta['upper_chrom']}:{rep_upper_bp}-{meta['lower_chrom']}:{rep_lower_bp}_bpGroup_{bp_group_counter}_pathGroup_0" # Path group 默认为 0
                read_alignment_info = {read_id: self.read_alignments[read_id]} if read_id in self.read_alignments else {}
                upper_data = {'direction': meta['upper_dir'], 'path': meta['upper_path'], 'chrom': meta['upper_chrom'], 'breakpoint': rep_upper_bp, 'score': total_score / 2}
                lower_data = {'direction': meta['lower_dir'], 'path': meta['lower_path'], 'chrom': meta['lower_chrom'], 'breakpoint': rep_lower_bp, 'score': total_score / 2}
                final_result[result_key] = {
                    'upper': upper_data, 'lower': lower_data, 'score': total_score,
                    'supporting_reads': 1, 'read_ids': list(unique_read_ids),
                    'bp_group': bp_group_counter, 
                    'path_group': 0,
                    'read_alignment_info': read_alignment_info
                }
            bp_group_counter += 1

        return final_result
    
    def format_fusion_results(self, clusters, gene_names):
        final_fusion = {}
        if not clusters:
            return final_fusion
        for cluster_key, cluster_data in clusters.items():
            bp_group = cluster_data.get('bp_group', -1)
            path_group = cluster_data.get('path_group', -1)
            norm_key = f"{':'.join(gene_names)}_bpGroup_{bp_group}_pathGroup_{path_group}"
            if norm_key not in final_fusion:
                final_fusion[norm_key] = {}
            upper_gene_key = gene_names[0]
            upper_data = cluster_data['upper']
            final_fusion[norm_key][upper_gene_key] = {
                "Path": [upper_data['direction']] + upper_data['path'],
                "Reads": cluster_data['read_ids'],
                "Chrom": upper_data['chrom'],
                "Breakpoint": upper_data['breakpoint'],
                "SupportingReads": cluster_data['supporting_reads'],
                "ConfidenceScore": cluster_data['upper']['score'] * 2,
                "BreakpointGroup": bp_group,
                "PathGroup": path_group,
                "ReadAlignmentInfo": cluster_data.get('read_alignment_info', {})
            }

            lower_gene_key = gene_names[1]
            lower_data = cluster_data['lower']
            final_fusion[norm_key][lower_gene_key] = {
                "Path": [lower_data['direction']] + lower_data['path'],
                "Reads": cluster_data['read_ids'],
                "Chrom": lower_data['chrom'],
                "Breakpoint": lower_data['breakpoint'],
                "SupportingReads": cluster_data['supporting_reads'],
                "ConfidenceScore": cluster_data['lower']['score'] * 2,
                "BreakpointGroup": bp_group,
                "PathGroup": path_group,
                "ReadAlignmentInfo": cluster_data.get('read_alignment_info', {})
            }
        return final_fusion

def proc_fuison_cluster(cluster, bp_cluster_distance=15):
    final_fusion = {}

    for fusion_key, gene_data in cluster.items():
        gene_names = fusion_key.split(':')
        if len(gene_names) != 2: 
            continue
        
        aggregator = FusionAggregator()
        upper_nodes = gene_data.get(0, [])
        lower_nodes = gene_data.get(1, [])
        node_pairs = defaultdict(dict)
        for node in upper_nodes:
            if hasattr(node, 'read_id'):
                 node_pairs[node.read_id][0] = node
        for node in lower_nodes:
            if hasattr(node, 'read_id'):
                 node_pairs[node.read_id][1] = node
        processed_read_ids = set()
        for read_id, nodes in node_pairs.items():
            if read_id not in processed_read_ids and 0 in nodes and 1 in nodes:
                aggregator.process_read_pair(nodes[0], nodes[1])
                processed_read_ids.add(read_id)
        if not aggregator.candidates:
            continue
        
        clustered_results = aggregator._cluster(bp_cluster_distance)

        if not clustered_results:
            continue
        fusion_result = aggregator.format_fusion_results(clustered_results, gene_names)
        for key, value in fusion_result.items():
            final_fusion[key] = value 

    return final_fusion

class Gene:
    def __init__(self, name, chrom, breakpoint, strand, start=0, end=0):
        self.name = name
        self.chrom = chrom
        self.breakpoint = breakpoint
        self.strand = strand
        self.start = start
        self.end = end
    def __str__(self):
        return f"{self.name}({self.chrom}:{self.breakpoint})"

class FusionCandidate:
    def __init__(self, upstream_gene, downstream_gene, read_alignment_info, supporting_reads=None, group_id=0, path_group_id=0):
        self.upstream_gene = upstream_gene
        self.downstream_gene = downstream_gene
        self.supporting_reads = supporting_reads or []
        self.group_id = group_id
        self.path_group_id = path_group_id
        self.confidence_level = 1
        self.upper_consensus_sequence = ""
        self.lower_consensus_sequence = ""
        self.consensus_sequence = ""
        self.minigraph_upper_alignments = []
        self.minigraph_lower_alignments = []
        self.minimap2_full_alignments = [] 
        self.totalreads = 0
        self.cverror = 0
        self.read_alignment_info = read_alignment_info
    def get_id(self):
        return f"{self.upstream_gene.name}_{self.downstream_gene.name}_group{self.group_id}_path{self.path_group_id}"
    def clone(self):
        candidate_clone = FusionCandidate(
            upstream_gene=self.upstream_gene,
            downstream_gene=self.downstream_gene,
            read_alignment_info=self.read_alignment_info.copy() if self.read_alignment_info else {}, # 拷贝字典
            supporting_reads=list(self.supporting_reads),
            path_group_id=self.path_group_id
        )
        candidate_clone.confidence_level=self.confidence_level 
        return candidate_clone
    def get_bp_group_id(self):
        return f"{self.upstream_gene.name}_{self.downstream_gene.name}_group{self.group_id}"
    def get_read_count(self):
        return len(self.supporting_reads)
    def add_upper_consensus_sequence(self, sequence): self.upper_consensus_sequence = sequence
    def add_lower_consensus_sequence(self, sequence): self.lower_consensus_sequence = sequence
    def add_consensus_sequence(self, sequence): self.consensus_sequence = sequence # 新增
    def add_minigraph_upper_alignment(self, alignment): self.minigraph_upper_alignments.append(alignment)
    def add_minigraph_lower_alignment(self, alignment): self.minigraph_lower_alignments.append(alignment)
    def add_minimap2_full_alignment(self, alignment): self.minimap2_full_alignments.append(alignment)
    def get_unique_key(self):
        return tuple([
            f"{self.upstream_gene.name}:{self.upstream_gene.breakpoint}:{self.upstream_gene.strand}",
            f"{self.downstream_gene.name}:{self.downstream_gene.breakpoint}:{self.downstream_gene.strand}"
        ])
    def merge_with(self, other_fusion):
        all_reads = set(self.supporting_reads).union(set(other_fusion.supporting_reads))
        self.supporting_reads = list(all_reads)
        return self
    def get_modify_bp(self, gene_positions, pos_key_up, pos_key_down, up_bp, lo_bp,farest_distance=30):
        uFexons_str = gene_positions[pos_key_up]["Fexons"]
        uBexons_str = gene_positions[pos_key_up]["Bexons"]
        lFexons_str = gene_positions[pos_key_down]["Fexons"]
        lBexons_str = gene_positions[pos_key_down]["Bexons"]
        result_up = find_nearest_breakpoint(up_bp, uBexons_str if self.upstream_gene.strand == '>' else uFexons_str, farest_distance)
        modify_upper_bp = result_up if result_up is not None else up_bp
        result_lo = find_nearest_breakpoint(lo_bp, lFexons_str if self.downstream_gene.strand == '>' else lBexons_str, farest_distance)
        modify_lower_bp = result_lo if result_lo is not None else lo_bp
        return modify_upper_bp, modify_lower_bp
    def validate(self, gene_positions, fusion_annot_dict, exon_boundary_distance=50):
        def get_smallest_pos(node):
            pos = 0
            if node.seq[0] == '>':
                pos = node.path[0].split(':')[-1].split('-')[0]
            elif node.seq[0] == '<':
                pos = node.path[-1].split(':')[-1].split('-')[0]
            return int(pos)
        def get_biggest_pos(node):
            pos = 0
            if node.seq[0] == '>':
                pos = node.path[-1].split(':')[-1].split('-')[-1]
            elif node.seq[0] == '<':
                pos = node.path[0].split(':')[-1].split('-')[-1]
            return int(pos)
        def is_near(a, b, max_distance):
            if a is None or b is None:
                return False
            return np.abs(int(a) - int(b)) <= max_distance
        def compute_sam_bp(Sam_aln,position):
            sam_bp = -1
            if position == 'upper':
                if Sam_aln.is_reverse:
                    sam_bp = Sam_aln.start
                else:
                    sam_bp = Sam_aln.end
            elif position == 'lower':
                if Sam_aln.is_reverse:
                    sam_bp = Sam_aln.end
                else:
                    sam_bp = Sam_aln.start
            return sam_bp
        def get_intergenic_protein_coding_genes(key1, key2, gene_positions):
            gene1,gene2 = gene_positions[key1], gene_positions[key2]
            if gene1["chrom"] != gene2["chrom"]:
                return []
            start = min(gene1["start"], gene2["start"])
            end = max(gene1["end"], gene2["end"])
            if gene1["strand"] == gene2["strand"]:
                return [
                g for g in gene_positions.values()
                if g["chrom"] == gene1["chrom"]
                and g["start"] > start and g["end"] < end
                and g["type"] == "protein_coding"
                and g["strand"] == gene1["strand"]
                and g["gene_name"] not in (gene1["gene_name"], gene2["gene_name"])
                ]
            else:
                return [
                    g for g in gene_positions.values()
                    if g["chrom"] == gene1["chrom"]
                    and g["start"] > start and g["end"] < end
                    and g["type"] == "protein_coding"
                    and g["gene_name"] not in (gene1["gene_name"], gene2["gene_name"])
                ]
        def check_node_segments_in_gene_overlap(node, expected_gene_name_str, aln_gene_name_str, node_aligned_chrom, gene_positions):
            key_expected_gene = f"{expected_gene_name_str}_{node_aligned_chrom}"
            key_aligned_gene = f"{aln_gene_name_str}_{node_aligned_chrom}"
            if key_expected_gene not in gene_positions or key_aligned_gene not in gene_positions:
                return False
            info_expected = gene_positions[key_expected_gene]
            info_aligned = gene_positions[key_aligned_gene]
            if not (info_expected["chrom"] == node_aligned_chrom and info_aligned["chrom"] == node_aligned_chrom):
                return False
            overlap_start = max(info_expected["start"], info_aligned["start"])
            overlap_end = min(info_expected["end"], info_aligned["end"])
            if overlap_start >= overlap_end:
                return False
            if not node.path:
                return False
            for segment_str in node.path:
                try:
                    parts = segment_str.split(':')
                    seg_chrom = parts[0]
                    seg_gene_name_in_path = parts[1]
                    seg_coords = parts[2].split('-')
                    seg_start = int(seg_coords[0])
                    seg_end = int(seg_coords[1])
                except (IndexError, ValueError) as e:
                    return False
                if seg_chrom != node_aligned_chrom or seg_gene_name_in_path != aln_gene_name_str:
                    return False
                if not (seg_start >= overlap_start and seg_end <= overlap_end):
                    return False
            return True
        def validate_position(gene_positions,upper_gene,lower_gene,upper_node,lower_node,distance=30):
            pos_key_up = f"{upper_gene.name}_{upper_gene.chrom}"
            pos_key_down = f"{lower_gene.name}_{lower_gene.chrom}"
            uFexons_str = gene_positions[pos_key_up]["Fexons"]
            uBexons_str = gene_positions[pos_key_up]["Bexons"]
            lFexons_str = gene_positions[pos_key_down]["Fexons"]
            lBexons_str = gene_positions[pos_key_down]["Bexons"]
            cupper_bp = int(compute_breakpoints_with_position_correct(upper_node,'upper'))
            clower_bp = int(compute_breakpoints_with_position_correct(lower_node,'lower'))
            result_up = find_nearest_breakpoint(cupper_bp, uBexons_str if upper_gene.strand == '>' else uFexons_str, distance)
            result_lo = find_nearest_breakpoint(clower_bp, lFexons_str if lower_gene.strand == '>' else lBexons_str, distance)
            if result_up is None or result_lo is None:
                return False
            else:
                return result_up,result_lo
        self.confidence_level = 0
        if len(self.minigraph_upper_alignments) < 1 or len(self.minigraph_lower_alignments) < 1:
            self.confidence_level = 3
            return False
        minimap2_check_passed = False
        minigraph_check_passed = True
        minimap2_upstream_hit, minimap2_downstream_hit = None, None
        high_mapq_count = 0
        if not self.minimap2_full_alignments:
            self.confidence_level = 3
            return True 
        candidate_hits = [aln for aln in self.minimap2_full_alignments if not aln.is_secondary and aln.mapq >= 1]
        upper_node, lower_node = self.minigraph_upper_alignments[0], self.minigraph_lower_alignments[0]
        upper_node_sp, lower_node_sp = get_smallest_pos(upper_node), get_smallest_pos(lower_node)
        upper_node_ep, lower_node_ep = get_biggest_pos(upper_node), get_biggest_pos(lower_node)
        up_minimap2_aln, lo_minimap2_aln = None, None
        for aln in candidate_hits:
            matched_this_aln = False
            if aln.ref_name == self.upstream_gene.chrom:
                if is_near(aln.start,upper_node_sp,500) or is_near(aln.end,upper_node_ep,500):
                    up_minimap2_aln = aln
                    minimap2_upstream_hit = aln
                    matched_this_aln = True
            if aln.ref_name == self.downstream_gene.chrom:
                if is_near(aln.start,lower_node_sp,500) or is_near(aln.end,lower_node_ep,500):
                    lo_minimap2_aln = aln
                    minimap2_downstream_hit = aln
                    matched_this_aln = True
            if matched_this_aln:
                high_mapq_count += 1
        if minimap2_upstream_hit and minimap2_downstream_hit and high_mapq_count >= 2:
            minimap2_check_passed = True
        if self.upstream_gene.chrom != upper_node.path[0].split(':')[0] or self.downstream_gene.chrom != lower_node.path[0].split(':')[0] or self.upstream_gene.strand != upper_node.seq[0] or self.downstream_gene.strand != lower_node.seq[0]:
            minigraph_check_passed = False
        upstream_gene_name = self.upstream_gene.name
        downstream_gene_name = self.downstream_gene.name
        pos_key_up = f"{upstream_gene_name}_{self.upstream_gene.chrom}"
        pos_key_down = f"{downstream_gene_name}_{self.downstream_gene.chrom}"
        aln_up_gene_name = upper_node.path[0].split(':')[1]
        aln_down_gene_name = lower_node.path[0].split(':')[1]
        up_gene_match = (upstream_gene_name == aln_up_gene_name) or (upstream_gene_name in aln_up_gene_name) or (aln_up_gene_name in upstream_gene_name)
        down_gene_match = (downstream_gene_name == aln_down_gene_name) or (downstream_gene_name in aln_down_gene_name) or (aln_down_gene_name in downstream_gene_name)
        if not up_gene_match:
            if self.upstream_gene.chrom == upper_node.path[0].split(':')[0]:
                if check_node_segments_in_gene_overlap(upper_node, upstream_gene_name, aln_up_gene_name, self.upstream_gene.chrom, gene_positions):
                    up_gene_match = True
        if not down_gene_match:
            if self.downstream_gene.chrom == lower_node.path[0].split(':')[0]:
                if check_node_segments_in_gene_overlap(lower_node, downstream_gene_name, aln_down_gene_name, self.downstream_gene.chrom, gene_positions):
                    down_gene_match = True
        if not (up_gene_match and down_gene_match):
            minigraph_check_passed = False

        final_upstream_bp, final_downstream_bp = self.upstream_gene.breakpoint, self.downstream_gene.breakpoint
        if minimap2_check_passed:
            up_minimap_bp = compute_sam_bp(up_minimap2_aln, 'upper')
            lo_minimap_bp = compute_sam_bp(lo_minimap2_aln, 'lower')
            mum_bp,mlm_bp = self.get_modify_bp(gene_positions, pos_key_up, pos_key_down, up_minimap_bp, lo_minimap_bp, 100)
        if minigraph_check_passed:
            position_check_result = validate_position(gene_positions,self.upstream_gene,self.downstream_gene,upper_node,lower_node)
            if position_check_result == False:
                minigraph_check_passed = False
            else:
                modify_upper_bp, modify_lower_bp = position_check_result
        if minigraph_check_passed:
            final_upstream_bp, final_downstream_bp = modify_upper_bp, modify_lower_bp
            check = check_fusion_pair(self.upstream_gene.name, self.downstream_gene.name,fusion_annot_dict)
            if check == "NORMAL_FUSION":
                self.confidence_level = 2
                logging.info(f"{self.get_id()}: The confidence level is set to 2 because it is a normal fusion verified by the database.")
            else:
                chrom_is_same = self.upstream_gene.chrom == self.downstream_gene.chrom
                upg_real_strand = gene_positions[pos_key_up]["strand"]
                downg_real_strand = gene_positions[pos_key_down]["strand"]
                real_strand_is_same = upg_real_strand == downg_real_strand
                need_validate = False
                if minimap2_check_passed:
                    if is_near(final_upstream_bp, up_minimap_bp, exon_boundary_distance) and is_near(final_downstream_bp, lo_minimap_bp, exon_boundary_distance):
                        self.confidence_level = 5
                        logging.info(f"{self.get_id()}: Confidence is set to 5, reason: minigraph_check_passed AND minimap2_check_passed. distance<={exon_boundary_distance}, breakpoint positions match.")
                    else:
                        self.confidence_level = 4
                        logging.info(f"{self.get_id()}: Confidence set to 4, reason: minigraph_check_passed AND minimap2_check_passed, but breakpoint positions do not match exactly.")
                else:
                    if chrom_is_same:
                        need_validate = True
                    else:
                        logging.info(f"{self.get_id()}: Confidence is set to 3, reason: minigraph_check_passed but minimap2_check not passed.")
                        self.confidence_level = 3
                if need_validate:
                    bp_distance = abs(final_upstream_bp - final_downstream_bp)
                    if bp_distance <= 100000:
                        intergenic_pc_genes = get_intergenic_protein_coding_genes(pos_key_up, pos_key_down, gene_positions)
                        found_sstrand_pc_gene = False
                        if intergenic_pc_genes:
                            for pc_gene_info in intergenic_pc_genes:
                                gene_name_part = pc_gene_info["gene_name"].split('_chr')[0]
                                if 'ENSG' not in pc_gene_info["gene_name"] and not all(gene in gene_name_part for gene in [self.upstream_gene.name, self.downstream_gene.name]):
                                    found_sstrand_pc_gene = True
                                    break
                            if found_sstrand_pc_gene:
                                self.confidence_level = 5
                                logging.info(f"{self.get_id()}: The confidence of the gene in the same direction is set to 5, because there is a pc gene (not ENSG)")
                            else:
                                if bp_distance > 50000:
                                    self.confidence_level = 4
                                    logging.info(f"{self.get_id()}: The confidence of the gene in the same direction is set to 4, there is no pc gene in the same direction, there is a pc gene in the opposite direction, but the breakpoint distance is greater than 50kbp.")
                                else:
                                    self.confidence_level = 2
                                    logging.info(f"{self.get_id()}: The confidence of the gene in the same direction is set to 2. The reason is: there is no pc gene in the same direction, but there are pcs in other directions, and the breakpoint distance is less than 50kbp.")
                        else:
                            self.confidence_level = 2
                            logging.info(f"{self.get_id()}: The confidence level is set to 2, because there is no pc gene between genes in the same direction on the same chromosome.")
                    else:
                        self.confidence_level = 5
                        logging.info(f"{self.get_id()}:The confidence level is set to 5, because the breakpoint distance between two genes on the same chromosome is greater than 100kbp.")

        else:
            if minimap2_check_passed:
                final_upstream_bp = mum_bp
                final_downstream_bp = mlm_bp
                self.confidence_level = 3
                logging.info(f"{self.get_id()}: Confidence is set to 3. minimap2_check_passed but minigraph_check not passed.")
            else:
                self.confidence_level = 1
        if self.confidence_level not in [1,2]:
            if gene_positions[pos_key_up]["type"] != 'protein_coding' or gene_positions[pos_key_down]["type"] != 'protein_coding':
                self.confidence_level = max(3, self.confidence_level-1)
        self.upstream_gene.breakpoint = final_upstream_bp
        self.downstream_gene.breakpoint = final_downstream_bp
        return self.confidence_level > 1
    def format_output(self):
        g1 = self.upstream_gene
        g2 = self.downstream_gene
        confidence_desc = {
            -1: "ERROR Fusion",
            0: "ERROR Fusion",
            1: "Impossible Fusion",
            2: "Potential Artifacts / Read-through",
            3: "Low Confidence Fusion",
            4: "Medium Confidence Fusion",
            5: "High Confidence Fusion"
        }
        fields = [
            g1.name,                                                  
            g2.name,                                                  
            f"{g1.chrom}:{g1.breakpoint}; {g2.chrom}:{g2.breakpoint}",
            f"{g1.chrom}{g1.strand}: {g1.start}, {g1.end}",           
            f"{g2.chrom}{g2.strand}: {g2.start}, {g2.end}",           
            str(self.get_read_count()),                               
            f"flag = {self.confidence_level}: {confidence_desc[self.confidence_level]}",
            ";".join(self.supporting_reads)
        ]
        return "\t".join(fields)

def parse_cigar_for_query_details(cigar_string):
    if cigar_string == '*':
        return 0, 0, 0, 0, 0, 0, 0, 0, 0
    all_ops_match = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    if not all_ops_match:
        return 0, 0, 0, 0, 0, 0, 0, 0, 0
    all_ops = [(int(length), op) for length, op in all_ops_match]
    h_start_len = 0
    h_end_len = 0
    if all_ops[0][1] == 'H':
        h_start_len = all_ops[0][0]
    if len(all_ops) > 1 and all_ops[-1][1] == 'H':
        h_end_len = all_ops[-1][0]
    elif len(all_ops) == 1 and all_ops[0][1] == 'H':
        h_end_len = 0 
    seq_ops_start_idx = 1 if h_start_len > 0 else 0
    seq_ops_end_idx = len(all_ops) - 1 if h_end_len > 0 else len(all_ops)
    if seq_ops_start_idx >= seq_ops_end_idx:
         return 0, 0, 0, 0, 0, 0, 0, h_start_len, h_end_len
    seq_ops = all_ops[seq_ops_start_idx:seq_ops_end_idx]
    q_len_cigar = 0
    r_len_cigar = 0
    q_core_len = 0
    s_start_len = 0
    s_end_len = 0
    if not seq_ops:
        return 0, 0, 0, 0, 0, 0, 0, h_start_len, h_end_len

    if seq_ops[0][1] == 'S':
        s_start_len = seq_ops[0][0]

    if len(seq_ops) > 1 and seq_ops[-1][1] == 'S':
        s_end_len = seq_ops[-1][0]
    elif len(seq_ops) == 1 and seq_ops[0][1] == 'S':
        s_start_len = seq_ops[0][0]
        s_end_len = 0
    current_q_pos_in_seq_segment = 0
    temp_q_aln_start = -1
    for length, op_char in seq_ops:
        if op_char in ('M', 'I', 'S', '=', 'X'):
            q_len_cigar += length
            if op_char in ('M', 'I', '=', 'X'):
                q_core_len += length
                if temp_q_aln_start == -1: 
                    temp_q_aln_start = current_q_pos_in_seq_segment
            current_q_pos_in_seq_segment += length

        if op_char in ('M', 'D', 'N', '=', 'X'):
            r_len_cigar += length

    q_aln_start = s_start_len
    q_aln_end = q_len_cigar - s_end_len
    if q_core_len == 0:
        q_aln_start = s_start_len
        q_aln_end = s_start_len
    return q_len_cigar, r_len_cigar, q_aln_start, q_aln_end, s_start_len, s_end_len, q_core_len, h_start_len, h_end_len

def run_command(command, stdin=None):
    process = Popen(shlex.split(command), stdin=stdin, stdout=PIPE, stderr=PIPE,
                    bufsize=8388608, universal_newlines=True)
    return process

class FusionProcessor:
    def __init__(self, gene_positions, reference_genome_fasta, reference_graph_gfa, min_supporting_reads=2, threads=8, max_breakpoint_distance=200, temp_dir=None):
        self.gene_positions = gene_positions
        self.threads = threads
        self.min_supporting_reads = min_supporting_reads
        self.reference_genome_fasta = reference_genome_fasta
        self.reference_graph_gfa = reference_graph_gfa
        self.max_breakpoint_distance = max_breakpoint_distance
        self.temp_dir = temp_dir
        if self.temp_dir and not os.path.exists(self.temp_dir):
            os.makedirs(self.temp_dir)
        self.alignments_upper = None
        self.alignments_lower = None
        self.totalreads = 0
        self.count = 0
        self.minigraph_split_alignments_upper = defaultdict(list)
        self.minigraph_split_alignments_lower = defaultdict(list)
        self.minimap2_full_alignments_sam = defaultdict(list)
    
    def _create_fusion_candidates_from_data(self, final_fusion):
        fusion_candidates = []
        for fusion_id, fusion_data in final_fusion.items():
            match = re.match(r"(.+)_bpGroup_(\d+)_pathGroup_(\d+)", fusion_id)
            if not match:
                continue
            base_fusion_id, bp_group_id, path_group_id = match.groups()
            gene_names = base_fusion_id.split(':')
            if len(gene_names) != 2:
                continue
            upstream_name, downstream_name = gene_names
            upstream_info = None
            downstream_info = None
            supporting_reads = set()
            read_alignment_info = {}
            for gene_key, gene_info in fusion_data.items():
                if upstream_name in gene_key:
                    upstream_info = gene_info
                    if "ReadAlignmentInfo" in gene_info and gene_info["ReadAlignmentInfo"]:
                        read_alignment_info = gene_info["ReadAlignmentInfo"]
                elif downstream_name in gene_key:
                    downstream_info = gene_info
                if "Reads" in gene_info and isinstance(gene_info["Reads"], list):
                    supporting_reads.update(gene_info["Reads"])
            if not upstream_info or not downstream_info:
                continue
            upstream_chrom = upstream_info.get('Chrom', '')
            upstream_bp = int(upstream_info.get('Breakpoint', 0))
            upstream_strand = upstream_info.get('Path', [''])[0]
            downstream_chrom = downstream_info.get('Chrom', '')
            downstream_bp = int(downstream_info.get('Breakpoint', 0))
            downstream_strand = downstream_info.get('Path', [''])[0]
            pos_key_up = f"{upstream_name}_{upstream_chrom}"
            pos_key_down = f"{downstream_name}_{downstream_chrom}"
            upstream_start = self.gene_positions.get(pos_key_up, {}).get('start', 0)
            upstream_end = self.gene_positions.get(pos_key_up, {}).get('end', 0)
            downstream_start = self.gene_positions.get(pos_key_down, {}).get('start', 0)
            downstream_end = self.gene_positions.get(pos_key_down, {}).get('end', 0)
            upstream_gene = Gene(
                name=upstream_name,
                chrom=upstream_chrom,
                breakpoint=upstream_bp,
                strand=upstream_strand,
                start=upstream_start,
                end=upstream_end
            )
            downstream_gene = Gene(
                name=downstream_name,
                chrom=downstream_chrom,
                breakpoint=downstream_bp,
                strand=downstream_strand,
                start=downstream_start,
                end=downstream_end
            )
            fusion_candidate = FusionCandidate(
                upstream_gene=upstream_gene,
                downstream_gene=downstream_gene,
                read_alignment_info=read_alignment_info,
                supporting_reads=list(supporting_reads),
                group_id=bp_group_id,
                path_group_id=path_group_id
            )
            fusion_candidates.append(fusion_candidate)

        return fusion_candidates
    
    def _process_fusion_chunk(self, chunk, sequences, output_fasta_split, output_fasta_full, output_lock, min_cseq_length):
        local_successful_count = 0
        local_processed_count = 0
        for fusion in chunk:
            fusion_id = fusion.get_id()
            read_ids = fusion.supporting_reads
            if len(read_ids) < 1:
                continue
            logging.debug(f"Handle fusion {fusion_id} (supports reads: {len(read_ids)})")
            local_processed_count += 1
            upper_sequences_for_msa = []
            lower_sequences_for_msa = []
            valid_read_count_for_msa = 0
            for read_id in read_ids:
                if read_id not in sequences: continue
                try:
                    full_seq = str(sequences[read_id])
                except KeyError: continue
                if read_id in fusion.read_alignment_info:
                    upper_info = fusion.read_alignment_info[read_id].get('upper', {})
                    lower_info = fusion.read_alignment_info[read_id].get('lower', {})
                    upper_start = upper_info.get('start', -1)
                    upper_end = upper_info.get('end', -1)
                    lower_start = lower_info.get('start', -1)
                    lower_end = lower_info.get('end', -1)
                    upper_seq_full_part = full_seq[:upper_end]
                    lower_seq_full_part = full_seq[lower_start:]
                    if upper_seq_full_part and lower_seq_full_part and \
                       upper_info.get('strand') in ['>', '<'] and \
                       lower_info.get('strand') in ['>', '<']:
                        upper_sequences_for_msa.append(upper_seq_full_part)
                        lower_sequences_for_msa.append(lower_seq_full_part)
                        valid_read_count_for_msa += 1
                    else:
                         logging.warning(f"Warning: Read '{read_id}' (for {fusion_id}) has invalid upstream/downstream sequence or strand information extracted, skipping.")
                else:
                     logging.warning(f"Warning: Read '{read_id}' (for {fusion_id}) is missing alignment information.")
            if valid_read_count_for_msa < 1:
                continue
            try:
                aligner = pa.msa_aligner()
                upper_result = aligner.msa(upper_sequences_for_msa, out_cons=True, out_msa=False)
                upper_consensus = upper_result.cons_seq[0] if upper_result.cons_seq else ""
                lower_result = aligner.msa(lower_sequences_for_msa, out_cons=True, out_msa=False)
                lower_consensus = lower_result.cons_seq[0] if lower_result.cons_seq else ""
                if upper_consensus and lower_consensus:
                    if len(upper_consensus) < min_cseq_length or len(lower_consensus) < min_cseq_length:
                         logging.info(f"Skipping {fusion_id}: consensus sequence too short ({len(upper_consensus)}, {len(lower_consensus)} < {min_cseq_length})")
                         continue
                    full_consensus = upper_consensus + lower_consensus
                    num_reads_in_msa = valid_read_count_for_msa
                    upper_header = f">{fusion_id}_{num_reads_in_msa}_upper"
                    upper_output = f"{upper_header}\n" + upper_consensus + "\n"
                    lower_header = f">{fusion_id}_{num_reads_in_msa}_lower"
                    lower_output = f"{lower_header}\n" + lower_consensus + "\n"
                    full_header = f">{fusion_id}_{num_reads_in_msa}_full"
                    full_output = f"{full_header}\n" + full_consensus + "\n"
                    with output_lock:
                        output_fasta_split.write(upper_output)
                        output_fasta_split.write(lower_output)
                        output_fasta_full.write(full_output)
                    try:
                        fusion.add_upper_consensus_sequence(upper_consensus)
                        fusion.add_lower_consensus_sequence(lower_consensus)
                        fusion.consensus_sequence = full_consensus
                    except AttributeError:
                        fusion.upper_consensus_sequence = upper_consensus
                        fusion.lower_consensus_sequence = lower_consensus
                        fusion.consensus_sequence = full_consensus
                    local_successful_count += 1
                    logging.debug(f"Generate consensus sequence: {fusion_id}, based on {num_reads_in_msa} sequences")
                else:
                    logging.warning(f"WARNING: {fusion_id} failed to generate upstream or downstream consensus sequence")

            except Exception as e:
                logging.exception(f"Error: Exception occurred while processing {fusion_id} (MSA input sequence number: {valid_read_count_for_msa}): {str(e)}")

        return local_successful_count, local_processed_count
    
    def generate_consensus_sequences(self, fusion_candidates, fasta_sequences, output_fasta_split, output_fasta_full, min_cseq_length=30, thread=4):
        start_time = time.time()
        print(f"Starting to generate fusion gene consensus sequence (using {thread} thread)...")
        sequences = fasta_sequences
        self.totalreads = len(sequences)
        output_lock = threading.Lock()
        total_successful_count = 0
        total_processed_count = 0
        num_candidates = len(fusion_candidates)
        if num_candidates == 0:
            for out_file in [output_fasta_split, output_fasta_full]:
                try:
                    with open(out_file, 'w') as f: pass
                except IOError as e:
                    logging.error(f"Unable to create empty output file '{out_file}':{e}")
            return 0

        actual_threads = min(thread, num_candidates)
        chunk_size = math.ceil(num_candidates / actual_threads) if actual_threads > 0 else 1
        chunks = [fusion_candidates[i:i + chunk_size] for i in range(0, num_candidates, chunk_size)]
        print(f"Split {num_candidates} candidates into {len(chunks)} chunks, using {actual_threads} threads")
        try:
            with open(output_fasta_split, 'w') as f_split, open(output_fasta_full, 'w') as f_full:
                with ThreadPoolExecutor(max_workers=actual_threads) as executor:
                    future_to_chunk_index = {
                        executor.submit(self._process_fusion_chunk, chunk, sequences, f_split, f_full, output_lock, min_cseq_length): i
                        for i, chunk in enumerate(chunks)
                    }
                    for future in as_completed(future_to_chunk_index):
                        chunk_index = future_to_chunk_index[future]
                        try:
                            result = future.result()
                            local_successful, local_processed = result
                            total_successful_count += local_successful
                            total_processed_count += local_processed
                        except Exception as exc:
                            logging.error(f"Exception occurred while executing chunk {chunk_index}: {exc}", exc_info=True)
        except IOError as e:
            logging.error(f"ERROR: Cannot open or write output file '{output_fasta_split}' or '{output_fasta_full}': {e}")
            return total_successful_count
        except Exception as e:
             logging.exception(f"An unexpected error occurred during multithreading: {e}")
             return total_successful_count
        self.count = total_processed_count
        elapsed_time = time.time() - start_time
        print("-" * 30)
        print(f"Consensus sequence generation completed!")
        print(f"Time taken: {elapsed_time:.2f} seconds")
        print(f"Split sequence output files: {output_fasta_split}")
        print(f"Full sequence output file: {output_fasta_full}")
        print("-" * 30)
        return total_successful_count
    
    def align_full_consensus_minimap2(self, consensus_fasta_full, sam_output_file, threads=32):
        command = f'minimap2 -ax splice -m20 -s20 --secondary=no -t {threads} {self.reference_genome_fasta} {consensus_fasta_full}'
        minimap2_process = run_command(command)
        alignments_by_fusion_sam = defaultdict(list)
        try:
            with open(sam_output_file, 'w') as samfile:
                for line in minimap2_process.stdout:
                    samfile.write(line)
                    parsed_sam = parse_sam_line(line)
                    if parsed_sam:
                        query_parts = parsed_sam.query_name.split('_')
                        if len(query_parts) >= 5 and query_parts[-1] == 'full':
                             fusion_id = f"{query_parts[0]}_{query_parts[1]}_{query_parts[2]}_{query_parts[3]}"
                             alignments_by_fusion_sam[fusion_id].append(parsed_sam)

            stderr_output = minimap2_process.stderr.read()
            minimap2_process.wait()
            if minimap2_process.returncode != 0:
                logging.error(f"Minimap2 comparison failed, return code: {minimap2_process.returncode}")
                logging.error(f"Minimap2 stderr:\n{stderr_output}")
            else:
                print(f"Minimap2 alignments completed, parsing {len(alignments_by_fusion_sam)} fused SAM alignments")
                if stderr_output:
                    logging.debug(f"Minimap2 stderr:\n{stderr_output}")

        except Exception as e:
            logging.exception(f"Error processing Minimap2 output or writing SAM file: {e}")
            return alignments_by_fusion_sam

        return alignments_by_fusion_sam
    
    def align_split_consensus_minigraph(self, consensus_fasta_split, gaf_output_file, threads=32):
        command = f'minigraph -cxlr -k12 -w1 -m30,20 -t{threads} {self.reference_graph_gfa} {consensus_fasta_split}'
        minigraph_process = run_command(command)
        alignments_by_fusion_upper_gaf = defaultdict(list)
        alignments_by_fusion_lower_gaf = defaultdict(list)
        try:
            with open(gaf_output_file, 'w') as gaffile:
                for line in minigraph_process.stdout:
                    gaffile.write(line)
                    parsed_gaf = parse_readinfo(line.strip())
                    if parsed_gaf:
                        read_id = parsed_gaf.read_id
                        if read_id.endswith("_upper") or read_id.endswith("_lower"):
                            query_parts = read_id.split('_')
                            if len(query_parts) >= 5:
                                fusion_id = f"{query_parts[0]}_{query_parts[1]}_{query_parts[2]}_{query_parts[3]}"
                                if read_id.endswith("_upper"):
                                    alignments_by_fusion_upper_gaf[fusion_id].append(parsed_gaf)
                                else:
                                    alignments_by_fusion_lower_gaf[fusion_id].append(parsed_gaf)

            stderr_output = minigraph_process.stderr.read()
            minigraph_process.wait()
            if minigraph_process.returncode != 0:
                logging.error(f"Minigraph comparison failed, return code: {minigraph_process.returncode}")
                logging.error(f"Minigraph stderr:\n{stderr_output}")
            else:
                print(f"Minigraph alignment completed, a total of {len(alignments_by_fusion_upper_gaf) + len(alignments_by_fusion_lower_gaf)} upstream/downstream gaf alignment results parsed")
                if stderr_output:
                    logging.debug(f"Minigraph stderr:\n{stderr_output}")

        except Exception as e:
            logging.exception(f"Error processing Minigraph output or writing gaf file: {e}")
            return alignments_by_fusion_upper_gaf, alignments_by_fusion_lower_gaf

        return alignments_by_fusion_upper_gaf, alignments_by_fusion_lower_gaf

    def validate_fusion_candidates(self, fusion_candidates,fusion_annot_dict):
        validated_fusions = []
        bp_group_candidates = defaultdict(list)
        for fusion in fusion_candidates:
            bp_group_id = fusion.get_bp_group_id()
            bp_group_candidates[bp_group_id].append(fusion)
        print(f"Starting validation of {len(fusion_candidates)} candidates (from {len(bp_group_candidates)} breakpoint groups)...")
        for bp_group_id, group_fusions in bp_group_candidates.items():
            valid_path_groups_in_bp_group = []
            for fusion in group_fusions:
                fusion_id = fusion.get_id()
                fusion.minigraph_upper_alignments = []
                fusion.minigraph_lower_alignments = []
                fusion.minimap2_full_alignments = []
                if fusion_id in self.minigraph_split_alignments_upper:
                    for alignment in self.minigraph_split_alignments_upper[fusion_id]:
                        fusion.add_minigraph_upper_alignment(alignment)
                if fusion_id in self.minigraph_split_alignments_lower:
                    for alignment in self.minigraph_split_alignments_lower[fusion_id]:
                        fusion.add_minigraph_lower_alignment(alignment)
                if fusion_id in self.minimap2_full_alignments_sam:
                    for alignment in self.minimap2_full_alignments_sam[fusion_id]:
                        fusion.add_minimap2_full_alignment(alignment)
                if fusion.validate(self.gene_positions,fusion_annot_dict):
                    valid_path_groups_in_bp_group.append(fusion)

            validated_fusions.extend(valid_path_groups_in_bp_group)
        print(f"Validation completed: {len(validated_fusions)} fusion events passed validation (confidence > 1)")
        return validated_fusions

    def align_process(self,consensus_fasta_split,consensus_fasta_full,gaf_output_minigraph,sam_output_minimap2):
        align_threads = max(1, self.threads // 2)
        print(f"Allocate {align_threads} threads to Minigraph and Minimap2 respectively")
        future_minigraph = None
        future_minimap2 = None
        results_minigraph = (defaultdict(list), defaultdict(list)) # (upper, lower)
        results_minimap2 = defaultdict(list)
        with ThreadPoolExecutor(max_workers=2) as executor:
            future_minigraph = executor.submit(
                self.align_split_consensus_minigraph,
                consensus_fasta_split, gaf_output_minigraph, align_threads
            )

            future_minimap2 = executor.submit(
                self.align_full_consensus_minimap2,
                consensus_fasta_full, sam_output_minimap2, align_threads
            )
            try:
                results_minigraph = future_minigraph.result()
                self.minigraph_split_alignments_upper, self.minigraph_split_alignments_lower = results_minigraph
            except Exception as e:
                logging.exception("Minigraph comparison task failed!")
            try:
                results_minimap2 = future_minimap2.result()
                self.minimap2_full_alignments_sam = results_minimap2
            except Exception as e:
                logging.exception("Minimap2 comparison task failed!")
    
    def post_process_fusion_results(self, all_fusion_candidates, validated_fusions, include_isoforms=False):
        gene_pair_to_fusions = defaultdict(list)
        for fusion in all_fusion_candidates:
            gene_pair = tuple(sorted([
                (fusion.upstream_gene.name, fusion.upstream_gene.chrom),
                (fusion.downstream_gene.name, fusion.downstream_gene.chrom)
            ]))
            gene_pair_to_fusions[gene_pair].append(fusion)
        final_fusions = []
        processed_fusions_overall = set()
        for gene_pair, fusions in gene_pair_to_fusions.items():
            validated = [f for f in fusions if f in validated_fusions]
            reliable_fusions = [f for f in validated if f.confidence_level == 5]
            suspected_fusions = [f for f in validated if f.confidence_level == 4]
            lowconfi_fusions = [f for f in validated if f.confidence_level == 3]
            verylowconfi_fusions = [f for f in validated if f.confidence_level == 2]
            validated = reliable_fusions + suspected_fusions + lowconfi_fusions + verylowconfi_fusions
            if reliable_fusions:
                if include_isoforms:
                    processed_in_group = set()
                    reliable_fusions.sort(key=lambda f: f.get_read_count(), reverse=True)
                    for reliable_fusion in reliable_fusions:
                        if reliable_fusion in processed_in_group or reliable_fusion in processed_fusions_overall:
                            continue
                        current_support_reads_set = set(reliable_fusion.supporting_reads)
                        absorbed_in_this_round = set()
                        for other_fusion in validated:
                            if (other_fusion == reliable_fusion or
                                other_fusion in processed_in_group or
                                other_fusion in processed_fusions_overall):
                                continue
                            if (other_fusion.upstream_gene.breakpoint == reliable_fusion.upstream_gene.breakpoint and
                                other_fusion.downstream_gene.breakpoint == reliable_fusion.downstream_gene.breakpoint):
                                current_support_reads_set.update(other_fusion.supporting_reads)
                                absorbed_in_this_round.add(other_fusion)
                        reliable_fusion.supporting_reads = list(current_support_reads_set)
                        final_fusions.append(reliable_fusion)
                        processed_in_group.add(reliable_fusion)
                        processed_in_group.update(absorbed_in_this_round)
                        processed_fusions_overall.add(reliable_fusion)
                        processed_fusions_overall.update(absorbed_in_this_round)
                else:
                    breakpoint_groups = defaultdict(list)
                    for fusion in reliable_fusions:
                        bp_key = (fusion.upstream_gene.breakpoint, fusion.downstream_gene.breakpoint)
                        breakpoint_groups[bp_key].append(fusion)
                    best_bp_group = None
                    max_reads = -1
                    for bp_key, group in breakpoint_groups.items():
                        all_reads = set()
                        for fusion in group:
                            all_reads.update(fusion.supporting_reads)

                        if len(all_reads) > max_reads:
                            max_reads = len(all_reads)
                            best_bp_group = group
                    if best_bp_group:
                        base_fusion = best_bp_group[0]
                        all_reads = set(base_fusion.supporting_reads)
                        for fusion in best_bp_group[1:]:
                            all_reads.update(fusion.supporting_reads)
                        for fusion in validated:
                            if fusion not in best_bp_group:
                                all_reads.update(fusion.supporting_reads)
                        base_fusion.supporting_reads = list(all_reads)
                        final_fusions.append(base_fusion)
                        processed_fusions_overall.update(validated)

            elif any([suspected_fusions, lowconfi_fusions, verylowconfi_fusions]):
                highest_confi_fusions = None
                if suspected_fusions:
                    highest_confi_fusions = suspected_fusions
                elif lowconfi_fusions:
                    highest_confi_fusions = lowconfi_fusions
                else:
                    highest_confi_fusions = verylowconfi_fusions
                best_fusion = max(highest_confi_fusions, key=lambda f: f.get_read_count())
                all_reads = set(best_fusion.supporting_reads)
                for other_fusion in validated:
                    if other_fusion != best_fusion:
                        all_reads.update(other_fusion.supporting_reads)
                best_fusion.supporting_reads = list(all_reads)
                fusion_gene_name = [best_fusion.upstream_gene.name, best_fusion.downstream_gene.name]
                final_fusions.append(best_fusion)
                processed_fusions_overall.update(validated)
        for fusion in final_fusions:
            if fusion.confidence_level == 2:
                continue
            else:
                fusion_gene_name = [fusion.upstream_gene.name, fusion.downstream_gene.name]
                if 'ENSG' in fusion_gene_name[0] or 'ENSG' in fusion_gene_name[1]:
                    fusion.confidence_level = max(3, fusion.confidence_level - 1)
        return final_fusions

    def write_results(self, fusions, output_file):
        with open(output_file, 'w') as f:
            header = "#gene1 name\tgene2 name\tbreak points\tgene1 position\tgene2 position\tsupport num\trank class\tsupporting reads information\n"
            f.write(header)
            for fusion in fusions:
                if fusion.get_read_count() >= self.min_supporting_reads:
                    f.write(fusion.format_output() + "\n")

        print(f"Results written to: {output_file}")
        
    def process_fusion_data(self, final_fusion, fasta_sequences, fusion_annot_dict, output_dir):
        consensus_fasta_split = os.path.join(self.temp_dir, "consensus_split.fasta")
        consensus_fasta_full = os.path.join(self.temp_dir, "consensus_full.fasta")
        gaf_output_minigraph = os.path.join(self.temp_dir, "alignment_split_minigraph.gaf")
        sam_output_minimap2 = os.path.join(self.temp_dir, "alignment_full_minimap2.sam")
        final_output_csv = os.path.join(output_dir, "final_validated_fusions.csv")
        validated_output_csv = os.path.join(self.temp_dir, "validated_fusions_raw.csv") 

        print("Creating ensemble candidates from preliminary clustering results...")
        fusion_candidates = self._create_fusion_candidates_from_data(final_fusion)
        print(f"{len(fusion_candidates)} fusion candidates created")
        if not fusion_candidates:
            logging.warning("There are no valid fusion candidates and the process terminates.")
            return [], [], []

        print("Generate consensus sequences (split and full)...")
        successful_consensus = self.generate_consensus_sequences(fusion_candidates, fasta_sequences, consensus_fasta_split, consensus_fasta_full,min_cseq_length=min_cseq_length, thread=self.threads)
        if successful_consensus == 0:
            logging.warning("If no consensus sequence is successfully generated, the process terminates.")
            for f in [consensus_fasta_split, consensus_fasta_full]:
                if os.path.exists(f) and os.path.getsize(f) == 0: os.remove(f)
            return fusion_candidates, [], []
        
        print("Parallel execution of Minigraph (split) and Minimap2 (full) alignments...")
        self.align_process(consensus_fasta_split,consensus_fasta_full,gaf_output_minigraph,sam_output_minimap2)

        print("Verify fusion candidates...")
        validated_fusions = self.validate_fusion_candidates(fusion_candidates,fusion_annot_dict)
        # self.write_results(validated_fusions, validated_output_csv)

        print("Post-processing and final screening...")
        final_validate_fusions = self.post_process_fusion_results(fusion_candidates, validated_fusions, include_isoforms=False)
        print(f"Processing is complete, and {len(final_validate_fusions)} final fusion events are screened out.")

        self.write_results(final_validate_fusions, final_output_csv)
        # return fusion_candidates, validated_fusions, final_validate_fusions

def align_and_process_reads(fasta_path: str, gfa_path: str, output_gaf_path: str, threads: int = 8, k: int = 12, w: int = 1, m: str = '30,20'):
    print(f"Starting alignment for {fasta_path}, output will be saved to {output_gaf_path}")
    command = f'minigraph -cxlr -k{k} -w{w} -m{m} -t{threads} {gfa_path} {fasta_path}'
    print(f"Executing Minigraph command: {command}")
    
    minigraph_process = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE, text=True)
    fusion_cluster = defaultdict(list)
    current_read_id = None
    current_group = []

    print("Reading GAF stream, writing to file, and processing on-the-fly...")
    try:
        with open(output_gaf_path, 'w') as gaf_file:
            for line in minigraph_process.stdout:
                gaf_file.write(line)
                try:
                    info = parse_readinfo(line.strip())
                    read_id = info.read_id
                    if info.path:
                        gene_name = info.path[0].split(':')[1]
                    if current_read_id is not None and read_id != current_read_id:
                        if len(current_group) > 1:
                            fusion_cluster[current_read_id].extend(current_group)
                        current_group = []
                    current_read_id = read_id
                    current_group.append(info)
                except Exception as e:
                    logging.warning(f"Could not parse GAF line: '{line.strip()}'. Error: {e}")
        if current_group and len(current_group) > 1:
            fusion_cluster[current_read_id].extend(current_group)
    except IOError as e:
        logging.error(f"Could not write to GAF file {output_gaf_path}: {e}")
    
    stderr_output = minigraph_process.stderr.read()
    minigraph_process.wait()
    if minigraph_process.returncode != 0:
        logging.error(f"Minigraph process failed with return code {minigraph_process.returncode}")
        logging.error(f"Minigraph stderr:\n{stderr_output}")
    else:
        print("Minigraph process completed successfully.")
    print(f"Identified {len(fusion_cluster)} potential fusion reads from alignment.")
    return fusion_cluster

def _load_reference_data(fusion_annot_lib_path: str, fasta_file_path: str):
    print(f"Starting to load fusion annotation library: {fusion_annot_lib_path}")
    fusion_annot_dict = parse_fusion_file(fusion_annot_lib_path)
    print(f"Finished loading fusion annotation library.")
    print(f"Starting to index FASTA file: {fasta_file_path}")
    fasta_sequences = Fasta(fasta_file_path)
    print(f"Finished indexing FASTA file. {len(fasta_sequences.keys())} sequences indexed.")
    return fusion_annot_dict, fasta_sequences

def parallel_initial_setup(fasta_path: str, gfa_path: str, output_gaf_path: str, fusion_annot_lib: str, total_threads: int, k: int = 12, w: int = 1, m: str = '30,20'):
    align_threads = max(1, total_threads - 1)
    print(f"Allocating {align_threads} threads for alignment and 1 thread for data loading.")
    fusion_cluster = defaultdict(list)
    fusion_annot_dict = {}
    fasta_sequences = None
    with ThreadPoolExecutor(max_workers=2) as executor:
        future_alignment = executor.submit(
            align_and_process_reads,
            fasta_path,
            gfa_path,
            output_gaf_path,
            threads=align_threads,
            k=k,
            w=w,
            m=m
        )
        future_loading = executor.submit(
            _load_reference_data,
            fusion_annot_lib,
            fasta_path
        )
        try:
            print("Waiting for data loading task to complete...")
            fusion_annot_dict, fasta_sequences = future_loading.result()
            print("Data loading task completed.")
        except Exception as e:
            logging.error("An exception occurred during data loading.", exc_info=True)
            future_alignment.cancel()
        try:
            print("Waiting for alignment task to complete...")
            fusion_cluster = future_alignment.result()
            print("Alignment task completed.")
            print("-" * 30)
        except Exception as e:
            logging.error("An exception occurred during alignment.", exc_info=True)
    return fusion_cluster, fusion_annot_dict, fasta_sequences

def first_process(fusion_cluster,gene_positions,ref_path,output_dir,fasta_sequences,skip_chrM=True,only_keep_pcgene=True):
    gene2_cluster = defaultdict(gene2_dict)
    gene2_set = set()
    for key in fusion_cluster.keys():
        if skip_chrM:
            if any(node.path[0].startswith("chrM") for node in fusion_cluster[key]):
                continue
        sorted_cluster = sorted(fusion_cluster[key], key=lambda node: node.read_start_offset)
        has_significant_overlap = check_overlaps(sorted_cluster,overlap = 15)
        genes_similarity = check_genes_similarity(sorted_cluster)
        # genes_similarity = False
        if has_significant_overlap or genes_similarity:
            continue
        merged_cluster = merge_fusion_records(sorted_cluster)
        length_is_short = False
        for record in merged_cluster:
            if record.block_len <= min_block_len:
                length_is_short = True
                break
        if length_is_short or not merged_cluster:
            continue
        fusion_gene_name = [node.path[0].split(':')[1] for node in merged_cluster]

        fusion_len = len(merged_cluster)
        if fusion_len == 2:
            fusion_key = ':'.join(fusion_gene_name)
            gap = merged_cluster[-1].read_start_offset - merged_cluster[0].read_end_offset
            if gap > max_gap_length:
                continue
            if not validate_fusion(gene_positions,merged_cluster[0],merged_cluster[-1],only_keep_pcgene):
                continue
            for i in range(2):
                gene2_cluster[fusion_key][i].append(merged_cluster[i])
        elif fusion_len > 2:
            merged_cluster,truncation_point = truncate_to_first_two_genes(merged_cluster)
            if len(merged_cluster) <= 1:
                continue
            if len(merged_cluster) == 2 and truncation_point is not None:
                fusion_key = f"{fusion_gene_name[0]}:{fusion_gene_name[1]}"
                gap = merged_cluster[-1].read_start_offset - merged_cluster[0].read_end_offset
                if gap > max_gap_length:
                    continue
                if not validate_fusion(gene_positions,merged_cluster[0],merged_cluster[-1],only_keep_pcgene):
                    continue
                for i, record in enumerate(merged_cluster):
                    updated_record = record._replace(read_length=truncation_point)
                    gene2_cluster[fusion_key][i].append(updated_record)
            if len(merged_cluster) == 3:
                all_genes = set()
                for record in merged_cluster:
                    all_genes.add(record.path[0].split(':')[1])
                if len(all_genes) != 3:
                    continue
                ab_gap = merged_cluster[1].read_start_offset - merged_cluster[0].read_end_offset
                if -15 < ab_gap < 15 :
                    ab_fusion_name = [fusion_gene_name[0], fusion_gene_name[1]]
                    ab_fusion_key = ':'.join(ab_fusion_name)
                    ab_len = merged_cluster[1].read_end_offset
                    if not validate_fusion(gene_positions,merged_cluster[0],merged_cluster[1],only_keep_pcgene): continue
                    for i in range(2):
                        updated_record = merged_cluster[i]._replace(read_length=ab_len)
                        gene2_cluster[ab_fusion_key][i].append(updated_record)
                bc_gap = merged_cluster[2].read_start_offset - merged_cluster[1].read_end_offset
                if -15 < bc_gap < 15 :
                    bc_fusion_name = [fusion_gene_name[1], fusion_gene_name[2]]
                    bc_fusion_key = ':'.join(bc_fusion_name)
                    bc_len = merged_cluster[2].read_end_offset
                    if not validate_fusion(gene_positions,merged_cluster[1],merged_cluster[2],only_keep_pcgene): continue
                    for i in range(1,3):
                        updated_record = merged_cluster[i]._replace(read_length=bc_len)
                        gene2_cluster[bc_fusion_key][i-1].append(updated_record)
    keys_to_delete = set()
    checked_pairs = set()
    for key in list(gene2_cluster.keys()):
        geneA, geneB = key.split(':')
        unordered_pair = tuple(sorted([geneA, geneB]))
        if unordered_pair in checked_pairs:
            continue
        checked_pairs.add(unordered_pair)
        key1 = f"{geneA}:{geneB}"
        key2 = f"{geneB}:{geneA}"
        count1 = len(gene2_cluster.get(key1, [[], []])[0])
        count2 = len(gene2_cluster.get(key2, [[], []])[0])
        total_count = count1 + count2
        if total_count < min_supporting_reads:
            if key1 in gene2_cluster:
                keys_to_delete.add(key1)
            if key2 in gene2_cluster:
                keys_to_delete.add(key2)
        else:
            for i in range(2):
                if key1 in gene2_cluster and gene2_cluster[key1][i]:
                    parts = gene2_cluster[key1][i][0].path[0].split(':')
                    gene2_set.add((parts[1], parts[0]))
                if key2 != key1 and key2 in gene2_cluster and gene2_cluster[key2][i]:
                    parts = gene2_cluster[key2][i][0].path[0].split(':')
                    gene2_set.add((parts[1], parts[0]))
    total_reference = os.path.join(output_dir, "total_reference.fasta")
    build_total_reference_fasta(gene2_set, gene_positions, ref_path, total_reference, extension=GENE_EXTENSION)
    validate_gene2_cluster = validate_fusion_genes_with_mappy(gene2_cluster,total_reference,fasta_sequences)
    keys_to_delete = set()
    for key in list(validate_gene2_cluster.keys()):
        geneA, geneB = key.split(':')
        unordered_pair = tuple(sorted([geneA, geneB]))
        if unordered_pair in checked_pairs:
            continue
        checked_pairs.add(unordered_pair)
        key1 = f"{geneA}:{geneB}"
        key2 = f"{geneB}:{geneA}"
        count1 = len(validate_gene2_cluster.get(key1, [[], []])[0])
        count2 = len(validate_gene2_cluster.get(key2, [[], []])[0])
        total_count = count1 + count2
        if total_count < min_supporting_reads:
            if key1 in validate_gene2_cluster:
                keys_to_delete.add(key1)
            if key2 in validate_gene2_cluster:
                keys_to_delete.add(key2)
    for key in keys_to_delete:
        del validate_gene2_cluster[key]
    return validate_gene2_cluster

def main(readfile, ref_path, ref_dir, output_dir, threads, graph_align_params,process_params):
    global min_supporting_reads, min_block_len, max_gap_length, min_cseq_length, bp_cluster_distance, MAX_GAP_TWO_ALIGNMENTS, MAX_GAP_ONE_ALIGNMENT, skip_chrM, only_keep_pcgene
    min_supporting_reads, min_block_len, max_gap_length, min_cseq_length, bp_cluster_distance, MAX_GAP_TWO_ALIGNMENTS, MAX_GAP_ONE_ALIGNMENT, skip_chrM, only_keep_pcgene = process_params
    graph_aln_k, graph_aln_w, graph_aln_m = graph_align_params
    
    temp_dir = os.path.join(output_dir, "middle_files")
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir, exist_ok=True)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    setup_logging(os.path.join(temp_dir, "fusion_detect_middle.log"))
    logging.info("Fusion detection started.")
    fusion_annot_lib = os.path.join(ref_dir, 'fusion_annot_lib')
    gene_pos_file = os.path.join(ref_dir, 'gene_pos.txt')
    comprehensive_gfa = os.path.join(ref_dir, 'comprehensive.gfa')
    only_pc_gfa = os.path.join(ref_dir, 'only_pc.gfa')
    output_gaf = os.path.join(output_dir, "primary_minigraph_alignment.gaf")
    fasta_path = convert_fastq_to_fasta_if_needed(readfile)
    if not fasta_path:
        print("FASTA file preparation failed. Exiting.")
        return
    fusion_cluster, fusion_annot_dict, fasta_sequences = parallel_initial_setup(
        fasta_path=fasta_path,
        gfa_path=comprehensive_gfa,
        output_gaf_path=output_gaf,
        fusion_annot_lib=fusion_annot_lib,
        total_threads=threads,
        k = graph_aln_k,
        w = graph_aln_w,
        m = graph_aln_m
    )
    if not fasta_sequences or not fusion_annot_dict:
        print("Initial data loading or alignment failed. Exiting pipeline.")
        return
    if not fusion_cluster:
        print("No valid fusion data. Exiting..")
        return
    gene_positions = load_gene_positions(gene_pos_file)
    print(f"Total genes loaded: {len(gene_positions)}")
    gene2_cluster = first_process(fusion_cluster,gene_positions,ref_path,temp_dir,fasta_sequences)
    final_fusion = proc_fuison_cluster(gene2_cluster, bp_cluster_distance)
    max_breakpoint_distance = 200
    processor = FusionProcessor(
        gene_positions,
        ref_path,
        only_pc_gfa,
        min_supporting_reads,
        threads,
        max_breakpoint_distance,
        temp_dir
    )

    processor.process_fusion_data(final_fusion, fasta_sequences, fusion_annot_dict, output_dir)
    
