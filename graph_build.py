import os
import re
import time
import logging
import traceback
import argparse
import pandas as pd
from collections import defaultdict
from pyfaidx import Fasta
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

def setup_logging(log_path='graph_build.log'):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s', 
        datefmt='%Y-%m-%d %H:%M:%S',  
        filename=log_path, 
        filemode='w'  
    )

def parse_attributes(attr_str):
    attrs = {}
    for attr in attr_str.split(';'):
        if '=' in attr:
            key, value = attr.split('=')
            attrs[key.strip()] = value.strip()
    if 'Dbxref' in attrs:
        dbxref_values = attrs['Dbxref'].split(',')
        dbxref_dict = {}
        for dbxref in dbxref_values:
            parts = dbxref.split(':')
            if len(parts) == 2:
                identifier, value = parts
                dbxref_dict[identifier.strip()] = value.strip()
        attrs['Dbxref'] = dbxref_dict
    else:  
        attrs['Dbxref'] = {}

    return attrs

def init_chromosome_graph():
    chrom_graph = defaultdict(Graph)
    return chrom_graph

def split_gff_into_subregions(gff_path):
    starttime = time.time()
    subregions = []
    current_subregion = []
    current_metadata = None
    columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

    with open(gff_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                if line.startswith('##sequence-region'):
                    if current_subregion:
                        df = pd.DataFrame(current_subregion, columns=columns).assign(
                            start=lambda df: df['start'].astype(int),
                            end=lambda df: df['end'].astype(int)
                        )
                        subregions.append((df, current_metadata))
                        current_subregion = []
                    parts = line.strip().split()
                    current_metadata = {'chrom': parts[1], 'start': int(parts[2]), 'end': int(parts[3])}
                continue
            row = line.strip().split('\t')
            if len(row) < 9:
                continue
            if row[2] not in ['gene', 'transcript', 'exon']:
                continue
            current_subregion.append(row)
        if current_subregion:
            df = pd.DataFrame(current_subregion, columns=columns).assign(
                start=lambda df: df['start'].astype(int),
                end=lambda df: df['end'].astype(int)
            )
            subregions.append((df, current_metadata))
    subregions.sort(key=lambda subregion: len(subregion[0]), reverse=True)
    endtime = time.time()
    logging.info(f'Completed split_gff_into_subregions, Execution time: {endtime - starttime:.2f} seconds')
    print(f'Completed split_gff_into_subregions, Execution time: {endtime - starttime:.2f} seconds')
    return subregions

def convert_node_id(node_id):
    match = re.match(r'(.+?)_(?:gene|pseu)-(.+?)_exon_(\d+)_(\d+)', node_id)
    if match:
        chrom, gene, start, end = match.groups()
        return f"{chrom}:{gene}:{start}-{end}"
    else:
        return node_id

class Node:
    def __init__(self, node_id, row, attrs, parent_id=None, sequence=None, path_num=None, seqid=None):
        self.node_id = node_id
        self.row = row
        self.attrs = attrs
        self.parent_id = parent_id
        self.sequence = sequence
        self.path_num = path_num
        self.seqid = seqid
        if row is not None and attrs is not None:
            self.node_attrs = self._create_node_attrs()

    def _create_node_attrs(self):
        node_attrs = {
            'type': self.row.type,
            'start': self.row.start,
            'end': self.row.end,
            'ID': self.attrs['ID'],
            'GeneID': [],
            'strand': self.row.strand,
            'seq': self.sequence,
            'path_nums': [self.path_num],
            'gene_name': self.attrs.get('gene_name', ''),
            'transcript_id': self.attrs.get('transcript_id', ''),
            'repeat_node_id': False,
            'repeat_node_parent_id': [self.parent_id],
            'pre_exon': {},
            'next_exon': {}
        }
        if self.path_num is not None:
            node_attrs['pre_exon'][self.path_num] = None
            node_attrs['next_exon'][self.path_num] = None
        if 'Parent' in self.attrs.keys():
            node_attrs['Parent'] = self.attrs['Parent']
        if 'GeneID' in self.attrs['Dbxref'].keys():
            node_attrs['GeneID'].append(self.attrs['Dbxref']['GeneID'])
        return node_attrs

class Graph:
    def __init__(self, name=None):
        self.name = name
        self.nodes = {} 
        self.edges = {}  
        self.path_num = None
        self.gene_sequence = None
        self.gene_start = None 
        self.gene_end = None   
        self.strand = None     
        self.gene_type = None  # Added gene_type field

    def add_node(self, node):
        node_id = node.node_id
        if node.attrs:
            attrs = node.node_attrs
            if node_id not in self.nodes:
                self.nodes[node_id] = attrs
            else:
                self.nodes[node_id]['path_nums'].append(node.path_num)
                self.nodes[node_id]['repeat_node_id'] = True

    def add_chrome_node(self, node):
        node_id = node.node_id
        self.nodes[node_id] = node_id

    def add_edge(self, source, target):
        assert source in self.nodes and target in self.nodes, "Both nodes must exist in the graph."
        if target not in self.edges.get(source, []):
            self.edges.setdefault(source, []).append(target)

    def get_exon_edges(self):
        exon_edges = set()
        for source, targets in self.edges.items():
            for target in targets:
                if (self.nodes[source].get('type') == 'exon' and
                        self.nodes[target].get('type') == 'exon'):
                    exon_edges.add((source, target))
        return list(exon_edges)

    def add_nodes_from(self, nodes_data):
        for node_id, attrs in nodes_data:
            if node_id not in self.nodes:
                self.nodes[node_id] = attrs

    def add_edges_from(self, edges_data):
        for source, target in edges_data:
            self.add_edge(source, target)

    def clear(self):
        self.nodes.clear()
        self.edges.clear()

def netw_to_gfa(chr_name, subGraph, output_dir):
    all_nodes = list(subGraph.nodes.keys())
    if not all_nodes:
        logging.info(f"Warning: The graph has no nodes:,{subGraph.name}")
        return

    gene_type = subGraph.name.split('_')[0]
    subGraph_name = subGraph.name
    if gene_type == 'pseu':
        gene_type = 'pseudogene'
        subGraph_name = subGraph_name.replace('pseu', 'pseudogene')
    
    GFA_filename = f"{subGraph_name}.gfa"
    exon_nodes = [node_id for node_id, attrs in subGraph.nodes.items() if attrs.get('type') == 'exon']
    gene_node = next((node_id for node_id, attrs in subGraph.nodes.items() if attrs.get('type') == 'gene'), None)
    gene_name = "unknown_gene"
    
    if gene_node:
        gene_attrs = subGraph.nodes[gene_node]
        gene_name = gene_attrs.get('ID', gene_attrs.get('gene_name', gene_name))
        gene_start = gene_attrs.get('start')
        gene_sequence = subGraph.gene_sequence
    
    exon_nodes.sort(key=lambda node: int(node.split('_')[-1]))
    strand = subGraph.nodes[all_nodes[0]].get('strand')
    
    # Determine the correct output subfolder based on gene type and if it's PC or not
    if subGraph.is_protein_coding:
        output_folder = os.path.join(output_dir, 'pc')
    else:
        output_folder = os.path.join(output_dir, 'non_pc')
        
    os.makedirs(output_folder, exist_ok=True)
    gfa_path = os.path.join(output_folder, GFA_filename)

    try:
        with open(gfa_path, "w") as gfa_file:
            if len(all_nodes) == 1:
                node_id = all_nodes[0]
                sequence = subGraph.nodes[node_id].get('seq', '')
                current_type = subGraph.nodes[node_id].get('type')
                start = subGraph.nodes[node_id].get('start')
                end = subGraph.nodes[node_id].get('end')
                if subGraph.nodes[node_id].get('type') == 'pseudogene':
                    gene_name = subGraph.nodes[node_id].get('ID').replace("gene", "pseu")
                else:
                    gene_name = subGraph.nodes[node_id].get('ID')
                node_id = f"{chr_name}_{gene_name}_{current_type}_{start}_{end}"
                node_id_new = convert_node_id(node_id)
                gfa_file.write(f"S\t{node_id_new}\t{sequence}\n")
                return

            if exon_nodes:
                # S行
                for node_id in exon_nodes:
                    sequence = subGraph.nodes[node_id].get('seq', '')
                    node_id_new = convert_node_id(node_id)
                    gfa_file.write(f"S\t{node_id_new}\t{sequence}\n")

                exon_edges = subGraph.get_exon_edges()
                if strand == '-':
                    swapped_edges = [(v, u) for u, v in exon_edges]
                elif strand == '+':
                    swapped_edges = exon_edges
                else:
                    logging.info(f"Warning: Unknown strand '{strand}', skipping edge sorting.")
                    swapped_edges = exon_edges

                # L
                for source_id, target_id in swapped_edges:
                    source_strand = '+'
                    target_strand = '+'
                    source_id_new = convert_node_id(source_id)
                    target_id_new = convert_node_id(target_id)
                    gfa_file.write(f"L\t{source_id_new}\t{source_strand}\t{target_id_new}\t{target_strand}\t0M\n")

                # P
                if len(exon_edges) >= 1:
                    path_num = subGraph.path_num
                    gene_direction = '+'
                    path_record = {}
                    for path in range(1, path_num + 1):
                        if path not in path_record:
                            path_record[path] = []
                        for node_id in exon_nodes:
                            if subGraph.nodes[node_id]['repeat_node_id']:
                                if path in subGraph.nodes[node_id]['path_nums']:
                                    path_record[path].append(node_id)
                            else:
                                if path == subGraph.nodes[node_id]['path_nums'][0]:
                                    path_record[path].append(node_id)

                    def sort_path_record(path_record, direction='+'):
                        sorted_path_record = {}
                        for path, nodes in path_record.items():
                            if direction == '+':
                                sorted_nodes = sorted(nodes, key=lambda x: int(x.split('_')[-2]))
                            elif direction == '-':
                                sorted_nodes = sorted(nodes, key=lambda x: int(x.split('_')[-2]), reverse=True)
                            else:
                                raise ValueError("Direction must be '+' or '-'")
                            sorted_path_record[path] = sorted_nodes
                        return sorted_path_record

                    sorted_path_record = sort_path_record(path_record, direction=gene_direction)
                    for path, nodes in sorted_path_record.items():
                        segment_names = ','.join([convert_node_id(node) for node in nodes])
                        if nodes:
                            path_id = "_".join(nodes[0].split("_")[:2])
                        else:
                            path_id = "unknown"
                        gfa_file.write(f"P\t{path_id}_{path}\t{segment_names}\n")

    except IOError as e:
        logging.info(f"Error writing to GFA file '{gfa_path}': {e}")

def build_gencode_graph(subregion, fna_file_path, output_dir):
    df = subregion[0]
    metadata = subregion[1]
    start_time = time.time()
    gene_sequence_cache = {} 
    CDS_exon_nodeid_to_id = {}
    CDS_exon_id_to_nodeid = {}
    gene_positions = {}
    fasta = Fasta(fna_file_path, as_raw=True) 
    chromosome = metadata['chrom']
    if chromosome not in fasta:
        raise ValueError(f"Chromosome {chromosome} not found in ref file.")
    chromosome_seq = fasta[chromosome]
    chromosome_subgraph = Graph()
    root_id = chromosome
    chromosome_subgraph.name = root_id
    chromosome_subgraph.add_chrome_node(Node(root_id, row=None, attrs=None))

    try:
        previous_row = None
        next_row = df.iloc[1] if len(df) > 1 else None
        for i, current_row in enumerate(df.itertuples()):
            next_row = df.iloc[i + 1] if i < len(df) - 1 else None
            try:
                if i > 0:
                    pre_attrs = parse_attributes(previous_row.attributes)
                    pre_type = previous_row.type
                    pre_ID = pre_attrs['ID']
                else:
                    pre_type = None
                    pre_ID = None
                if next_row is not None:
                    next_attrs = parse_attributes(next_row.attributes)
                    next_GeneID = next_attrs.get('gene_id', None)
                else:
                    next_GeneID = None
                attrs = parse_attributes(current_row.attributes)
                current_type = current_row.type
                now_GeneID = attrs.get('gene_id', None)
                strand = current_row.strand
                if current_type == 'gene' and now_GeneID not in gene_sequence_cache:
                    try:
                        gene_sequence = chromosome_seq[current_row.start - 1:current_row.end]
                        gene_sequence_cache[now_GeneID] = gene_sequence
                    except KeyError as e:
                        logging.error(f"Error in extracting sequence for GeneID {now_GeneID}: {e}")
                if now_GeneID:
                    if current_type == 'gene':
                        gene_start = current_row.start
                        gene_end = current_row.end
                        graph = Graph()
                        gene_node_id = attrs['ID']
                        gene_sequence = gene_sequence_cache.get(now_GeneID, '')
                        graph.gene_sequence = gene_sequence
                        graph.gene_start = current_row.start 
                        graph.gene_end = current_row.end     
                        graph.strand = strand                
                        gene_name = attrs.get('gene_name', gene_node_id)
                        
                        # Check if this is a protein-coding gene
                        is_protein_coding = attrs.get('gene_type', '') == 'protein_coding'
                        graph.is_protein_coding = is_protein_coding
                        
                        if 'pseudogene' in attrs.get('gene_type', 'gene'):
                            gene_type = 'pseu'
                        else:
                            gene_type = 'gene'
                        graph.name = f"{gene_type}_{gene_name}_{chromosome}_{current_row.start}_{current_row.end}"
                        gene_positions[f"{gene_name}_{chromosome}"] = {
                            'start': current_row.start,
                            'end': current_row.end,
                            'type': attrs.get('gene_type', ''),
                            'Fexons': set(),
                            'Bexons': set(),
                            'strand': strand  
                        }
                        path_num = 0
                        graph.path_num = 0

                        graph.add_node(Node(gene_node_id, current_row, attrs, sequence=gene_sequence, seqid=chromosome))
                        if chromosome_subgraph:
                            chromosome_subgraph.add_node(Node(gene_node_id, current_row, attrs))
                            chromosome_subgraph.add_edge(root_id, gene_node_id)
                    elif current_type == 'transcript':
                        node_id = attrs['ID']
                        transcript_id = attrs.get('transcript_id', node_id)
                        parent_node_id = attrs.get('Parent', gene_node_id)
                        graph.add_node(Node(node_id, current_row, attrs, parent_id=parent_node_id))
                        graph.add_edge(parent_node_id, node_id)
                        path_num += 1  
                    elif current_type == 'exon':
                        node_id = f"{current_row.seqid}_{gene_type}-{gene_name}_{current_type}_{current_row.start}_{current_row.end}"
                        CDS_exon_nodeid_to_id[node_id] = attrs['ID']
                        CDS_exon_id_to_nodeid[attrs['ID']] = node_id
                        parent_node_id = attrs['Parent']
                        exon_start, exon_end = current_row.start, current_row.end
                        gene_positions[f"{gene_name}_{chromosome}"]['Fexons'].add(exon_start)
                        gene_positions[f"{gene_name}_{chromosome}"]['Bexons'].add(exon_end)
                        # 提取外显子序列
                        if now_GeneID in gene_sequence_cache:
                            gene_sequence = gene_sequence_cache[now_GeneID]
                            exon_seq = gene_sequence[exon_start - gene_start: exon_end - gene_start + 1]
                        else:
                            logging.info(f"GeneID {now_GeneID} not in gene_sequence_cache")
                            exon_seq = None
                        # 添加外显子节点到图
                        graph.add_node(
                            Node(node_id, current_row, attrs, parent_id=parent_node_id, sequence=exon_seq, path_num=path_num)
                        )
                        # 添加路径
                        if pre_type == 'exon':
                            pre_exon_node_id = CDS_exon_id_to_nodeid[pre_ID]
                            graph.add_edge(pre_exon_node_id, node_id)
                            graph.nodes[gene_node_id]['pre_exon'][path_num] = pre_exon_node_id
                            graph.nodes[pre_exon_node_id]['next_exon'][path_num] = node_id
                        else:
                            graph.add_edge(parent_node_id, node_id)
                    elif current_type in ['CDS', 'start_codon', 'stop_codon', 'three_prime_UTR', 'five_prime_UTR','stop_codon_redefined_as_selenocysteine']:
                        # 处理其他类型的节点（如 CDS、UTR）
                        pass
                    else:
                        pass
                else:
                    logging.warning(f"Skipping row with no gene_id: {current_row}")
                    graph = None
                if graph and now_GeneID and now_GeneID != next_GeneID:
                    graph.path_num = path_num
                    # Always process all genes and write to appropriate directory
                    netw_to_gfa(root_id, graph, output_dir)
                    graph.clear()
                    CDS_exon_nodeid_to_id.clear()
                    CDS_exon_id_to_nodeid.clear()
                    gene_sequence_cache.clear()
                
                previous_row = current_row 

            except KeyError as e:
                logging.error(f"KeyError in processing row: {e}")
            except Exception as e:
                logging.error(f"Unexpected error in graphthread: {e}\n{traceback.format_exc()}")
    except Exception as e:
        logging.error(f"An error occurred in mm_thread: {e}")
    finally:
        end_time = time.time()
        logging.info(f'Complete region: {df.iloc[0, 0]}, Runtime: {end_time - start_time} seconds')
        print(f'Complete region: {df.iloc[0, 0]}, Runtime: {end_time - start_time} seconds')
    return chromosome, chromosome_subgraph, gene_positions

def write_gene_positions(gene_positions, output_dir):
    gene_pos_path = os.path.join(output_dir, 'gene_pos.txt')
    with open(gene_pos_path, 'w') as f:
        f.write("Gene_Name\tStrand\tType\tStart\tEnd\tFexons\tBexons\n")
        for gene_key, position in gene_positions.items():
            Fexons_str = ';'.join(sorted(map(str, position['Fexons'])))
            Bexons_str = ';'.join(sorted(map(str, position['Bexons'])))
            f.write(f"{gene_key}\t{position['strand']}\t{position['type']}\t{position['start']}\t{position['end']}\t{Fexons_str}\t{Bexons_str}\n")

    logging.info(f"Successfully wrote {len(gene_positions)} gene position information to {gene_pos_path}")

def process_subregions(subregions, chromosome_graph, fna_file_path, output_dir, num_threads):
    start_time = time.time()
    logging.info(f'Start building subgraphs, number of threads:{num_threads}')
    all_gene_positions = {}

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(build_gencode_graph, subregion, fna_file_path, output_dir) for subregion in subregions]
        for future in as_completed(futures):
            try:
                chromosome, graph, gene_positions = future.result()
                all_gene_positions.update(gene_positions)
                node_items = list(graph.nodes.items())
                edge_items = [(source, target) for source, targets in graph.edges.items() for target in targets]
                chromosome_graph[chromosome].add_nodes_from(node_items)
                chromosome_graph[chromosome].add_edges_from(edge_items)
            except Exception as exc:
                logging.error(f'The task generates an exception: {exc}')
    
    # Always write gene positions for all genes
    write_gene_positions(all_gene_positions, output_dir)
    
    end_time = time.time()
    logging.info(f'Finished building subgraph, took: {end_time - start_time} seconds')

def read_and_filter_content(file_path):
    try:
        with open(file_path, 'r') as file:
            return [line for line in file if line.startswith('S') or line.startswith('L')]
        
    except IOError as e:
        logging.info(f"Error reading file {file_path}: {e}")
        return []

def merge_gfa_files(input_folder, output_file, num_threads=8):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    gfa_files = [
        entry.path for entry in os.scandir(input_folder)
        if entry.is_file() and entry.name.endswith('.gfa')
    ]
    filtered_contents = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = executor.map(read_and_filter_content, gfa_files)
        for content in results:
            filtered_contents.extend(content)
    try:
        with open(output_file, 'w') as output_file_handle:
            output_file_handle.writelines(filtered_contents)
        logging.info(f"Successfully created GFA file: {output_file}")
    except IOError as e:
        logging.info(f"Error writing to file {output_file}: {e}")

def main(gff_path,ref_path,num_threads,output_dir,log_file):
    os.makedirs(output_dir, exist_ok=True)
    
    # Create PC and non-PC subdirectories
    pc_dir = os.path.join(output_dir, 'pc')
    non_pc_dir = os.path.join(output_dir, 'non_pc')
    os.makedirs(pc_dir, exist_ok=True)
    os.makedirs(non_pc_dir, exist_ok=True)
    
    # Setup logging with custom path
    log_file = os.path.join(output_dir, log_file)
    setup_logging(log_file)

    # gff_path = parsed_args.gff_path
    fna_file_path = ref_path
    # num_threads = parsed_args.num_threads

    start_time = time.time()
    logging.info('Start processing gene annotation data')
    print('Start processing gene annotation data')
    # Initialize chromosome_graph
    chromosome_graph = init_chromosome_graph()
    logging.info('Initialization of chromosome_graph completed')

    # Split GFF file
    subregions = split_gff_into_subregions(gff_path)
    logging.info(f'The annotation file is split into {len(subregions)} subregions')
    print(f'The annotation file is split into {len(subregions)} subregions')
    # Process subregions - now passes output_dir to the function
    process_subregions(subregions, chromosome_graph, fna_file_path, output_dir, num_threads)
    logging.info('Sub-region processing completed')
    print('Sub-region processing completed')
    
    # Create output paths for merged GFA files
    pc_only_gfa = os.path.join(output_dir, 'only_pc.gfa')
    comprehensive_gfa = os.path.join(output_dir, 'comprehensive.gfa')
    
    # Merge PC GFA files to create PC-only file
    merge_gfa_files(pc_dir, pc_only_gfa, num_threads)
    logging.info('PC GFA file merge completed')
    print(f"Successfully created GFA file: {pc_only_gfa}")
    
    # Merge PC and non-PC GFA files for comprehensive file
    # First, merge PC files
    pc_contents = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        pc_files = [entry.path for entry in os.scandir(pc_dir) if entry.is_file() and entry.name.endswith('.gfa')]
        results = executor.map(read_and_filter_content, pc_files)
        for content in results:
            pc_contents.extend(content)
    
    # Then merge non-PC files
    non_pc_contents = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        non_pc_files = [entry.path for entry in os.scandir(non_pc_dir) if entry.is_file() and entry.name.endswith('.gfa')]
        results = executor.map(read_and_filter_content, non_pc_files)
        for content in results:
            non_pc_contents.extend(content)
    
    # Write comprehensive GFA file
    try:
        with open(comprehensive_gfa, 'w') as output_file:
            output_file.writelines(pc_contents + non_pc_contents)
        logging.info(f"Successfully created GFA file: {comprehensive_gfa}")
        print(f"Successfully created GFA file: {comprehensive_gfa}")
    except IOError as e:
        logging.info(f"Error writing to comprehensive GFA file: {e}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    logging.info(f'graph_build total runtime: {elapsed_time:.2f} seconds')
    print(f'graph_build total runtime: {elapsed_time:.2f} seconds')
    
    # clean middle files
    print('Deleting middle files')
    os.system('rm -rf ' + pc_dir)
    os.system('rm -rf ' + non_pc_dir)
    
    
    
def main_argparse(args=None):
    parser = argparse.ArgumentParser(description='Process gene annotations to build graph representations')
    parser.add_argument(
        '--gff_path',
        type=str,
        required=True,
        help='Path to gff3/gtf annotation file'
    )
    parser.add_argument(
        '--fna_file_path',
        type=str,
        required=True,
        help='Path to reference genome fasta file'
    )
    parser.add_argument(
        '--num_threads',
        type=int,
        default=8,
        help='Number of threads to use, default is 8'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        default='output',
        help='Directory to store output files'
    )
    parser.add_argument(
        '--log_file',
        type=str,
        default='graph_build.log',
        help='Path to log file'
    )

    if args is not None:
        parsed_args = parser.parse_args(args)
    else:
        parsed_args = parser.parse_args()

    # Set up output directories
    output_dir = parsed_args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    # Create PC and non-PC subdirectories
    pc_dir = os.path.join(output_dir, 'pc')
    non_pc_dir = os.path.join(output_dir, 'non_pc')
    os.makedirs(pc_dir, exist_ok=True)
    os.makedirs(non_pc_dir, exist_ok=True)
    
    # Setup logging with custom path
    log_file = os.path.join(output_dir, parsed_args.log_file)
    setup_logging(log_file)

    gff_path = parsed_args.gff_path
    fna_file_path = parsed_args.fna_file_path
    num_threads = parsed_args.num_threads

    start_time = time.time()
    logging.info('Start processing gene annotation data')

    # Initialize chromosome_graph
    chromosome_graph = init_chromosome_graph()
    logging.info('Initialization of chromosome_graph completed')

    # Split GFF file
    subregions = split_gff_into_subregions(gff_path)
    logging.info(f'The GFF file is split into {len(subregions)} subregions')

    # Process subregions - now passes output_dir to the function
    process_subregions(subregions, chromosome_graph, fna_file_path, output_dir, num_threads)
    logging.info('Sub-region processing completed')

    # Create output paths for merged GFA files
    pc_only_gfa = os.path.join(output_dir, 'only_pc.gfa')
    comprehensive_gfa = os.path.join(output_dir, 'comprehensive.gfa')
    
    # Merge PC GFA files to create PC-only file
    merge_gfa_files(pc_dir, pc_only_gfa, num_threads)
    logging.info('PC GFA file merge completed')
    
    # Merge PC and non-PC GFA files for comprehensive file
    # First, merge PC files
    pc_contents = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        pc_files = [entry.path for entry in os.scandir(pc_dir) if entry.is_file() and entry.name.endswith('.gfa')]
        results = executor.map(read_and_filter_content, pc_files)
        for content in results:
            pc_contents.extend(content)
    
    # Then merge non-PC files
    non_pc_contents = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        non_pc_files = [entry.path for entry in os.scandir(non_pc_dir) if entry.is_file() and entry.name.endswith('.gfa')]
        results = executor.map(read_and_filter_content, non_pc_files)
        for content in results:
            non_pc_contents.extend(content)
    
    # Write comprehensive GFA file
    try:
        with open(comprehensive_gfa, 'w') as output_file:
            output_file.writelines(pc_contents + non_pc_contents)
        logging.info(f"Successfully created comprehensive GFA file: {comprehensive_gfa}")
    except IOError as e:
        logging.info(f"Error writing to comprehensive GFA file: {e}")

    # Finish timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'graph_build running time: {elapsed_time:.2f} seconds')
    logging.info(f'graph_build total runtime: {elapsed_time:.2f} seconds')
    
    # clean middle files
    os.system('rm -rf ' + pc_dir)
    os.system('rm -rf ' + non_pc_dir)

if __name__ == "__main__":
    main_argparse()
    
