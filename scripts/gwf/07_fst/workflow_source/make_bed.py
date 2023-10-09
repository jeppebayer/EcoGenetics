#!/bin/env python

import sys

def parse_fasta(fasta_file: str):
    """Parses `FASTA` file returning all sequence names and lengths paired in a list of dictionaries.
    
    ::
    
        return [{'sequence_name': str, 'sequence_length': int}, ...]
    
    :param str fasta_file:
        Sequence file in `FASTA` format.
    """
    fasta_list = []
    seq_name = None
    length = 0
    with open(fasta_file, 'r') as fasta:
        for entry in fasta:
            entry = entry.strip()
            if entry.startswith(">"):
                if seq_name:
                    fasta_list.append({'sequence_name': seq_name, 'sequence_length': length})
                    length = 0
                entry = entry.split(" ", 1)
                seq_name = entry[0][1:]
            else:
                length += len(entry)
        fasta_list.append({'sequence_name': seq_name, 'sequence_length': length})
    return fasta_list

def partition_chrom(parse_fasta: list, size: int = 500000):
    """
    Partitions `FASTA` file parsed with **parse_fasta**.
    
    Uses the list of dictionaries from **parse_fasta** to creates a list of dictionaries
    containing with partition number, sequence name, start and end postion (0 based).
    By default the partition size is 500kbs.
    
    ::
    
        return [{'num': int, 'region': str, 'start': int, 'end': int}, ...]
    
    :param list parse_fasta:
        List of dictionaries produced by **parse_fasta**.
    :param int size:
        Size of partitions. Default 500kb.
    """
    chrom_partition = []
    num = 1
    for chrom in parse_fasta:
        whole_chunks = chrom['sequence_length'] // size
        partial_chunk = chrom['sequence_length'] - whole_chunks * size
        start = 0
        for chunk in range(whole_chunks):
            end = start + size
            chrom_partition.append({'num': num, 'region': chrom['sequence_name'], 'start': start, 'end': end})
            start = end
            num += 1
        chrom_partition.append({'num': num, 'region': chrom['sequence_name'], 'start': start, 'end': start + partial_chunk})
        num += 1
    return chrom_partition

def make_bed(partitions: list, output_directory: str = '.'):
    """
    Creates :format:`BED` file for each partition from **partition_chrom**.

    Uses the partitions from **partition_chrom** to create a series of :format:`BED` files for each partition.
    The :format:`BED`files will be named: 'num'.bed, where 'num' corresponds to the partitions number.
    Each :format:`BED`files contains one line with three tab-separated columns: 'chromosome', 'start', 'end', with coordinates being 0-based.
    Returns list of all files created.

    ::

        return ['1.bed', '2.bed',... 'n.bed']

    :param list partitions:
        List of dictionaries produced by **partitions_chrom**
    :param str output_directory:
        Directory for output :format:`BED`files.
    """
    file_list = []
    for part in partitions:
        with open('{dir_path}/{num}.bed'.format(dir_path=output_directory, num=part['num']), 'w') as bed_file:
            bed_file.write('{chrom}\t{start}\t{end}'.format(chrom=part['region'], start=part['start'], end=part['end']))
        file_list.append('{dir_path}/{num}.bed'.format(dir_path=output_directory, num=part['num']))
    return file_list

fasta_file = sys.argv[1]
size = sys.argv[2]
output_directory = sys.argv[3]

make_bed(partitions=partition_chrom(parse_fasta=parse_fasta(fasta_file), size=size), output_directory=output_directory)