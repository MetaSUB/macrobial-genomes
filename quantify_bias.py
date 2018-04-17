#! /usr/bin/env python3

import click
import gzip

from glob import glob
from json import dumps
from sys import stdout
from os.path import dirname


GENOME_PATTERN = dirname(__file__) + '/100M_base_limited_genomes/*.100M'
M8_PATTERN = dirname(__file__) + '/blast_tabulars/*.m8'
GENOME_INFO = dirname(__file__) + 'macrobe_genomes.csv'


def get_real_names():
    name_tbl = {}
    with open(GENOME_INFO) as gf:
        for line in gf:
            tkns = line.split(',')
            common_name = tkns[0]
            gname = tkns[2].split('.')[0]
            name_tbl[gname] = common_name
    return name_tbl


def parse_fasta_len(fastaf):
    tot_bases = 0
    with gzip.open(fastaf, mode='r') as ff:
        for line in ff:
            if line[0] != '>':
                tot_bases += len(line)
    return tot_bases


def get_raw_genome_lengths():
    length_tbl = {}
    for fastaf in glob(GENOME_PATTERN):
        gname = fastaf.split('/')[-1].split('.')[0]
        length_tbl[gname] = parse_fasta_len(fastaf)
    return length_tbl


def parse_m8_line(line):
    tkns = line.strip().split('\t')
    perc_id = float(tkns[2])
    qstart = int(tkns[6])
    qend = int(tkns[7])
    qlen = max(qend, qstart) - min(qend, qstart)
    return perc_id, qlen


def len_of_similar_stretches(m8fname, sim_cutoff):
    tot_bases = 0
    with open(m8fname) as m8f:
        for line in m8f:
            perc_id, qlen = parse_m8_line(line)
            if perc_id >= sim_cutoff:
                tot_bases += qlen
    return tot_bases


def get_similar_stretches(sim_cutoff):
    sim_tbl = {}
    for m8fname in glob(M8_PATTERN):
        gname = m8fname.split('/')[-1].split('.')[0]
        sim_tbl[gname] = len_of_similar_stretches(m8fname, sim_cutoff)
    return sim_tbl


def get_dissimilar_lengths(sim_cutoff=0.9):
    dis_tbl = {}
    len_tbl = get_raw_genome_lengths()
    sim_tbl = get_similar_stretches(sim_cutoff)
    for gname, qlen in len_tbl.items():
        dis_len = sim_tbl[gname]
        dis_tbl[gname] = dis_len
    return dis_tbl


@click.command()
@click.option('-s', '--similarity', default=0.9)
def main(similarity):
    dis_tbl = get_dissimilar_lengths(sim_cutoff=similarity)
    name_tbl = get_real_names()
    out_tbl = {}
    for gname, dis in dis_tbl.items():
        common_name = name_tbl[gname]
        out_tbl[gname] = {
            'common_name': common_name,
            'effective_length': dis,
        }
    stdout.write(dumps(out_tbl))


if __name__ == '__main__':
    main()
