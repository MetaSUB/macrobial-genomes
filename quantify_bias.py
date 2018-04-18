#! /usr/bin/env python3

import click
import gzip

from glob import glob
from json import dumps
from sys import stdout
from os.path import dirname


GENOME_PATTERN = '100M_base_limited_genomes/*.100M.unzip.fa'
M8_PATTERN = 'blast_tabulars/*.m8'
GENOME_INFO = 'macrobe_genomes.csv'


def get_real_names():
    name_tbl = {}
    with open(GENOME_INFO) as gf:
        for line in gf:
            tkns = line.split(',')
            common_name = tkns[0]
            gname = tkns[2].split('.')[0].strip()
            name_tbl[gname] = common_name
    return name_tbl


def parse_fasta_len(fastaf):
    tot_bases = 0
    chrs = []
    with open(fastaf) as ff:
        for line in ff:
            line = line.strip()
            if line[0] != '>':
                tot_bases += len(line)
            else:
                chr_name = line[1:].split()[0]
                chrs.append(chr_name)
    return tot_bases, chrs


def get_raw_genome_lengths():
    length_tbl = {}
    chr_tbl = {}
    for fastaf in glob(GENOME_PATTERN):
        gname = fastaf.split('/')[-1].split('.')[0]
        tot_bases, chrs = parse_fasta_len(fastaf)
        length_tbl[gname] = tot_bases
        chr_tbl[gname] = chrs
    return length_tbl, chr_tbl


def parse_m8_line(line):
    tkns = line.strip().split('\t')
    perc_id = float(tkns[2])
    qstart = int(tkns[6])
    qend = int(tkns[7])
    s, e = min(qend, qstart), max(qend, qstart)
    return perc_id, s, e


def len_of_similar_stretches(m8fname, sim_cutoff):
    sections = []
    with open(m8fname) as m8f:
        for line in m8f:
            perc_id, start, end = parse_m8_line(line)
            if perc_id >= sim_cutoff:
                sections.append((start, end))
    sections = sorted(sections, key=lambda el: el[0])

    tot_bases = 0
    cur_start, cur_end = sections[0]
    for start, end in sections[1:]:
        if start <= cur_end:  # section overlaps
            cur_end = end
        else:  # sections do not overlap
            tot_bases += cur_end - cur_start
            cur_start, cur_end = start, end
    tot_bases += cur_end - cur_start

    return tot_bases


def get_similar_stretches(sim_cutoff):
    sim_tbl = {}
    for m8fname in glob(M8_PATTERN):
        gname = m8fname.split('/')[-1].split('.')[0]
        sim_tbl[gname] = len_of_similar_stretches(m8fname, sim_cutoff)
    return sim_tbl


def get_dissimilar_lengths(len_tbl, sim_cutoff=0.9):
    dis_tbl = {}
    sim_tbl = get_similar_stretches(sim_cutoff)
    for gname, qlen in len_tbl.items():
        dis_len = qlen - sim_tbl[gname]
        dis_tbl[gname] = dis_len
    return dis_tbl


@click.command()
@click.option('-s', '--similarity', default=0.9)
def main(similarity):
    len_tbl, chr_tbl = get_raw_genome_lengths()
    dis_tbl = get_dissimilar_lengths(len_tbl, sim_cutoff=similarity)
    name_tbl = get_real_names()
    out_tbl = {}
    for gname, dis in dis_tbl.items():
        common_name = name_tbl[gname]
        out_tbl[gname] = {
            'common_name': common_name,
            'raw_length': len_tbl[gname],
            'effective_length': dis,
            'chrs': chr_tbl[gname],
        }
    stdout.write(dumps(out_tbl))


if __name__ == '__main__':
    main()
