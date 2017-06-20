#!/usr/bin/python

from Bio import SeqIO,pairwise2
import re
from sys import argv
import argparse
import primer3

script, gb, left, right = argv
left, right = int(left), int(right)
#parser = argparse.ArgumentParser(description="Placeholder. Program to define primers in bacterial recombination.")

seqrec = SeqIO.read(gb,'genbank')
seqrec = str(seqrec.seq)
seqincludedregion = [left,right]

primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': seqrec,
        'SEQUENCE_INCLUDED_REGION': seqincludedregion
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_TASK': 'pick_cloning_primers',
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
    })