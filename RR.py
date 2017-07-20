#!/usr/bin/python

from Bio import SeqIO,pairwise2,Seq
import re
from sys import argv
import argparse
import primer3

script, gb, left, right = argv
left, right = int(left), int(right)
#parser = argparse.ArgumentParser(description="Placeholder. Program to define primers in bacterial recombination.")

#read in the sequence(s)
seqrec = SeqIO.read(gb,'genbank')
seqrec = str(seqrec.seq)
seqincludedregion = [left,right]

inner = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': seqrec,
        'SEQUENCE_INCLUDED_REGION': seqincludedregion
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_TASK': 'pick_cloning_primers',
        'PRIMER_PICK_ANYWAY': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 26,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 65.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 5,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
    })

for k in sorted(inner.keys()):
     print('%s\t%s' % (k,inner[k]))

seqleftflank = [left-1000,left]
seqrightflank = [right, right + 1000]

leftflank = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': seqrec,
        'SEQUENCE_INCLUDED_REGION': seqleftflank
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_TASK': 'pick_cloning_primers',
        'PRIMER_PICK_ANYWAY': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 26,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 65.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 5,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
    })

rightflank = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': seqrec,
        'SEQUENCE_INCLUDED_REGION': seqleftflank
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_TASK': 'pick_cloning_primers',
        'PRIMER_PICK_ANYWAY': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 26,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 65.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 5,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
    })
    
