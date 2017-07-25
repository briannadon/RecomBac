#!/usr/bin/python

from Bio import SeqIO,pairwise2,Seq
import re
from sys import argv
import argparse
import primer3

parser = argparse.ArgumentParser(description="Placeholder. Program to define primers in bacterial recombination.")
parser.add_argument('-s', '--sequence', required=True, type=open, help = 'The (currently only genbank supported) sequence to replace.')
parser.add_argument('-ft', '--filetype', required=True, type=str, help = 'The filetype.  Genbank only supported at the moment.')
parser.add_argument('-l', '--leftbound', required=True, type=int, help = 'The left coordinate (1-based) of the sequence to be removed or replaced')
parser.add_argument('-r', '--rightbound', required=True, type=int, help = 'The right coordinate (1-based) of the sequence to be removed or replaced')
parser.add_argument('-lsize', '--leftflanksize', nargs='?', required = False, type=int, default = 1200, help = 'The size of the left flank. Default 1200.')
parser.add_argument('-rsize', '--rightflanksize', nargs = '?', required = False, type=int, default = 1000, help = 'The size of the right flank.  Default 1000.')
parser.add_argument('-o', '--output', nargs= '?', required = False, type = argparse.FileType('w'), help = 'The optional output file.  If excluded, output is to stdout.')
args = parser.parse_args()

#read in the sequence(s)
seqrec = SeqIO.read(args.sequence,args.filetype)
seqrec = str(seqrec.seq)


#designing primers that bind the left flank
leftflank = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'Left Flank',
        'SEQUENCE_TEMPLATE': seqrec,
        'SEQUENCE_INCLUDED_REGION': [args.leftbound - args.leftflanksize, args.leftflanksize],
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
        'PRIMER_NUM_RETURN': 1,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 5,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
    })

if args.output:
     args.output.write('LEFT PRIMER PAIR\n************************\n\n')
     for k in sorted(leftflank.keys()):
          args.output.write('%s\t%s\n' % (k,leftflank[k]))
     args.output.write('\n\n')
else:
     print('LEFT PRIMER PAIR\n************************\n\n')
     for k in sorted(leftflank.keys()):
          print('%s\t%s' % (k,leftflank[k]))
     print('\n\n')

rightflank = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'Right Flank',
        'SEQUENCE_TEMPLATE': seqrec,
        'SEQUENCE_INCLUDED_REGION': [args.rightbound, args.rightflanksize]
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
        'PRIMER_NUM_RETURN': 1,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 5,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
    })
    
if args.output:
     args.output.write('RIGHT PRIMER PAIR\n************************\n\n')
     for k in sorted(rightflank.keys()):
          args.output.write('%s\t%s\n' % (k,rightflank[k]))
     args.output.write('\n\n')
else:
     print('RIGHT PRIMER PAIR\n************************\n\n')
     for k in sorted(rightflank.keys()):
          print('%s\t%s' % (k,rightflank[k]))
     print('\n\n')