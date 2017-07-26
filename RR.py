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

def print_primers(lprimer, rprimer):
     for k in sorted(lprimer.keys()):
          print('%s\t%s' % (k,lprimer[k]))
     for k in sorted(rprimer.keys()):
          print('%s\t%s' % (k,rprimer[k]))
          
def write_primers(lprimer, rprimer, wfile):
     for k in sorted(lprimer.keys()):
          args.output.write('%s\t%s\n' % (k,lprimer[k]))
     for k in sorted(rprimer.keys()):
          args.output.write('%s\t%s\n' % (k,rprimer[k]))
     args.output.write('\n\n')
     

#designing primers that bind the left flank
lflanklprimer = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'Left Flank',
        'SEQUENCE_TEMPLATE': seqrec,
        'SEQUENCE_INCLUDED_REGION': [args.leftbound - args.leftflanksize, args.leftflanksize],
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 0,
        'PRIMER_PICK_ANYWAY': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 15,
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

lflankrprimer = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'Left Flank',
        'SEQUENCE_TEMPLATE': seqrec,
        'SEQUENCE_INCLUDED_REGION': [args.leftbound - args.leftflanksize, args.leftflanksize],
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_TASK': 'pick_cloning_primers',
        'PRIMER_PICK_LEFT_PRIMER': 0,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_PICK_ANYWAY': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 15,
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
     args.output.write('LEFT FLANK PRIMERS\n***********************\n')
     write_primers(lflanklprimer, lflankrprimer, args.output)
else:
     print('LEFT FLANK PRIMERS\n***********************\n')
     print_primers(lflanklprimer, lflankrprimer)

rflanklprimer = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'Right Flank',
        'SEQUENCE_TEMPLATE': seqrec,
        'SEQUENCE_INCLUDED_REGION': [args.rightbound, args.rightflanksize],
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_TASK': 'pick_cloning_primers',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 0,
        'PRIMER_PICK_ANYWAY': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 15,
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

rflankrprimer = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'Right Flank',
        'SEQUENCE_TEMPLATE': seqrec,
        'SEQUENCE_INCLUDED_REGION': [args.rightbound, args.rightflanksize],
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 0,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_PICK_ANYWAY': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 15,
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
     args.output.write('\n\nRIGHT FLANK PRIMERS\n***********************\n')
     write_primers(lflanklprimer, lflankrprimer, args.output)
else:
     print('\n\nRIGHT FLANK PRIMERS\n***********************\n')
     print_primers(lflanklprimer, lflankrprimer)
     

