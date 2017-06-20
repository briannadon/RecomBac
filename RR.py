#!/usr/bin/python

from Bio import SeqIO,pairwise2
import re
from sys import argv
import argparse

script, gb, left, right = argv
left, right = int(left), int(right)
#parser = argparse.ArgumentParser(description="Placeholder. Program to define primers in bacterial recombination.")

region = SeqIO.read(gb,'genbank')

cutseq = str(region.seq)[(left-1):(right+1)]

print(cutseq)
