import sys
from parser import parser

if __name__ == '__main__':
    p = parser()
    p.parse_fasta_file(filename='TestFiles\\reads.txt')