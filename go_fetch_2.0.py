import os
import sys
import shutil
import argparse
from Bio import Entrez, SeqIO
import subprocess

# argparse
argparse_description = '''

'''

argparse_usage = '''

'''

# argparse
parser = argparse.ArgumentParser(description=argparse_description, usage=argparse_usage)
parser.add_argument('--taxonomy',     help='Taxonomy of rank to search for e.g. "Arabidopsis"', type=str, required=True)
#parser.add_argument('--target',       help='Target sequence type.', choices=['chloroplast', 'mitochondrion', 'ribosomal'], required=True)
#parser.add_argument('--db',           help='Databse to search/download from. Either refseq (NCBI) or genbank (INSDC). Default=refseq.', choices=['refseq', 'genbank'], required=False, default='refseq')
#parser.add_argument('--name',         help='Database name. Required with --download.', required=False)
#parser.add_argument('--min',          help='Minimum number of target sequences to download. Required with --download.', type=int, required=False)
#parser.add_argument('--max',          help='Maximum number of target sequences to download. Must be larger than --min. Required with --download.',       type=int, required=False)
#parser.add_argument('--output',       help='Output directory. Required with --download.', required=False)
parser.add_argument('--email',        help='Email for Entrez.', required=True)
#parser.add_argument('--overwrite',    help='Overwrite output directory.', action='store_true', required=False)

args = parser.parse_args()

# additional dependency checks
#try:
#    subprocess.call(["samtools"], stderr=subprocess.DEVNULL)
#except FileNotFoundError:
#    sys.exit("Error: samtools not in path")


### set email
Entrez.email = args.email

### functions

# check if taxonomy exists
def taxid_exists(taxid):
    handle = Entrez.efetch(db="taxonomy", id=taxid)
    record = Entrez.read(handle)
    if record:
        return True
    else:
        return False
assert taxid_exists(3701) == True
assert taxid_exists(123456789) == False

# 

