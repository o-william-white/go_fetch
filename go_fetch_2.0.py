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

# get scientific name from taxonomic id
def get_scientific_name(taxid):
    handle = Entrez.efetch(db="Taxonomy", id=taxid)
    record = Entrez.read(handle)
    return record[0]["ScientificName"]

assert get_scientific_name("3702") == "Arabidopsis thaliana"

# get lineage from taxid
def get_lineage(taxid):
    # efetch
    handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(handle)
    # get lineage
    lineage = records[0]["Lineage"].split("; ")[::-1]
    # return lineage
    return lineage

assert get_lineage(3701) == ['Camelineae', 'Brassicaceae', 'Brassicales', 'malvids', 'rosids', 'Pentapetalae', 'Gunneridae', 'eudicotyledons', 'Mesangiospermae', 'Magnoliopsida', 'Spermatophyta', 'Euphyllophyta', 'Tracheophyta', 'Embryophyta', 'Streptophytina', 'Streptophyta', 'Viridiplantae', 'Eukaryota', 'cellular organisms']

# generate search term
# target = 'chloroplast', 'mitochondrion', 'robosomal'
# db can be refseq or genbank
def search_term(taxid, target, db):
    # add taxid to term
    term = f"{taxid}[Organism]"
    if target == 'chloroplast' or target == 'mitochondrion':
        # genbank is part of the International Nucleotide Sequence Database Collaboration (INSDC) along with the European Nucleotide Archive and the DNA Data Bank of Japan (DDBJ)
        # refseq are derived from genbank but not part of
        term += f" AND {target}[Title] AND complete genome[Title]"
    if target == 'robosomal':
        term += f" AND 28S[Title] AND 18S[Title] AND 5.88S[Title]"
    if db == 'refseq':
        term += f" AND refseq[filter]"
    if db == 'genbank':
        term += f" AND ddbj_embl_genbank[filter]"
    return term 

print(search_term(get_scientific_name(3701), "chloroplast", "refseq"))


