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
python3 go_fetch_2.0.py --taxonomy "Arabidopsis lyrata" --target chloroplast --db genbank --min 2 --max 22 --email o.william.white@gmail.com
'''

# argparse
parser = argparse.ArgumentParser(description=argparse_description, usage=argparse_usage)
parser.add_argument('--taxonomy',     help='Taxonomy of rank to search for e.g. "Arabidopsis"', type=str, required=True)
parser.add_argument('--target',       help='Target sequence type.', choices=['chloroplast', 'mitochondrion', 'ribosomal'], required=True)
parser.add_argument('--db',           help='Databse to search/download from. Either refseq (NCBI) or genbank (INSDC). Default=refseq.', choices=['refseq', 'genbank'], required=False, default='refseq')
#parser.add_argument('--name',         help='Database name. Required with --download.', required=False)
parser.add_argument('--min',          help='Minimum number of target sequences to download.', type=int, required=False)
parser.add_argument('--max',          help='Maximum number of target sequences to download. Must be larger than --min.', type=int, required=False)
#parser.add_argument('--output',       help='Output directory. Required with --download.', required=False)
parser.add_argument('--email',        help='Email for Entrez.', required=True)
#parser.add_argument('--overwrite',    help='Overwrite output directory.', action='store_true', required=False)

args = parser.parse_args()

# additional dependency checks
#try:
#    subprocess.call(["samtools"], stderr=subprocess.DEVNULL)
#except FileNotFoundError:
#    sys.exit("Error: samtools not in path")

# additional checks

if args.max <= args.min:
    sys.exit('Error: --max must be larger than --min')


### set email
Entrez.email = args.email

### functions

# get taxonomic id from scientific name
def get_taxonomic_id(taxonomy):
    handle = Entrez.esearch(db="Taxonomy", term=f"{taxonomy}[Scientific Name]")
    record = Entrez.read(handle)
    return str(record["IdList"][0])
assert get_taxonomic_id("Arabidopsis thaliana") == "3702"

# check if taxonomy id exists
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
def get_scientific_name(taxonomy):
    handle = Entrez.efetch(db="Taxonomy", id=taxonomy)
    record = Entrez.read(handle)
    return record[0]["ScientificName"]
assert get_scientific_name("3702") == "Arabidopsis thaliana"

# check if taxonomy name exists
def scientific_name_exists(taxonomy):
    handle = Entrez.esearch(db="Taxonomy", term=f"{taxonomy}[Scientific Name]")
    record = Entrez.read(handle)
    if int(record["Count"]) >= 1:
        return True
    else:
        return False
assert scientific_name_exists("Arabidopsis") == True
assert scientific_name_exists("NotArabidopsis") == False

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

def get_children(taxonomy):
    handle = Entrez.esearch(db='taxonomy', term=f"{taxonomy}[next level]", retmax=9999)
    record = Entrez.read(handle)
    children = []
    for i in record["IdList"]:
        children.append(get_scientific_name(i))
    return(children)
assert get_children("Arabidopsis lyrata") == ['Arabidopsis petraea subsp. umbrosa', 'Arabidopsis petraea subsp. septentrionalis', 'Arabidopsis lyrata subsp. lyrata', 'Arabidopsis lyrata subsp. petraea']


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

# count the number of sequences on ncbi using the term generated
def entrez_esearch(input_term):
    # esearch
    handle = Entrez.esearch(db='Nucleotide', term=input_term, retmax=999)
    record = Entrez.read(handle)
    return record["IdList"]


### main

# get and check ncbi taxonomy_id and scientific_name
try:
    int(args.taxonomy)
    print("Taxonomy input is a numeric NCBI ID.")
    # check if taxonomy id exists
    if not taxid_exists(args.taxonomy):
        sys.exit(f"Taxonomy ID {args.taxonomy} does not exist. Please check for the correct taxonomy id on https://www.ncbi.nlm.nih.gov/taxonomy")
    else:
        print("Taxonomy ID exists on NCBI")
    taxonomy_id = args.taxonomy
    taxonomy_name = get_scientific_name(taxonomy_id)
except:
    print("Taxonomy input is a scientific name.")
    if not scientific_name_exists(args.taxonomy):
        sys.exit(f"Scientific name {args.taxonomy} does not exist. Please check for the correct taxonomy id on https://www.ncbi.nlm.nih.gov/taxonomy")
    taxonomy_name = args.taxonomy
    taxonomy_id = get_taxonomic_id(taxonomy_name)

print("NCBI Id: " + taxonomy_id)
print("Scientific name: " + taxonomy_name + "\n")

# get lineage
lineage = get_lineage(taxonomy_id)
lineage.insert(0, taxonomy_name)

# get children
children = get_children(taxonomy_name)

# print phylogeny
spacer = ""
for l in lineage[::-1]:
    print(f"{spacer}{l}")
    if spacer == "":
        spacer = "|_" + spacer
    else: 
        spacer = "  " + spacer
for c in children: 
    print(f"{spacer}{c}")



####### working progress... #######

def recursive_search(taxonomy, lineage, target, db, min_th, max_th):

    # define search term
    term = search_term(taxonomy, target, db)
    print(f"\nUsing search term = {term}")

    # esearch
    idlist = entrez_esearch(term)

    # count sequences
    count = len(idlist)
    print(f"{count} sequences found")

    # does the number of sequences meet the minimum threshold
    if count >= min_th:
        print("Minimum threhold reached")

        # if maximum not exceeded, download all
        if count <= max_th: 
            print("Maximum threshold not exceeded. Downloading {count} sequences")
            print(f"EFETCH")
    
        # if maximum exceeded, download subsample
        else:
            print("Maximum threshold exceeded. Subsampling {count} sequences from children")
            print(f"SUBSAMPLE")
    
    else: 
        print("Minimum threhold not reached, moving up within lineage")
        taxonomy = lineage[lineage.index(taxonomy)+1]
        print(f"Next level is {taxonomy}")


recursive_search(taxonomy_name, lineage, args.target, args.db, args.min, args.max)

################################
