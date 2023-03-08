import os
import sys
import shutil
import argparse
from Bio import Entrez, SeqIO
import subprocess

# argparse
argparse_description = '''
Simple python script to count or fetch reference sequences from NCBI and format for GetOrganelle (i.e. seed and gene format).

The script requires a single taxonomic rank (--taxonomy) or taxonomic lineage (--lineage) as input.

If the --download option is specified, the user must also specify --min and --max thresholds for the number of sequnces to download.

If the --download option is used with --taxonomy, the script will download sequences for this rank.

If the --download option is used with --lineage, the script will iterate through the lineage and find the rank at which a minimum number of references is available to download.

Sequences are downloaded from NCBI using using Biopython.

If the target sequence is an organelle genome, tandem repeats are masked using trf and annotated genes are extrated using GetOrganelle util script get_annotated_regions_from_gb.py.

Please report any issues to o.william.white@gmail.com
'''
argparse_usage = '''

# count sequences chloroplast sequences for Arabidopsis
python go_fetch.py --taxonomy "Arabidopsis" --target chloroplast --email o.william.white@gmail.com

# count sequences chloroplast sequences across lineage for Arabidopsis thaliana
python go_fetch.py --lineage "Brassicales,Brassicaceae,Camelineae,Arabidopsis,Arabidopsis thaliana" --target chloroplast --email o.william.white@gmail.com

# download chloroplast sequences for Arabidopsis
python go_fetch.py --taxonomy Arabidopsis --target chloroplast --download --min 2 --max 10 --output example --name Arabidopsis_go_db --email o.william.white@gmail.com

# download chloroplast sequences from lineage for Arabidopsis thaliana
python go_fetch.py --lineage "Brassicales,Brassicaceae,Camelineae,Arabidopsis,Arabidopsis thaliana" --target chloroplast --download --min 2 --max 10 --output example --name Arabidopsis_go_db --email o.william.white@gmail.com --overwrite

'''

# argparse
parser = argparse.ArgumentParser(description=argparse_description, usage=argparse_usage)
parser.add_argument('--taxonomy',     help='Taxonomy of rank to search for e.g. "Arabidopsis"', type=str, required=False)
parser.add_argument('--lineage',      help='Comma separated list of lineage to search across e.g. "Brassicales,Brassicaceae,Camelineae,Arabidopsis,Arabidopsis thaliana"', type=str, required=False)
parser.add_argument('--target',       help='Target sequence type.', choices=['chloroplast', 'mitochondrion', 'ribosomal'], required=True)
parser.add_argument('--download',     help='Only count the number of sequences on NCBI', action='store_true', required=False)
parser.add_argument('--name',         help='Database name. Required with --download.', required=False)
parser.add_argument('--min',          help='Minimum number of target sequences to download. Required with --download.', type=int, required=False)
parser.add_argument('--max',          help='Maximum number of target sequences to download. Must be larger than --min. Required with --download.',       type=int, required=False)
parser.add_argument('--output',       help='Output directory. Required with --download.', required=False)
parser.add_argument('--email',        help='Email for Entrez.', required=True)
parser.add_argument('--overwrite',    help='Overwrite output directory.', action='store_true', required=False)

args = parser.parse_args()

# additional argparse checks
if args.taxonomy != None and args.lineage != None:
    sys.exit('Error: Specify --taxonomy or --lineage, not both')

if args.taxonomy == None and args.lineage == None:
    sys.exit('Error: Must specify --taxonomy or --lineage')

if args.download and (args.name == None or args.output == None or args.min == None or args.min == None):
    sys.exit('Error: Must specify --name, --output, --min and --max when using --download')

if args.download and args.max <= args.min:
    sys.exit('Error: --max must be larger than --min')

def is_lineage(lineage):
    if len(lineage.split(",")) > 1:
        return True
    else:
        return False

if args.taxonomy != None and is_lineage(args.taxonomy):
    sys.exit('Error: --taxonomy looks like a lineage')

if args.lineage != None and is_lineage(args.lineage) == False:
    sys.exit('Error: --lineage looks like a single taxonomic rank')

# set email
Entrez.email = args.email

# functions

# create dir and overwrite if specified
def create_dir(path, target, overwrite):
    if os.path.exists(path):
        if overwrite == True:
            shutil.rmtree(path)
        else:
            sys.exit(f'{path} already exists. Remove or use --overwrite')
    os.mkdir(path)
    os.mkdir(f'{path}/fasta')
    if target != 'ribosomal':
        os.mkdir(f'{path}/genbank')

# check if taxonomy exists on ncbi
def taxonomy_exists(taxonomy):
    if taxonomy == 'NA':
        return False
    else:
        handle = Entrez.esearch(db='Taxonomy', term=taxonomy + '[Scientific Name]')
        record = Entrez.read(handle)
        count = int(record['Count'])
        if count != 0:
            return True
        else:
            return False

# esearch
def target_esearch(rank, target, max_n):
    # define search term
    if target == 'chloroplast' or target == 'mitochondrion':
        search_term = rank + '[Organism] AND ' + target + '[Title] AND complete genome[Title] AND refseq[filter]'
    else:
        if target == 'ribosomal':
            search_term = rank + '[Organism] AND 28S[Title] AND refseq[filter]'
    # esearch
    handle = Entrez.esearch(db='Nucleotide', term=search_term, retmax=max_n)
    record = Entrez.read(handle)
    return record

# efetch
def target_efetch(id, format, output_dir):
    handle = Entrez.efetch(db='Nucleotide', id=id, rettype=format, retmode='text')
    seq_record = SeqIO.read(handle, format)
    output_path = f'{output_dir}/{seq_record.id}.{format}'
    print(f'    Downloading {seq_record.id}.{format}')
    SeqIO.write(seq_record, output_path, format)

# print counts for taxonomy and return record ids
def taxonomy_count(taxonomy, target):
    # check if taxonomy exists in ncbi
    if not taxonomy_exists(taxonomy):
        # print counts
        print(f'{taxonomy}\tNA')
        # return None
        return None
    else:
        record = target_esearch(taxonomy, target, 1000000)
        counts = record["Count"]
        # print counts
        print(f'{taxonomy}\t{counts}')
        # return list of unique record ids
        target_ids = []
        for record_id in record['IdList']:
            if record_id not in target_ids:
                target_ids.append(record_id)
        return target_ids

# download taxonomy
def taxonomy_download(taxonomy, target, min_n, max_n, path):
    target_ids = taxonomy_count(taxonomy, target)
    if target_ids == None:
        sys.exit(f'Error: No {target} sequence files available for {taxonomy}')
    elif len(target_ids) < min_n:
        sys.exit(f'Error: Only {len(target_ids)} {target} sequence files available for {taxonomy}. Does not meet minimum theshold of {min_n}')
    elif len(target_ids) >= min_n:
        print(f'\nMinimum threshold of {min_n} met for {taxonomy}')
        # print the number of sequences that will be downloaded
        if len(target_ids) > max_n:
            print(f'\nDownloading the first {max_n} sequences')
            target_ids = target_ids[0:max_n]
        else:
            print(f'\nDownloading {len(target_ids)} sequences')
        # download sequences
        for target_id in target_ids:
            target_efetch(id=target_id, format='fasta', output_dir=f'{path}/fasta')
            if target != 'ribosomal':
                target_efetch(id=target_id, format='gb', output_dir=f'{path}/genbank')

# download lineage
def lineage_download(lineage, target, min_n, max_n, path):
    threshold_met = False
    target_ids_lineage = []
    lineage_root = lineage[-1]
    for taxonomy in lineage:
        if threshold_met == False:
            target_ids = taxonomy_count(taxonomy, target)
            if target_ids != None:
                target_ids_lineage.extend(target_ids)
                # catch error if root of lineage reached abd there are not enough records on ncbi
                if taxonomy == lineage_root and len(target_ids) < min_n:
                    sys.exit(f'Error: Only {len(target_ids)} {target} sequence files available for {taxonomy} at root of lineage. Does not meet minimum theshold of {min_n}')
                if len(target_ids) >= min_n:
                    print(f'\nMinimum threshold of {min_n} met for {taxonomy}')
                    print(target_ids)
                    threshold_met = True
                    # print the number of sequences that will be downloaded
                    if len(target_ids) > max_n:
                        print(f'\nDownloading the first {max_n} sequences')
                        target_ids = target_ids[0:max_n]
                        print(target_ids)
                    else:
                        print(f'\nDownloading {len(target_ids)} sequences')
                    # download sequences
                    for target_id in target_ids:
                        target_efetch(id=target_id, format='fasta', output_dir=f'{path}/fasta')
                        if target != 'ribosomal':
                            target_efetch(id=target_id, format='gb', output_dir=f'{path}/genbank')

# cat files
def cat_files(input_dir, file_ending, output_file):
    # cmd to concatenate files
    cmd_cat = ['cat']
    for f in os.listdir(input_dir):
        if f.endswith(file_ending):
             cmd_cat.append(f'{input_dir}/{f}')
    # subprocess run
    result_cat = subprocess.run(cmd_cat, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # capture stdout to text file
    with open(output_file, 'w') as log:
        log.write(f'{result_cat.stdout}\n')

# format seed
def format_seed(target, path):
    print('\nFormating seed database')
    if target == 'chloroplast' or target == 'mitochondrion':
        print('    Running trf')
        # get pwd
        pwd = os.getcwd()
        # cd to dir containing fasta files
        os.chdir(f'{path}/fasta')
        # iterate through fasta files
        for fasta in os.listdir():
            if fasta.endswith('.fasta'):
                # define cmd
                cmd_trf = f'trf {fasta} 2 7 7 80 10 50 500 -f -d -m -h'
                print(f'    {cmd_trf}')
                # subprocess run
                result_trf = subprocess.run(cmd_trf.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                # capture stdout and stderr to text file
                with open(f'{fasta}.log', 'w') as log:
                    log.write(f'stdout:\n{result_trf.stdout}\n')
                    log.write(f'stderr:\n{result_trf.stderr}\n')
        # change back to orginal dir
        os.chdir(pwd)
        # cat masked fastas
        cat_files(f'{path}/fasta', '.mask',  f'{path}/seed.fasta')
    else:
        cat_files(f'{path}/fasta', '.fasta', f'{path}/seed.fasta')

# format gene
def format_gene(target, path):
    print('\nFormating gene database')
    if target == 'chloroplast' or target == 'mitochondrion':
        print('    Running get_annotated_regions_from_gb.py')
        cmd_gar = ['get_annotated_regions_from_gb.py']
        for gb in os.listdir(f'{path}/genbank/'):
            if gb.endswith('.gb'):
                cmd_gar.append(f'{path}/genbank/{gb}')
        cmd_gar.extend(['-o', f'{path}/annotated_regions', '-t', 'CDS', '--mix'])
        # subprocess run
        result_gar = subprocess.run(cmd_gar, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # capture stdout and stderr to text file
        with open(f'{path}/annotated_regions/get_annotated_regions_from_gb.log', 'w') as log:
            log.write(f'stdout:\n{result_gar.stdout}\n')
            log.write(f'stderr:\n{result_gar.stderr}\n')
        # cp gene file to sample dir
        subprocess.run(['cp', f'{path}/annotated_regions/gene/gene.fasta', f'{path}/'])
    else:
        # cp seed file and adjust sample names
        subprocess.run(['cp', f'{path}/seed.fasta', f'{path}/gene.fasta'])
        # remove spaces
        subprocess.run(['sed', '-i', '-e', 's: :_:g',       f'{path}/gene.fasta'])
        # adjust sample names
        subprocess.run(['sed', '-i', '-e', 's:>:>28S - :g', f'{path}/gene.fasta'])

# main

# create dir if download true
if args.download:
    # make main output dir if not already present
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    # define output path
    output_path = f'{args.output}/{args.name}/'
    # create directory for results
    create_dir(path=output_path, target=args.target, overwrite=args.overwrite)

# taxonomy count
if args.taxonomy != None and args.download == False:
    print(f'\nCounting the number of {args.target} sequences for {args.taxonomy}:\n')
    taxonomy_count(args.taxonomy, args.target)

# taxonomy download
if args.taxonomy != None and args.download == True:
    print(f'\nDownloading {args.target} sequences for {args.taxonomy}:\n')
    taxonomy_download(args.taxonomy, args.target, args.min, args.max, output_path)
    # format seed
    format_seed(args.target, output_path)
    # format gene
    format_gene(args.target, output_path)

# lineage
if args.lineage != None:

    # format lineage as list
    lineage = args.lineage.split(',')
    # get root
    lineage_root = lineage[0]
    # reverse so lowest rank first
    lineage.reverse()

    # lineage count
    if args.download == False:
        print(f'\nCounting the number of {args.target} sequences across the lineage:')
        print(f'   {";".join(lineage)}\n')
        for taxonomy in lineage:
            taxonomy_count(taxonomy, args.target)

    # lineage download
    if args.download:
        print(f'\nDownloading {args.target} sequences across the lineage:')
        print(f'   {";".join(lineage)}\n')
        lineage_download(lineage, args.target, args.min, args.max, output_path)
        # format seed
        format_seed(args.target, output_path)
        # format gene
        format_gene(args.target, output_path)

# Complete :)
print('\ngo_fetch.py finished!\n')
