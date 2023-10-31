# Go fetch

Simple python script to fetch organelle or ribosomal reference sequences from NCBI for a given taxonomy. 

The script was originally designed to format reference data for GetOrganelle (i.e. seed and gene format). You can specify if you want sequences formatted for GetOrganelle using the `--getorganelle` parameter.

### Dependencies

The script requires getorganelle, biopython and trf. These can be installed in a conda environment:
```
conda create -n go_fetch -c bioconda getorganelle biopython trf
```

### Input

Below is a summary of the main input parameters.

```
  -h, --help            show this help message and exit
  --taxonomy TAXONOMY   Taxonomy of rank to search for e.g. "Arabidopsis"
  --target TARGET       Target sequence type. Options = [chloroplast, mitochondrion, ribosomal].
  --db DB               Database to search. Options = [refseq, genbank].
  --min MIN             Minimum number of target sequences to download.
  --max MAX             Maximum number of target sequences to download. 
  --seed SEED           Seed used for subsampling.
  --output OUTPUT       Output directory.
  --overwrite           Overwrite output directory.
  --getorganelle        Format seed and gene database for get organelle.
  --email EMAIL         Email for Entrez.
  --api API             API key for NCBI. Optional.
```

### How it works

The user provides a taxonomy as a taxonomic name (e.g. "Arabidopsis thaliana") or NCBI taxonomy ID (e.g. "3702"). However, **it is reccomended to use NCBI taxonomy IDs**, as some unrelated groups share names, for example the name Drosophila is used for both a genus of fruit flies and basidiomycete fungi. 

The script will count the number or target sequences available for the user given taxonomy. If there are not enough sequence available based on the `--min` parameter, the script will count the number of sequences available for the parent rank, until the minimum threshold is reached. 

If the maximum number of sequences is exceeded, based on the `--max` parameter, the script will subsample the children lineages, aiming for an even sampling across children. 

### Example usage

Below are simple examples of how to use the script. 
```
# Arabidopsis chloroplast
python3 go_fetch.py \
   --taxonomy 3701 \
   --target mitochondrion 
   --db genbank z
   --min 5 --max 10 \
   --output arabidopsis_chloroplast \
   --overwrite \
   --getorganelle \
   --email user_email@example.com

# Arabidopsis ribosomal
python3 go_fetch.py \
   --taxonomy "Arabidopsis" \
   --target ribosomal \
   --db genbank \
   --min 5 --max 10 \
   --output arabidopsis_ribosomal \
   --overwrite \
   --getorganelle \
   --email user_email@example.com

# Drosophila mitochondrion
python3 go_fetch.py \
   --taxonomy 7215 \
   --target mitochondrion \
   --db genbank \
   --min 10 --max 50 \
   --output drosophila_mitochondrion \
   --overwrite \
   --getorganelle \
   --email user_email@example.com

# Drosophila ribosomal
python3 go_fetch.py \
   --taxonomy 7215 \
   --target ribosomal \
   --db genbank \
   --min 10 --max 50 \
   --output drosophila_ribosomal \
   --overwrite \
   --getorganelle \
   --email user_email@example.com

```

