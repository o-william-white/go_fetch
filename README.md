# go fetch

Simple python script to count or fetch reference sequences from NCBI and format for GetOrganelle (i.e. seed and gene format). 

The script requires a single taxonomic rank (--taxonomy) or taxonomic lineage (--lineage) as input.

If the --download option is specified, the user must also specify --min and --max thresholds for the number of sequnces to download. 

If the --download option is used with --taxonomy, the script will download sequences for this rank. 

If the --download option is used with --lineage, the script will iterate through the lineage and find the rank at which a minimum number of references is available to download. 

Sequences are downloaded from NCBI using using Biopython. 

If the target sequence downloaded is an organelle genome, tandem repeats are masked using trf and annotated genes are extrated using GetOrganelle util script get_annotated_regions_from_gb.py. 

Please report any issues to o.william.white@gmail.com

### Dependencies

The script requires getorganelle, biopython and trf. These can be installed in a conda environment:
```
conda create -n go_fetch -c bioconda getorganelle biopython trf
```

### Usage

Below are simple examples of how to use the script. 
```
# count sequences chloroplast sequences for Arabidopsis
python go_fetch.py --taxonomy "Arabidopsis" --target chloroplast --email o.william.white@gmail.com

# count sequences chloroplast sequences across lineage for Arabidopsis thaliana
python go_fetch.py --lineage "Brassicales,Brassicaceae,Camelineae,Arabidopsis,Arabidopsis thaliana" --target chloroplast --email o.william.white@gmail.com

# download chloroplast sequences for Arabidopsis
python go_fetch.py --taxonomy Arabidopsis --target chloroplast --download --min 2 --max 10 --output example --name Arabidopsis_go_db --email o.william.white@gmail.com

# download chloroplast sequences from lineage for Arabidopsis thaliana
python go_fetch.py --lineage "Brassicales,Brassicaceae,Camelineae,Arabidopsis,Arabidopsis thaliana" --target chloroplast --download --min 2 --max 10 --output example --name Arabidopsis_go_db --email o.william.white@gmail.com --overwrite
```
