# Contents
- [Contents](#contents)
- [Introduciton](#introduciton)
- [repo tree](#repo-tree)
- [Data Download](#data-download)
  - [Data Merging](#data-merging)
  - [Data Clustering](#data-clustering)

# Introduciton
This Project is for EPCD (Enterobacteriaceae Protein-CDS datasets) download.

# repo tree
``` bash
.
├── finalData # the finished data 4 train
├── readme.md
└── src # some code 4 download data
```
# Data Download
Downloading Data Using Scripts:


1. Extract `TaxPhylogeny_Name.zip` and `assembly_summary_genbank.zip`.  
2. Retrieve all sequencing links corresponding to *Enterobacteriaceae*.  
``` bash
perl Fetch_SpeHttp_Based_Phylogeny.pl Enterobacteriaceae assembly_summary_genbank.txt TaxPhylogeny_Name.txt  
```    
3. Run `split.sh` to split the file into 10 parts.(just accelerate the speed)
4. Download all CDS sequences for Enterobacteriaceae in batches.
``` bash
cd XX  
perl Fetch_Cds_BaseSpeHttp.pl split_0X  
```
Wait for the download to complete.\
Process completed.

## Data Merging  

```bash
for i in {0..10}; do cat $i/split_0$i.cds >> allInOneData/merged_output.fasta; done
```

## Data Clustering
use `MMSeqs2` do this, and the finalData is the out.