# RNAseq Pipeline

## Data Acquistion 
Obtain FASTQ files for data on Primary human airway smooth muscle (ASM) cell lines untreated or treated with dexamethasone (Dex in file names).
Data obtained from Himes et al., PLOS one, 2014. https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625#s2

```

```

Obtain a reference genome from ENSEMBL. Use the soft-masked primary assembly file for the full reference genome and for chromosome 22 (.dna_sm). Also pull a reference annotation (GTF file)

```
#Create a new directory and copy the two reference files into it.
mkdir /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/RNA_REF_DIR
cd RNA_REF_DIR
wget
wget
wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
gunzip Homo_sapiens.GRCh38.86.gtf.gz
```
## Run FASTQC
Check quality of the data with FASTQC 

```
#Reserve one compute node (-N1) for an hour (-t1:00:00) and specify --ntasks=8
salloc -N1 -t1:00:00 --ntasks=8

#Load Java then FASTQC
module load java fastqc

#Navigate to the directory housing the FASTQ files
cd /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA

#Run FASTQC on all fastq files
fastqc *.fastq.gz
```
#OPTIONAL Review FASTQC output by transfering the .html files (code to be done on computer directory NOT server)

```
#scp command follows format of USERNAME@teach.scinet.utoronto.ca:path_to_your_FASTQC_files
scp lcl_uotmmg3003s2058@teach.scinet.utoronto.ca:/scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/\*.html ./
```

## Alignment with HISAT2
Create or download a version of the human genome which has been indexed for your particular aligner (HISAT2 in this case). HISAT2 in particular already has multiple versions available at https://registry.opendata.aws/jhu-indexes
```
#Create a new directory in the scratch directory for the HISAT2 index
mkdir /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/RNA_REF_INDEX

#A HISAT2 index is already available in another account's directory on the server. Thus I will be transferring the index from there.
cd /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/RNA_REF_INDEX
cp /home/l/lcl_uotmmg3003/mmg3003starter/RNA_REF_INDEX/grch38_snptran.tar.gz ./
tar -xvf grch38_snptran.tar.gz

```
Align the quality checked data for each sample (the FASTQ files) to the reference human genome

```
#Create a new directory in the scratch directory
mkdir /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_HISAT2_ALIGNMENTS

#
```
