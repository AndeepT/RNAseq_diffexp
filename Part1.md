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
cd /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_HISAT2_ALIGNMENTS

#Reserve server compute node
salloc -N1 -t4:00:00 --ntasks=8

#load HISAT2
module load gcc/7.3.0 hisat2

#Run HISAT2
#hisat2 -p 8 --rg-id=HCC1395_tumor_rep1 --rg SM:HCC1395_tumor_rep1 --rg PL:ILLUMINA -x/path_to_hisat2_index --rna-strandness RF -1 /path_to_read1 -2 /path_to_read2 -S ./HCC1395_tumor_rep1.sam
#'-p 8' tells HISAT2 to use eight CPUs for alignments.
#'--rg-id specifies a read group ID that is a unique identifier. In this case we are looking at HCC1395_tumor_rep1.
#'--rg SM: specifies a read group sample name. This together with rg-id will allow you to determine which reads came from which sample in the merged bam. The rg-id and SM can be different, but in our case they are the same.
#'--rg PL:ILLUMINA' specifies the read group sequencing platform at Illumina '-x /path_to_hisat2_index/genome_snp_tran' The HISAT2 index filename prefix (minus the trailing .X.ht2). REPLACE THIS WITH THE PATH TO YOUR INDEX FILENAME PREFIX (file name minus the trailing .X.ht2)
#'--rna-strandness RF' specifies strandness of RNAseq library. We will specify RF since the TruSeq strand-specific library was used to make these libraries.
#'-1 /path_to/read1.fq.gz' The read 1 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed. REPLACE THIS WITH THE PATH TO YOUR READ1
#'-2 /path_to/read2.fq.gz' The read 2 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed. REPLACE THIS WITH THE PATH TO YOUR READ2
#'-S /path_to/output.sam' The output SAM format text file of alignments. Here I am specifying to put the SAM file in the current directory I'm in.

hisat2 -p 8 --rg-id=N052611_Dex --rg SM:N052611_Dex --rg PL:ILLUMINA -x/scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/RNA_REF_INDEX/grch38_snp_tran/genome_snp_tran --rna-strandness RF -1 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N052611_Dex_r1.fastq.gz -2 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N052611_Dex_r2.fastq.gz -S ./N052611_Dex.sam

hisat2 -p 8 --rg-id=N052611_untreated --rg SM:N052611_untreated --rg PL:ILLUMINA -x/scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/RNA_REF_INDEX/grch38_snp_tran/genome_snp_tran --rna-strandness RF -1 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N052611_untreated_r1.fastq.gz -2 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N052611_untreated_r2.fastq.gz -S ./N052611_untreated.sam

hisat2 -p 8 --rg-id=N080611_Dex --rg SM:N080611_Dex --rg PL:ILLUMINA -x/scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/RNA_REF_INDEX/grch38_snp_tran/genome_snp_tran --rna-strandness RF -1 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N080611_Dex_r1.fastq.gz -2 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N080611_Dex_r2.fastq.gz -S ./N080611_Dex.sam

hisat2 -p 8 --rg-id=N080611_untreated --rg SM:N080611_untreated --rg PL:ILLUMINA -x/scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/RNA_REF_INDEX/grch38_snp_tran/genome_snp_tran --rna-strandness RF -1 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N080611_untreated_r1.fastq.gz -2 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N080611_untreated_r2.fastq.gz -S ./N080611_untreated.sam

hisat2 -p 8 --rg-id=N61311_Dex --rg SM:N61311_Dex --rg PL:ILLUMINA -x/scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/RNA_REF_INDEX/grch38_snp_tran/genome_snp_tran --rna-strandness RF -1 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N61311_Dex_r1.fastq.gz -2 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N61311_Dex_r2.fastq.gz -S ./N61311_Dex.sam

hisat2 -p 8 --rg-id=N61311_untreated --rg SM:N61311_untreated --rg PL:ILLUMINA -x/scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/RNA_REF_INDEX/grch38_snp_tran/genome_snp_tran --rna-strandness RF -1 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N61311_untreated_r1.fastq.gz -2 /scratch/l/lcl_uotmmg3003/lcl_uotmmg3003s2058/AIRWAY_DATA/N61311_untreated_r2.fastq.gz -S ./N61311_untreated.sam


```
