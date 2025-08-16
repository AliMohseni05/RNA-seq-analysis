# RNA-seq 
## Arabidopsis thaliana
1. Install Tools (hisat2, sratoolkit, subread)
   │
2. Download SRA Data (SRR31323940, SRR31323925) → Split into FASTQ
   │
3. Download Reference Genome (FASTA) → Unzip → Build HISAT2 Index
   │
4. Align FASTQ to Genome (HISAT2) → SAM → Convert to BAM (samtools) → Sort BAM
   │
5. Download Annotation (GFF3/GTF) → Unzip
   │
6. Quantify Gene Expression:
      ├─ htseq-count → Counts (count.txt)
      └─ featureCounts → Counts (feature_counts.txt) better 

The script performs RNA-seq analysis for Arabidopsis thaliana using data from two SRA accessions (SRR31323940 and SRR31323925). Key steps include:

Software Installation: Installs hisat2, sratoolkit, and subread for alignment and quantification.

Data Download: Downloads and splits SRA data into FASTQ files for the two accessions.

Reference Genome: Downloads and unzips the Arabidopsis thaliana genome (FASTA) and annotation files (GFF3/GTF).

Alignment: Builds a HISAT2 index for the genome and aligns the FASTQ files to generate SAM/BAM files.

Quantification: Uses htseq-count and featureCounts to quantify gene expression from the aligned BAM files using the annotation file.
