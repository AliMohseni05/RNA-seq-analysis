┌─────────────────────────────────────────────────────────────────────┐
│                        RNA-Seq Analysis Workflow                    │
├─────────────────┬───────────────────┬───────────────────────────────┤
│   Input Data    │  Processing Steps │         Output Directories    │
├─────────────────┼───────────────────┼───────────────────────────────┤
│                 │                   │                               │
│  input/*.fastq  │ 1. Initial QC     │ output/1_initial_qc/          │
│                 │  (FastQC)         │   - FastQC HTML reports       │
│                 │                   │                               │
├─────────────────┼───────────────────┼───────────────────────────────┤
│                 │ 2. Trimming       │ output/2_trimmed_output/      │
│                 │  (Trimmomatic)    │   - Adapter-free FASTQs       │
│                 │                   │                               │
├─────────────────┼───────────────────┼───────────────────────────────┤
│                 │ 3. rRNA Filtering │ output/3_rRNA/                │
│ sortmerna_db/   │  (SortMeRNA)      │   ├── filtered/*.fastq        │
│  - rRNA DBs     │                   │   ├── aligned/rRNA_reads.fastq│
│  - Indexes      │                   │   └── logs/                   │
│                 │                   │                               │
├─────────────────┼───────────────────┼───────────────────────────────┤
│ genome/*.fa     │ 4. Alignment      │ output/4_aligned_sequences/   │
│ star_index/     │  (STAR)           │   ├── aligned_bam/*.bam       │
│                 │                   │   └── aligned_logs/           │
│                 │                   │                               │
├─────────────────┼───────────────────┼───────────────────────────────┤
│ annotation/*.gtf│ 5. Quantification │ output/5_final_counts/        │
│                 │  (featureCounts)  │   - Gene count matrix         │
│                 │                   │                               │
├─────────────────┼───────────────────┼───────────────────────────────┤
│                 │ 6. MultiQC        │ output/6_multiQC/            │
│                 │  (Aggregate QC)   │   - Summary of all steps      │
└─────────────────┴───────────────────┴───────────────────────────────┘


# Download 
#Download SRA file (replace SRRXXXXXXX with your SRA run accession)
prefetch SRR9879604
# Convert
#Convert to FASTQ (single-end example)
fastq-dump --split-files  SRRXXXXXXX

# QC 
fastqc /content/SRRXXXXXXX_1.fastq

# trimmomatic
trimmomatic SE -threads 4 \
              /content/SRRXXXXXXX_1.fastq /content/SRRXXXXXXX_trim.fastq \
              ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
              LEADING:20 \
              TRAILING:20 \
              SLIDINGWINDOW:4:20 \
              MINLEN:36 \
              > trimmomatic.log 2>&1
# QC 
fastqc /content/SRRXXXXXXX_1.fastq

# annotation
## genome of your plant in fastA and annotation
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/solanum_lycopersicum/dna/Solanum_lycopersicum.SL3.0.dna.toplevel.fa.gz
wget  https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/solanum_lycopersicum/Solanum_lycopersicum.SL3.0.61.chr.gff3.gz
### unzip files
gunzip /content/Solanum_lycopersicum.SL3.0.dna.toplevel.fa.gz
gunzip /content/Solanum_lycopersicum.SL3.0.61.chr.gff3.gz

## hisat2-build to make a genome index 
hisat2-build -p 8 /content/Solanum_lycopersicum.SL3.0.dna.toplevel.fa reference_index

# combined
hisat2 -x /content/reference_index/reference_index \

        -U /content/sra-tomato-4/SRR9879604_1.fastq\

        -S /content/combined_output.sam


!ls -lh /content/combined_output.sam  # Verify SAM file was created

!head /content/combined_output.sam    # Peek at alignment results


# Sort BAM by genomic coordinates,
## make bam file
samtools view -b /content/combined_output.sam > /content/combined_output.bam
## make bam sorted 
samtools sort SRR5025144.bam -o SRR5025144.sorted.bam

# count reads 
htseq-count -f bam -s no sample.sorted.bam /content/tomato.gtf > counts.txt

#troubleshooting
If htseq-count didn't read GFF3 file use GFT file 
#### dowload GTF file 
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-54/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.54.gtf.gz

#### Uncompress the file
gunzip Arabidopsis_thaliana.TAIR10.54.gtf.gz

## Convert GFF3 → GTF
!gffread Solanum_lycopersicum.SL3.0.61.chr.gff3 -T -o tomato.gtf
