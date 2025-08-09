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

fastq-dump --split-files SRR9879604

# QC 
fastqc /content/SRR9879604_1.fastq
# trimmomatic
trimmomatic SE -threads 4 \
              /content/SRR9879604_1.fastq /content/SRR9879604_trim.fastq \
              ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
              LEADING:20 \
              TRAILING:20 \
              SLIDINGWINDOW:4:20 \
              MINLEN:36 \
              > trimmomatic.log 2>&1
# QC 
fastqc /content/SRR9879604_1.fastq

# annotation
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/solanum_lycopersicum/dna/Solanum_lycopersicum.SL3.0.dna.toplevel.fa.gz

wget  https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/solanum_lycopersicum/Solanum_lycopersicum.SL3.0.61.chr.gff3.gz

gunzip /content/Solanum_lycopersicum.SL3.0.dna.toplevel.fa.gz

gunzip /content/Solanum_lycopersicum.SL3.0.61.chr.gff3.gz

hisat2-build -p 8 /content/Solanum_lycopersicum.SL3.0.dna.toplevel.fa reference_index


mkdir -p reference_index

mv *.ht2 reference_index

# combined
hisat2 -x /content/reference_index/reference_index \

        -U /content/sra-tomato-4/SRR9879604_1.fastq\

        -S /content/combined_output.sam


!ls -lh /content/combined_output.sam  # Verify SAM file was created

!head /content/combined_output.sam    # Peek at alignment results


# Sort BAM by genomic coordinates
samtools view -b /content/combined_output.sam > /content/combined_output.bam

samtools sort -@ 4 -o sample.sorted.bam /content/combined_output.bam

  

# Index the sorted BAM

!samtools index sample.sorted.bam




htseq-count -f bam -s no sample.sorted.bam /content/tomato.gtf > counts.txt

# Convert GFF3 → GTF

!gffread Solanum_lycopersicum.SL3.0.61.chr.gff3 -T -o tomato.gtf

  

# Now run htseq-count with the GTF

!htseq-count \

  -f bam \

  -s no \

  -r pos \

  tomato.sorted.bam \

  tomato.gtf \

  > counts.txt
