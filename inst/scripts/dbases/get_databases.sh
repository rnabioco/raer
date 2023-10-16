#! /usr/bin/env bash

# download genome sequences and transcript annotations

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
gunzip GRCm38.primary_assembly.genome.fa.gz
samtools faidx GRCm38.primary_assembly.genome.fa

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
gunzip gencode.v37.annotation.gtf.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
gunzip gencode.vM25.annotation.gtf.gz

# make STAR indexes

STAR \
  --runMode genomeGenerate  \
  --runThreadN 12 \
  --genomeDir star/GRCm38 \
  --genomeFastaFiles GRCm38.primary_assembly.genome.fa

STAR \
  --runMode genomeGenerate  \
  --runThreadN 12 \
  --genomeDir star/GRCh38 \
  --genomeFastaFiles GRCh38.primary_assembly.genome.fa 
