# make small test data bam files
# bam files from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99249
# were processed with editing pipeline from vicens project to generate bam
# files
# see /beevol/home/riemondy/dev/ullr/test-data for bam files
fasta=/beevol/home/rbilab/ref/genome/human/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.UCSC.fa

# first make a tiny (< 100k) fasta file reference
bedtools getfasta -fi $fasta -bed test_regions.bed  -name > human.fasta

## build star index
STAR --runMode genomeGenerate --genomeDir star/ --genomeFastaFiles human.fasta --genomeSAindexNbases 4

## extract reads in fasta file regions
samtools view  -M -L test_regions.bed SRR5564269_dedup_sorted.bam \
    | cut -f 1 \
    | uniq > reads_to_get.txt

samtools view -b -N reads_to_get.txt SRR5564269_dedup_sorted.bam \
    | samtools sort -n \
    | samtools fastq -1 SRR5564269_1.fastq.gz -2 SRR5564269_2.fastq.gz

samtools view  -M -L test_regions.bed SRR5564277_dedup_sorted.bam \
    | cut -f 1 \
    | uniq > reads_to_get.txt

samtools view -b -N reads_to_get.txt SRR5564277_dedup_sorted.bam \
    | samtools sort -n \
    | samtools fastq -1 SRR5564277_1.fastq.gz -2 SRR5564277_2.fastq.gz

## remap

STAR \
    --readMapNumber 500 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMmode Full \
    --outFileNamePrefix SRR5564277/SRR5564277_ \
    --readFilesCommand "gunzip -c" \
    --genomeDir star/ \
    --readFilesIn SRR5564277_1.fastq.gz SRR5564277_2.fastq.gz

samtools index SRR5564277/SRR5564277_Aligned.sortedByCoord.out.bam

STAR \
    --readMapNumber 500 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMmode Full \
    --outFileNamePrefix SRR5564269/SRR5564269_ \
    --readFilesCommand "gunzip -c" \
    --genomeDir star/ \
    --readFilesIn SRR5564269_1.fastq.gz SRR5564269_2.fastq.gz

samtools calmd -b SRR5564269_Aligned.sortedByCoord.out.bam human.fasta > SRR5564277_Aligned.sortedByCoord.out.md.bam
samtools calmd -b SRR5564277_Aligned.sortedByCoord.out.bam human.fasta > SRR5564277_Aligned.sortedByCoord.out.md.bam

samtools index SRR5564269/SRR5564269_Aligned.sortedByCoord.out.md.bam
samtools index SRR5564277/SRR5564277_Aligned.sortedByCoord.out.md.bam
