## Examples

```
### Bowtie2
./bowtie2-build shared/hxb2.fasta shared/hxb2
./bowtie2 -x shared/hxb2 -1 ../FASTQ_File/Raw/Run04/10_S10_L001_R1_001.fastq.gz -2 ../FASTQ_File/Raw/Run04/10_S10_L001_R2_001.fastq.gz -S 10_S10_L001_R1_001.bowtie2.sam

### BWA MEM
./bwa index shared/hxb2.fasta
./bwa mem shared/hxb2.fasta ../FASTQ_File/Raw/Run04/10_S10_L001_R1_001.fastq.gz ../FASTQ_File/Raw/Run04/10_S10_L001_R2_001.fastq.gz -o 10_S10_L001_R1_001.bwamem.bam

# BAM to Consensus: https://samtools.github.io/bcftools/howtos/consensus-sequence.html

### ShoRAH
shorah -b 10_S10_L001_R1_001.bowtie2.sorted.bam -f shared/hxb2.fasta -r K03455.1:455-810 amplicon -m 0.5

### InDelFixer (Broken)
./indelfixer -i ../FASTQ_File/Raw/Run04/10_S10_L001_R1_001.fastq.gz -ir ../FASTQ_File/Raw/Run04/10_S10_L001_R2_001.fastq.gz -g shared/hxb2.fasta -illumina

### Tanoti (Broken)
./tanoti -r shared/hxb2.fasta -i ../FASTQ_File/Raw/Run04/10_S10_L001_R1_001.fastq.gz ../FASTQ_File/Raw/Run04/10_S10_L001_R2_001.fastq.gz -p 1 -o 10_S10_L001_R1_001.tanoti.sam
```
