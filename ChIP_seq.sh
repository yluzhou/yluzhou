# ChIP-seq analysis pipeline

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -c hg38.fa.gz > hg38.fa
bowtie2-build --threads 10 -f hg38.fa hg38
index="hg38"

trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 -o  trim/ ${raw_fq}
clean_fq="trim/${sample}_trimmed.fq.gz"

# alignment
bowtie2 -p 5 -x ${index} -U ${clean_fq} | samtools sort -O bam -@ 8 -o ${sample}.bam
samtools index ${sample}.bam

# remove PCR duplicate
sambamba markdup -r -p -t 8 ${sample}.bam ${sample}.rmdup.bam
samtools index ${sample}.rmdup.bam
samtools idxstats ${sample}.rmdup.bam | cut -f 1 | grep -v -E 'chrM|chrMT|MT' | xargs samtools view -b -f 2 -q 30 -@ 5 ${sample}.rmdup.bam > ${sample}.rmdup.rmchrM.bam
samtools index ${sample}.rmdup.rmchrM.bam
bedtools bamtobed -i ${sample}.rmdup.rmchrM.bam  > ${sample}.bed

# filter blacklist
bedtools intersect -v -a ${sample}.rmdup.rmchrM.bam -b hg38.blacklist.bed > ${sample}.blacklist_filtered.last.bam
samtools index ${sample}.blacklist_filtered.last.bam
bedtools bamtobed -i ${sample}.blacklist_filtered.last.bam > ${sample}.last.bed

# call peaks
macs3 callpeak -t ChIP.bam -c Control.bam -f BAM -n ${name} -g hs --outdir peaks/ -q 0.05 --keep-dup all

# generate bigwig file for visualization
bamCoverage --normalizeUsing RPKM -b ${sample}.blacklist_filtered.last.bam -o ${sample}.chip.bw
