# ChIP-seq analysis pipeline

trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 -o  trim/ ${fq1}

# alignment
bowtie2 -p 5 -x ${index} -U $id | samtools sort -O bam -@ 5 -o ./${sample}.bam

# remove PCR duplicate
sambamba markdup -r -p -t 12 $id $(basename -s .bam $id).rmdup.bam

# filter blacklist
bedtools intersect -v -a ${sample}.rmdup.rmchrM.bam -b hg38.blacklist.bed > ${sample}.blacklist_filtered.last.bam
samtools index ${sample}.blacklist_filtered.last.bam
bedtools bamtobed -i ${sample}.blacklist_filtered.last.bam > ${sample}.last.bed

# call peaks
macs3 callpeak -t ChIP.bam -c Control.bam -f BAMPE -n ${name} -g hs --outdir peaks/ -q 0.05 --keep-dup all

