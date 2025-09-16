#!/bin/bash

mkdir -p $4
while read line; do
  hisat2 -p 60 -x $2 -U $1/${line}.merge.fq.gz -S $4/${line}.merge.sam
  samtools view -S -b $4/${line}.merge.sam -o $4/${line}.bam
  samtools sort $4/${line}.bam -o $4/${line}.sort.bam
  samtools rmdup -s $4/${line}.sort.bam $4/${line}.sort.rmdup.bam
  samtools index $4/${line}.sort.rmdup.bam
  rm $4/*.sam
  rm $4/*_*.bam
  rm $4/${line}.bam
  rm $4/${line}.sort.bam
done < $3