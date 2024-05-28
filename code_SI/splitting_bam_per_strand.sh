

inp=$1
##annot=$2
name=$(basename ${inp})
out_name=${name%%.*}

fwd1=${out_name}_fwd1.bam
fwd2=${out_name}_fwd2.bam
rev1=${out_name}_rev1.bam
rev2=${out_name}_rev2.bam

samtools view -b -f 128 -F 16 ${inp} > ${fwd1}
samtools index ${fwd1}

samtools view -b -f 64 -F 32 ${inp} > ${fwd2}
samtools index ${fwd2}

samtools merge -f ${out_name}_fwd.bam ${fwd1} ${fwd2} 
samtools sort ${out_name}_fwd.bam >  ${out_name}_fwd_sorted.bam
samtools index ${out_name}_fwd_sorted.bam
rm -f ${fwd1} ${fwd2} ${fwd1}.bai ${fwd2}.bai  ${out_name}_fwd.bam


samtools view -b -f 144 ${inp} > ${rev1}
samtools index ${rev1}

samtools view -b -f 96 ${inp} > ${rev2}
samtools index ${rev2}

samtools merge -f ${out_name}_rev.bam ${rev1} ${rev2}
samtools sort ${out_name}_rev.bam >  ${out_name}_rev_sorted.bam
samtools index ${out_name}_rev_sorted.bam
rm -f ${rev1} ${rev2} ${rev1}.bai ${rev2}.bai  ${out_name}_rev.bam 

