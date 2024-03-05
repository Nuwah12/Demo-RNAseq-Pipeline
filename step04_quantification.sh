###############
# RNA-seq Demo Pipeline
# Step 4: Quantification
###############
if [[ $# -ne 4 ]]; then
        echo "Usage: ./step04_quantification <bam> <GTF> <outdir> <cpu>"
        exit 1
fi

bam=${1}
gtf=${2}
outdir=${3}
cpu=${4}

prefix=$(basename ${bam} .final.bam)

featureCounts \
	-p \
	-t exon \
	-g gene_id \
	-s 1 \
	-O \
	-T ${cpu} \
	-a ${gtf} \
	-o ${outdir}/${prefix}.counts \
	${bam}

##### RPM Normalization
if [[ -f ${outdir}/${prefix}.counts ]]; then
	rm ${outdir}/${prefix}.counts.normalized
fi
touch ${outdir}/${prefix}.counts.normalized
echo -e 'Gene\tcounts\tRPM\tRPKM' >> ${outdir}/${prefix}.counts.normalized
counts=`awk 'BEGIN {reads=0} {reads+=$7} END {print reads}' ${outdir}/${prefix}.counts`
echo ${counts} counts mapped to genome.
awk -v counts=${counts} 'BEGIN{OFS="\t"} NR>2 {rpm=($7/counts)*1000000; rpkm=(($7*1000000000)/(counts*$6)); print $1, $7, rpm, rpkm}' ${outdir}/${prefix}.counts >> ${outdir}/${prefix}.counts.normalized

