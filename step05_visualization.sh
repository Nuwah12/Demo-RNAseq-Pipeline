###############
# RNA-seq Demo Pipeline
# Step 5: Visualization
###############
if [[ $# -ne 4 ]]; then
        echo "Usage: ./step04_quantification <bam> <ChromSizes> <outdir> <cpu>"
        exit 1
fi

bam=${1}
chromsizes=${2}
outdir=${3}
cpu=${4}

prefix=$(basename ${bam} .final.bam)

READS=$(samtools view -c ${bam})
FACTOR=$(echo "scale=10; 1000000/${READS}" | bc -l)
bamToBed -i ${bam} -bed12 | bed12ToBed6 -i stdin | genomeCoverageBed -bg -i stdin -g ${chromsizes} -scale \${FACTOR} | sort -k1,1 -k2,2n > ${outdir}/${prefix}.tmp.bg
bedGraphToBigWig ${outdir}/${prefix}.tmp.bg ${chromsizes} ${outdir}/${prefix}.bw
rm ${outdir}/${prefix}.tmp.bg
