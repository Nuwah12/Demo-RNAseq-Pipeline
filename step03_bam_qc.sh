###############
# Demo RNA-seq pipeline
# Step 3: BAM file quality control
###############
if [[ $# -ne 4 ]]; then
        echo "Usage: ./step03_bam_qc.sh <bam> <ChromNoPatch> <outdir> <cpu>"
        exit 1
fi

bam=${1}
ChromNoPatch=${2}
outdir=${3}
cpu=${4}

prefix=$(basename ${bam} .raw.bam)

echo "$(date): Removing non-chromosome sequences..."
samtools view -h -L ${ChromNoPatch} ${bam} | samtools sort -m 50G - -o ${outdir}/${prefix}.noPatch.bam
echo "$(date): Marking duplicates with Picard..."
picard MarkDuplicates \
	AS=true \
	M=${outdir}/${prefix}.picard.metrics \
	O=${outdir}/${prefix}.dupmark.bam \
	I=${outdir}/${prefix}.noPatch.bam \
	TMP_DIR=/mnt/data0/noah/tmp \
	REMOVE_DUPLICATES=false \
	VALIDATION_STRINGENCY=SILENT
#rm ${outdir}/${prefix}.noPatch.bam
echo $(date): Sorting duplicate-marked bam...""
samtools sort -m 50G ${outdir}/${prefix}.dupmark.bam -o ${outdir}/${prefix}.dupmark.sorted.bam
#rm ${outdir}/${prefix}.dupmark.bam
echo "$(date): Removing duplicates, PCR artifacts, and those failing vendor quality checks..."
samtools view -b -F 1540 ${outdir}/${prefix}.dupmark.sorted.bam | samtools sort -m 50G - -o ${outdir}/${prefix}.final.bam
#rm ${outdir}/${prefix}.dupmark.sorted.bam
echo "$(date): Indexing final bam..."
bamtools index -in ${outdir}/${prefix}.final.bam
echo "$(date): DONE bam qc"
