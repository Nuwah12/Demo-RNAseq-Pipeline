###############
# Demo RNA-seq pipeline
# Step 2: Alignment with STAR
###############
if [[ $# -ne 5 ]]; then
        echo "Usage: ./step02_STAR_alignment.sh <trimmed fastq1> <trimmed fastq2> <STAR index> <outdir> <cpu>"
        exit 1
fi

fq1=${1}
fq2=${2}
STAR_index=${3}
outdir=${4}
cpu=${5}

prefix=$(basename ${fq1} .fastq.gz)

STAR \
	--genomeDir ${STAR_index} \
	--readFilesIn ${fq1} ${fq2} \
	--outFilterType BySJout \
	--readFilesCommand zcat \
	--runThreadN ${cpu} \
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--alignIntronMax 100000 \
	--outSAMstrandField intronMotif \
	--outFileNamePrefix ${outdir}/${prefix} \
	--outSAMunmapped Within \
	--chimSegmentMin 25 \
	--chimJunctionOverhangMin 25 --outStd SAM | samtools view -bS - | samtools sort -m 50G - -o ${outdir}/`basename $fq1 .fq.gz`.raw.bam

