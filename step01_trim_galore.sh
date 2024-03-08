###############
# RNA-seq demo pipeline
# Step 1: trim_galore
###############
# Parameter check
if [[ $# -ne 4 ]]; then
	echo "Usage: ./trim_galore.sh <fastq1> <fastq2> <outdir> <cpu>"
	exit 1
fi
# Define parameters
fq1=${1}
fq2=${2}
outdir=${3}
cpu=${4}
# Define a name for our sample
prefix=$(basename ${fq1} _R1.fastq.gz)
# Call the tool with appropriate parameters
trim_galore --fastqc \
	-j ${cpu} \
	-q 15 \
	--stringency 5 \
	-e 0.1 \
	--length 20 \
	--illumina \
	--paired ${fq1} ${fq2} \
	-o ${outdir}/${prefix}
