###############
# RNA-seq demo pipeline
# Step 1: trim_galore
###############

if [[ $# -ne 4 ]]; then
	echo "Usage: ./trim_galore.sh <fastq1> <fastq2> <outdir> <cpu>"
	exit 1
fi

fq1=${1}
fq2=${2}
outdir=${3}
cpu=${4}

prefix=$(basename ${fq1} .fastq.gz)

trim_galore --fastqc \
	-j ${cpu} \
	-q 15 \
	--stringency 5 \
	-e 0.1 \
	--length 20 \
	--illumina \
	--paired ${fq1} ${fq2} \
	-o ${outdir}/${prefix}
