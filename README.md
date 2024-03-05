# RNA-seq Demonstration Pipeline
### Bash scripts describing a simple RNA-seq pipeline, from adapter trimming to visualization on a genome browser.
#### Workflow steps:
1. Adapter trimming ([Trim Galore!](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md))
2. Alignment ([STAR](https://github.com/alexdobin/STAR))
3. Quality Control ([samtools](https://github.com/samtools/samtools), [Picard](https://broadinstitute.github.io/picard/), [bamtools](https://github.com/pezmaster31/bamtools))
4. Quantification ([Subread](https://subread.sourceforge.net/), awk)
5. Visualization ([bedtools](https://bedtools.readthedocs.io/en/latest/), [bedGraphToBigWig](https://www.encodeproject.org/software/bedgraphtobigwig/))

#### Detailed methods
**1. Adapter trimming is performed with Trim Galore! package.** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`-q 15` - Trim base calls that are below this threshold in addition to adapter sequences. \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`--stringency 5` - Minimum overlap of adapter sequence to be trimmed. \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`-e 0.1` - Max. allowed error rate (no. of errors divided by the length of the matching region) \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`--length 20` - Discard reads that become smaller than this after trimming \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`--illumina` - Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter `AGATCGGAAGAGC` \
**2. Alignment is performed with STAR aligner. This takes reads from our trimmed `.fastq` files and assigns them to regions on the genome.** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[parameters TBA] \
**3. Quality control of BAM file. It is important to document these steps, as they can affect reproducibility in the results, and there are many options to choose from.** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;a. Filter out any non-chromosomal sequences (contigs) and sort BAM \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;b. Mark duplicate reads (set FLAG) via Picard MarkDuplicates and sort \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;c. Remove the marked duplicates, as well as PCR artifacts, and reads failing vendor quality checks (`-F 1540`) and sort \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;d. Index the final BAM \
**4. Read counting (quantification) counts reads over gene regions, given a BAM file. We will use these to perform downstream statistical analyses.**
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
````
featureCounts \
        -p \ # Paired-end reads
        -t exon \ 
        -g gene_id \
        -s 1 \
        -O \
        -T <# of CPUs> \
        -a <GTF File> \
        -o <outfile> \
        <BAM File>
````
**5. We construct a `bigWig` file, which is a binary file format for viewing the data on a genome browser track.**

