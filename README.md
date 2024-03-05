# RNA-seq Demonstration Pipeline
### Bash scripts describing a simple RNA-seq pipeline, from adapter trimming to visualization via a genome browser.
#### Workflow steps:
1. Adapter trimming ([Trim Galore!](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md))
2. Alignment ([STAR](https://github.com/alexdobin/STAR))
3. Quality Control ([samtools](https://github.com/samtools/samtools), [Picard](https://broadinstitute.github.io/picard/), [bamtools](https://github.com/pezmaster31/bamtools))
4. Quantification ([Subread](https://subread.sourceforge.net/), awk)
5. Visualization
