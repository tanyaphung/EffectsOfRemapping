import os

configfile: "/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/process_rna_config.json"

rule all:
    input:
        expand("featureCounts_geneid/{sample}_XX_HISAT2_transcript_featurecounts.tsv", sample=config["all_rna_samples"])
    input:
        expand("featureCounts/{sample}_XX_HISAT2_transcript_featurecounts.tsv", sample=config["all_rna_samples"])

# featureCounts
rule xx_featurecounts:
    input:
        BAM = "/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/processed_bams/rna/{sample}.GRCh38.p12.genome.XXonly.sorted.merged.bam",
        GTF = "/scratch/tphung3/SmallProjects/CompareTPM/gencode.v26.annotation.gtf"
    output:
        COUNTS = "featureCounts/{sample}_XX_HISAT2_transcript_featurecounts.tsv"
    params:
        THREADS = 5
    shell:
        """
        featureCounts -T {params.THREADS} --primary -p -s 2 -t exon -g transcript_id -O -a {input.GTF} -o {output.COUNTS} {input.BAM}
        """

# featureCounts
rule xx_featurecounts_geneid:
    input:
        BAM = "/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/processed_bams/rna/{sample}.GRCh38.p12.genome.XXonly.sorted.merged.bam",
        GTF = "/scratch/tphung3/SmallProjects/CompareTPM/gencode.v26.annotation.gtf"
    output:
        COUNTS = "featureCounts_geneid/{sample}_XX_HISAT2_transcript_featurecounts.tsv"
    params:
        THREADS = 5
    shell:
        """
        featureCounts -T {params.THREADS} --primary -p -s 2 -t exon -g gene_id -O -a {input.GTF} -o {output.COUNTS} {input.BAM}
        """
