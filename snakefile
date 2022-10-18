import glob

############################
# Variables
############################



############################
# Setup Environement
############################
rule all:
  input:
    # Ensure peppro environement is available
    "env/peppro.singularity",

    # Merged fastq files
    expand("raw_data/D-NT2_20220610/{sample}_20220610000_R{direction}.fastq.gz", sample = [4,6], direction = [1,2]),

    # Processsed fastq files
    expand("results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R{direction}.trimmed_dedupped.fastq.gz.paired.fq", sample = [4,6], direction = [1,2]),

    # Aligment Files
    expand("results/aligned_reads/D-NT2_20220610/{sample}_20220610000.sam", sample = [4,6])

rule get_peppro_environement:
  output:
    "env/peppro.singularity"
  shell:
    "wget -O {output} http://big.databio.org/simages/peppro"

# rule download_raw_data:
#   output:
#     "raw_data/exp2_72h_UL87HF_Pro-seq"

##########################
# Process fastq files
##########################

rule concatenate_fastq:
  input:
    lambda wildcards: glob.glob('raw_data/D-NT2_20220610/{sample}_lane*_20220610000_S*_L00*_R{direction}_001.fastq.gz'.format(sample=wildcards.sample, direction=wildcards.direction))
  output:
    "raw_data/D-NT2_20220610/{sample}_20220610000_R{direction}.fastq.gz"
  shell:
    "cat {input} > {output}"

rule run_trim:
    input:
        "raw_data/D-NT2_20220610/{sample}_20220610000_R1.fastq.gz",
        "raw_data/D-NT2_20220610/{sample}_20220610000_R2.fastq.gz"
    output:
        "results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R1.trimmed.fastq.gz",
        "results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R2.trimmed.fastq.gz"
    shell:
	    """
        trim_galore --paired --small_rna --dont_gzip --output_dir results/processed_fastq/D-NT2_20220610/ {input}
        gzip results/processed_fastq/D-NT2_20220610/{wildcards.sample}_20220610000_R1_val_1.fq > results/processed_fastq/D-NT2_20220610/{wildcards.sample}_20220610000_R1.trimmed.fastq.gz
        gzip results/processed_fastq/D-NT2_20220610/{wildcards.sample}_20220610000_R2_val_2.fq > results/processed_fastq/D-NT2_20220610/{wildcards.sample}_20220610000_R2.trimmed.fastq.gz
        """

rule run_dedup:
    input:
        "results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R{direction}.trimmed.fastq.gz"
    output:
        "results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R{direction}.trimmed_dedupped.fastq.gz"
    shell:
	    "seqkit rmdup -s -o {output} {input}"

rule run_pair:
    input:
        "results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R1.trimmed_dedupped.fastq.gz",
        "results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R2.trimmed_dedupped.fastq.gz"
    output:
        "results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R1.trimmed_dedupped.fastq.gz.paired.fq",
        "results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R2.trimmed_dedupped.fastq.gz.paired.fq"
    shell:
        "fastq_pair {input}"

##########################
# Map Fastq Files
##########################

rule align_reads:
    input:
        f1 = "results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R1.trimmed_dedupped.fastq.gz.paired.fq",
        f2 = "results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R2.trimmed_dedupped.fastq.gz.paired.fq"
    output:
        "results/aligned_reads/D-NT2_20220610/{sample}_20220610000.sam"
    params:
        BOWTIE_INDEX = "~/GENOMES/hg38_moth_towne/townecombined",
        UMI_SIZE = 8
    shell:
        "bowtie -x {params.BOWTIE_INDEX} --trim5 {params.UMI_SIZE} --trim3 {params.UMI_SIZE} --fr --best --sam --fullref -1 {input.f1} -2 {input.f2} {output}"
