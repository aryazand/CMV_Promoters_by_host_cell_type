# --- Importing Configuration Files --- #

configfile: "snakemake_config.yml"

# --- Importing Python packages --- #

import glob

############################
# Variables
############################

dnt2_raw_runs = ["4_lane1_20220610000_S4_L001_R1_001.fastq.gz",
                  "4_lane1_20220610000_S4_L001_R2_001.fastq.gz",
                  "4_lane2_20220610000_S11_L002_R1_001.fastq.gz",
                  "4_lane2_20220610000_S11_L002_R2_001.fastq.gz",
                  "6_lane1_20220610000_S6_L001_R1_001.fastq.gz",
                  "6_lane1_20220610000_S6_L001_R2_001.fastq.gz",
                  "6_lane2_20220610000_S13_L002_R1_001.fastq.gz",
                  "6_lane2_20220610000_S13_L002_R2_001.fastq.gz"]

dnt2_raw_files = expand("raw_data/D-NT2_20220610/{run}", run = dnt2_raw_runs)

hff_raw_files = expand("raw_data/HFF_72hr/{SRR_ID}_R{direction}.fastq.gz", SRR_ID = ["SRR13848024", "SRR13848026"], direction = [1,2])

############################
# Setup Environement
############################
rule all:
  input:
    # Ensure peppro environement is available
    "env/peppro.singularity",
    dnt2_raw_files,
    hff_raw_files,

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

rule download_dnt2_data:
    output:
        "raw_data/D-NT2_20220610/{run}_001.fastq.gz"
    shell:
        "wget {config[dnt2_url]}/{output}"

rule download_hff_data:
    output:
        "raw_data/HFF_72hr/{SRR_ID}_1.fastq.gz",
        "raw_data/HFF_72hr/{SRR_ID}_2.fastq.gz"
    shell:
        "fasterq-qdump -e 5 -t sandbox -O raw_data/HFF_72hr {wildcards.SRR_ID}"

rule rename_hff_data:
    input:
        "raw_data/HFF_72hr/{SRR_ID}_{direction}.fastq.gz"
    output:
        "raw_data/HFF_72hr/{SRR_ID}_R{direction}.fastq.gz"
    shell:
        "gzip {input} > {output} && rm {input}"

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
        "raw_data/{experiment}/{sample}_R1.fastq.gz",
        "raw_data/{experiment}/{sample}_R2.fastq.gz"
    output:
        "results/processed_fastq/{experiment}/{sample}_R1_val_1.fq",
        "results/processed_fastq/{experiment}/{sample}_R2_val_2.fq"
    shell:
	    "trim_galore --paired --small_rna --dont_gzip --output_dir results/processed_fastq/D-NT2_20220610/ {input}"

rule rename_trimgalore_output:
    input:
        "results/processed_fastq/{experiment}/{sample}_R{direction}_val_{direction}.fq"
    output:
        "results/processed_fastq/{experiment}/{sample}_R{direction}.trimmed.fastq"
    shell:
        "mv {input} {output}"

rule run_dedup:
    input:
        "results/processed_fastq/{experiment}/{sample}.trimmed.fastq.gz"
    output:
        "results/processed_fastq/{experiment}/{sample}.trimmed_dedupped.fastq.gz"
    shell:
	    "seqkit rmdup -s -o {output} {input}"

rule run_pair:
    input:
        "results/processed_fastq/{experiment}/{sample}_R1.trimmed_dedupped.fastq.gz",
        "results/processed_fastq/{experiment}/{sample}_R2.trimmed_dedupped.fastq.gz"
    output:
        "results/processed_fastq/{experiment}/{sample}_R1.trimmed_dedupped.fastq.gz.paired.fq",
        "results/processed_fastq/{experiment}/{sample}_R2.trimmed_dedupped.fastq.gz.paired.fq"
    shell:
        "fastq_pair {input}"

##########################
# Map Fastq Files
##########################

rule align_reads:
    input:
        f1 = "results/processed_fastq/{experiment}/{sample}_R1.trimmed_dedupped.fastq.gz.paired.fq",
        f2 = "results/processed_fastq/{experiment}/{sample}_R2.trimmed_dedupped.fastq.gz.paired.fq"
    output:
        "results/aligned_reads/{experiment}/{sample}.sam"
    params:
        BOWTIE_INDEX = config["bowtie_index"],
        UMI_SIZE = 8
    shell:
        "bowtie -x {params.BOWTIE_INDEX} --trim5 {params.UMI_SIZE} --trim3 {params.UMI_SIZE} --fr --best --sam --fullref -1 {input.f1} -2 {input.f2} {output}"
