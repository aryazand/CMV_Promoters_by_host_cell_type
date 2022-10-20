# --- Importing Configuration Files --- #

configfile: "snakemake_config.yml"

# --- Importing Python packages --- #

import glob

# --- Set Input and Output Filenames --- #

dnt2_raw_runs = ["4_lane1_20220610000_S4_L001_R1_001.fastq.gz",
                  "4_lane1_20220610000_S4_L001_R2_001.fastq.gz",
                  "4_lane2_20220610000_S11_L002_R1_001.fastq.gz",
                  "4_lane2_20220610000_S11_L002_R2_001.fastq.gz",
                  "6_lane1_20220610000_S6_L001_R1_001.fastq.gz",
                  "6_lane1_20220610000_S6_L001_R2_001.fastq.gz",
                  "6_lane2_20220610000_S13_L002_R1_001.fastq.gz",
                  "6_lane2_20220610000_S13_L002_R2_001.fastq.gz"]

dnt2_raw_files = expand("raw_data/D-NT2_20220610/{run}", run = dnt2_raw_runs)
dnt2_processed_fastq = expand("results/processed_fastq/D-NT2_20220610/{sample}_20220610000_R{direction}.trimmed_dedupped.fastq.paired.fq", sample = [4,6], direction = [1,2])
dnt2_aligned = expand("results/aligned_reads/D-NT2_20220610/{sample}_20220610000.sam", sample = [4,6])

hff_raw_files = expand("raw_data/HFF_72hr/{SRR_ID}_R{direction}.fastq.gz", SRR_ID = ["SRR13848024", "SRR13848026"], direction = [1,2])
hff_processed_fastq = expand("results/processed_fastq/HFF_72hr/{SRR_ID}_R{direction}.trimmed_dedupped.fastq.paired.fq", SRR_ID = ["SRR13848024", "SRR13848026"], direction = [1,2])
hff_aligned = expand("results/aligned_reads/HFF_72hr/{SRR_ID}.sam", SRR_ID = ["SRR13848024", "SRR13848026"])

############################
# Setup Environement
############################
rule all:
  input:
    dnt2_raw_files,
    hff_raw_files,

    # Merged fastq files
    expand("raw_data/D-NT2_20220610/{sample}_20220610000_R{direction}.fastq.gz", sample = [4,6], direction = [1,2]),

    # Processsed fastq files
    dnt2_processed_fastq,
    hff_processed_fastq,

    # Aligment Files
    dnt2_aligned,
    hff_aligned

rule download_dnt2_data:
    output:
        "raw_data/D-NT2_20220610/{run}_001.fastq.gz"
    shell:
        "wget {config[dnt2_url]}/{output}"

rule download_hff_data:
    output:
        "raw_data/HFF_72hr/{SRR_ID}_1.fastq",
        "raw_data/HFF_72hr/{SRR_ID}_2.fastq"
    log:
    	out = "sandbox/download_hff_data.{SRR_ID}.out",
      	err = "sandbox/download_hff_data.{SRR_ID}.err"
    shell:
        "fasterq-dump -e 5 -t sandbox -O raw_data/HFF_72hr {wildcards.SRR_ID}"

rule rename_hff_data:
    input:
        "raw_data/HFF_72hr/{SRR_ID}_{direction}.fastq"
    output:
        "raw_data/HFF_72hr/{SRR_ID}_R{direction}.fastq.gz"
    shell:
        "gzip {input} > {output}"

##########################
# Process fastq files
##########################

rule concatenate_fastq:
  input:
    lambda wildcards: glob.glob('raw_data/D-NT2_20220610/{sample}_lane*_20220610000_S*_L00*_R{direction}_001.fastq.gz'.format(sample=wildcards.sample, direction=wildcards.direction))
  output:
    "raw_data/D-NT2_20220610/{sample}_20220610000_R{direction}.fastq.gz"
  log:
      out = "sandbox/concatenate_fastq.{sample}_{direction}.out",
      err = "sandbox/concatenate_fastq.{sample}_{direction}.err"
  shell:
    "cat {input} > {output} 2> {log.err} 1> {log.out}"

rule run_trim:
    input:
        "raw_data/{experiment}/{sample}_R1.fastq.gz",
        "raw_data/{experiment}/{sample}_R2.fastq.gz"
    output:
        "results/processed_fastq/{experiment}/{sample}_R1_val_1.fq",
        "results/processed_fastq/{experiment}/{sample}_R2_val_2.fq"
    log:
        out = "sandbox/run_trim.{experiment}_{sample}.out",
        err = "sandbox/run_trim.{experiment}_{sample}.err"
    shell:
	    "trim_galore --paired --small_rna --dont_gzip -j 5 --output_dir results/processed_fastq/{wildcards.experiment}/ {input} 2> {log.err} 1> {log.out}"

rule rename_trimgalore_output:
    input:
        "results/processed_fastq/{experiment}/{sample}_R{direction}_val_{direction}.fq"
    output:
        "results/processed_fastq/{experiment}/{sample}_R{direction}.trimmed.fastq"
    log:
        out = "sandbox/rename_trimgalore_output.{experiment}_{sample}_R{direction}.out",
        err = "sandbox/rename_trimgalore_output.{experiment}_{sample}_R{direction}.err"
    shell:
        "mv {input} {output}  2> {log.err} 1> {log.out}"

rule run_dedup:
    input:
        "results/processed_fastq/{experiment}/{sample}.trimmed.fastq"
    output:
        "results/processed_fastq/{experiment}/{sample}.trimmed_dedupped.fastq"
    log:
        out = "sandbox/run_dedup.{experiment}_{sample}.out",
        err = "sandbox/run_dedup.{experiment}_{sample}.err"
    shell:
	    "seqkit rmdup -s -o {output} {input} 2> {log.err} 1> {log.out}"

rule run_pair:
    input:
        "results/processed_fastq/{experiment}/{sample}_R1.trimmed_dedupped.fastq",
        "results/processed_fastq/{experiment}/{sample}_R2.trimmed_dedupped.fastq"
    output:
        "results/processed_fastq/{experiment}/{sample}_R1.trimmed_dedupped.fastq.paired.fq",
        "results/processed_fastq/{experiment}/{sample}_R2.trimmed_dedupped.fastq.paired.fq"
    log:
        out = "sandbox/run_pair.{experiment}_{sample}.out",
        err = "sandbox/run_pair.{experiment}_{sample}.err"
    shell:
        "fastq_pair {input} 2> {log.err} 1> {log.out}"

##########################
# Map Fastq Files
##########################

rule align_reads:
    input:
        f1 = "results/processed_fastq/{experiment}/{sample}_R1.trimmed_dedupped.fastq.paired.fq",
        f2 = "results/processed_fastq/{experiment}/{sample}_R2.trimmed_dedupped.fastq.paired.fq"
    output:
        "results/aligned_reads/{experiment}/{sample}.sam"
    log:
        out = "sandbox/align_reads.{experiment}_{sample}.out",
        err = "sandbox/align_reads.{experiment}_{sample}.err"
    params:
        BOWTIE_INDEX = config["bowtie_index"],
        UMI_SIZE = 8
    shell:
        "bowtie -x {params.BOWTIE_INDEX} --threads 4 --trim5 {params.UMI_SIZE} --trim3 {params.UMI_SIZE} --fr --best --sam --fullref -1 {input.f1} -2 {input.f2} {output} 2> {log.err} 1> {log.out}"
