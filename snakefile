import glob

############################
# Variables
############################



############################
# Setup Environement
############################
rule all:
  input:
    "env/peppro.singularity",
    expand("raw_data/D-NT2_20220610/{sample}_20220610000_R{direction}.fastq.gz", sample = [4,6], direction = [1,2]),
    expand("results/peppro_results/D-NT2_20220610/{sample}/PEPPRO_completed.flag", sample = [4,6])

rule get_peppro_environement:
  output:
    "env/peppro.singularity"
  shell:
    "wget -O {output} http://big.databio.org/simages/peppro"

# rule download_raw_data:
#   output:
#     "raw_data/exp2_72h_UL87HF_Pro-seq"

rule concatenate_fastq:
  input:
    lambda wildcards: glob.glob('raw_data/D-NT2_20220610/{sample}_lane*_20220610000_S*_L00*_R{direction}_001.fastq.gz'.format(sample=wildcards.sample, direction=wildcards.direction))
  output:
    "raw_data/D-NT2_20220610/{sample}_20220610000_R{direction}.fastq.gz"
  shell:
    "cat {input} > {output}"

############################
# Run Peppro
############################

rule run_peppro:
  input:
    "raw_data/D-NT2_20220610/{sample}_20220610000_R1.fastq.gz",
    "raw_data/D-NT2_20220610/{sample}_20220610000_R2.fastq.gz"
  output:
    "results/peppro_results/D-NT2_20220610/{sample}/PEPPRO_completed.flag"
  shell:
    """
    singularity exec env/peppro.singularity ~/peppro/pipelines/peppro.py \
    -P 5 \
    --single-or-paired paired \
    --genome townecombined \
    --genome-index ~/GENOME/hg38_moth_towne/ \
    --chrom-sizes ~/GENOME/hg38_moth_towne/townecombined.chrom.sizes \
    --sample-name {wildcards.sample} \
    --input {input} \
    --output-parent results
    """
