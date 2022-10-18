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
    expand("raw_data/D-NT2_20220610/{sample}_20220610000_R{direction}.fastq.gz", sample = [4,6], direction = [1,2])  

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
    glob.glob("raw_data/D-NT2_20220610/{sample}_lane*_20220610000_S*_L00*_R{direction}_001.fastq.gz")
  output:
    "raw_data/D-NT2_20220610/{sample}_20220610000_R{direction}.fastq.gz"
  shell:
    "cat {input} > {output}"

############################
# Run Peppro
############################
