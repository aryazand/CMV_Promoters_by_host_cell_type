############################
# Variables 
############################



############################
# Setup Environement
############################
rule all: 
    "env/peppro.singularity"

rule get_peppro_environement:
  output:
    "env/peppro.singularity"
  shell:
    "wget -O {output} http://big.databio.org/simages/peppro"
