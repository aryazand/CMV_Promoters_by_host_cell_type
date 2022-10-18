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

dnt2_raw_files = expand("raw_data/D-NT2_20220610/{run}", run = dnt2_raw_files)

hff_raw_files = expand("raw_data/HFF_72hr/{SRR_ID}_{direction}.fastq.gz", SRR_ID = [SRR13848024, SRR13848026], direction = [1,2])

############################
# Setup Environement
############################
rule all:
  input:
    "env/peppro.singularity",
    dnt2_raw_files,
    hff_raw_files
    expand("raw_data/D-NT2_20220610/{sample}_20220610000_R{direction}.fastq.gz", sample = [4,6], direction = [1,2])

rule get_peppro_environement:
  output:
    "env/peppro.singularity"
  shell:
    "wget -O {output} http://big.databio.org/simages/peppro"

rule download_dnt2_data:
    output:
        dnt2_raw_files
    shell:
        "wget -r --no-parent -nH -A [46]_*.fastq.gz {config[dnt2_url]}"

rule download_hff_data:
    output:
        "raw_data/HFF_72hr/{SRR_ID}_1.fastq.gz",
        "raw_data/HFF_72hr/{SRR_ID}_2.fastq.gz"
    shell:
        """
        fasterq-qdump -e 5 -t sandbox -O raw_data/HFF_72hr {wildcards.SRR_ID}
        gzip raw_data/HFF_72hr/{wildcards.SRR_ID}_1.fastq
        gzip raw_data/HFF_72hr/{wildcards.SRR_ID}_2.fastq
        """

rule concatenate_fastq:
  input:
    lambda wildcards: glob.glob('raw_data/D-NT2_20220610/{sample}_lane*_20220610000_S*_L00*_R{direction}_001.fastq.gz'.format(sample=wildcards.sample, direction=wildcards.direction))
  output:
    "raw_data/D-NT2_20220610/{sample}_20220610000_R{direction}.fastq.gz"
  shell:
    "cat {input} > {output}"
