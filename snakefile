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
dnt2_aligned = expand("results/aligned_reads/D-NT2_20220610/{sample}_20220610000_cmv.bed", sample = [4,6])
dnt2_5prime_coverage = expand("results/aligned_reads/D-NT2_20220610/{sample}_20220610000_5prime_coverage_{strand}_strand.tss", sample = [4,6], strand = ["pos", "minus"])

hff_raw_files = expand("raw_data/HFF_72hr/{SRR_ID}_R{direction}.fastq.gz", SRR_ID = ["SRR13848024", "SRR13848026"], direction = [1,2])
hff_processed_fastq = expand("results/processed_fastq/HFF_72hr/{SRR_ID}_R{direction}.trimmed_dedupped.fastq.paired.fq", SRR_ID = ["SRR13848024", "SRR13848026"], direction = [1,2])
hff_aligned = expand("results/aligned_reads/HFF_72hr/{SRR_ID}_cmv.bed", SRR_ID = ["SRR13848024", "SRR13848026"])
hff_5prime_coverage = expand("results/aligned_reads/HFF_72hr/{SRR_ID}_5prime_coverage_{strand}_strand.tss", SRR_ID = ["SRR13848024", "SRR13848026"], strand = ["pos", "minus"])

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
    dnt2_processed_fastq,
    hff_processed_fastq,

    # Aligment Files
    dnt2_aligned,
    hff_aligned,

    #5 prime coverage files
    dnt2_5prime_coverage,
    hff_5prime_coverage

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
        "bowtie -x {params.BOWTIE_INDEX} --threads 4 --trim5 {params.UMI_SIZE} --trim3 {params.UMI_SIZE} --chunkmbs 500 --fr --best --sam --fullref -1 {input.f1} -2 {input.f2} {output} 2> {log.err} 1> {log.out}"

rule sam_to_bam:
  input:
      "results/aligned_reads/{experiment}/{sample}.sam"
  output:
      "results/aligned_reads/{experiment}/{sample}_allgenomes.bam"
  log:
      out = "sandbox/sam_to_bam.{experiment}_{sample}.out",
      err = "sandbox/sam_to_bam.{experiment}_{sample}.err"
  shell:
      "samtools view --threads 5 -u {input} | samtools sort --threads 5 -o {output}"

rule index_bam:
  input:
      "results/aligned_reads/{experiment}/{sample}_allgenomes.bam"
  output:
      "results/aligned_reads/{experiment}/{sample}_allgenomes.bam.bai"
  log:
      out = "sandbox/index_bam.{experiment}_{sample}.out",
      err = "sandbox/index_bam.{experiment}_{sample}.err"
  shell:
      "samtools index {input} 2> {log.err} 1> {log.out}"

rule extract_only_cmv_alignment:
  input:
      bam = "results/aligned_reads/{experiment}/{sample}_allgenomes.bam",
      bai = "results/aligned_reads/{experiment}/{sample}_allgenomes.bam.bai"
  output:
      "results/aligned_reads/{experiment}/{sample}_cmv.bam"
  log:
      err = "sandbox/extract_only_cmv_alignment.{experiment}_{sample}.err"
  shell:
      "samtools view --threads 5 -b {input.bam} FJ616285.1 > {output} 2> {log.err}"

rule bam_to_bed:
  input:
      "results/aligned_reads/{experiment}/{sample}_cmv.bam"
  output:
      "results/aligned_reads/{experiment}/{sample}_cmv.bed"
  log:
      err = "sandbox/extract_only_cmv_alignment.{experiment}_{sample}.err"
  shell:
      "bedtools bamtobed -i {input} > {output} 2> {log.err}"

rule bam_to_5prime_coverage:
  input:
      "results/aligned_reads/{experiment}/{sample}_cmv.bed"
  output:
      plus_strand = "results/aligned_reads/{experiment}/{sample}_5prime_coverage_pos_strand.tss",
      minus_strand = "results/aligned_reads/{experiment}/{sample}_5prime_coverage_minus_strand.tss"
  log:
      err = "sandbox/bam_to_5prime_coverage.{experiment}_{sample}.err"
  params:
      ref_cmv_genome = "./raw_data/genomic_data/FJ616285.1.chrom.sizes"
  shell:
      """
      bedtools genomecov -5 -d -strand + -i {input} -g {params.ref_cmv_genome} > {output.plus_strand} 2> {log.err}
      bedtools genomecov -5 -d -strand - -i {input} -g {params.ref_cmv_genome} > {output.minus_strand} 2> {log.err}
      """
