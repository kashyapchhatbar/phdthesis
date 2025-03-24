import yaml

with open(config["chip_input"]) as file:
    chip_input = yaml.load(file, Loader=yaml.FullLoader)

rule fold_change_control:
    input: expand("results/bigwig_average/{sample}.bed", sample=chip_input.keys())

rule fold_change_control_cut:
    input: expand("results/bigwig_average_cut/{sample}.bed", sample=chip_input.keys())

rule fold_change_control_genes:
    input: expand("results/bigwig_average_genes/{sample}.bed", sample=chip_input.keys())    

rule bigwig_average_over_control_bed:
    input: 
        bw="results/log2_chip_control/{sample}.bw",
        quartiles="resources/old_assembly/mouse_quartiles.bed"
    output: "results/bigwig_average/{sample}.bed"
    conda: "../envs/ucsc_tools.yml"
    threads: 1
    shell: "bigWigAverageOverBed {input.bw} {input.quartiles} {output}"

rule bigwig_average_over_control_bed_cut:
    input: 
        bw="results/log2_chip_control/{sample}.bw",
        quartiles="resources/old_assembly/mouse_pgb.bed"
    output: "results/bigwig_average_cut/{sample}.bed"
    conda: "../envs/ucsc_tools.yml"
    threads: 1
    shell: "bigWigAverageOverBed {input.bw} {input.quartiles} {output}"

rule bigwig_average_over_control_bed_genes:
    input: 
        bw="results/log2_chip_control/{sample}.bw",
        quartiles="resources/old_assembly/mouse_genes.bed"
    output: "results/bigwig_average_genes/{sample}.bed"
    conda: "../envs/ucsc_tools.yml"
    threads: 1
    shell: "bigWigAverageOverBed {input.bw} {input.quartiles} {output}"