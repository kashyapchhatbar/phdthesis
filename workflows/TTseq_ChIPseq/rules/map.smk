rule star:
    input:
        fq1="fastq/{sample}_1.fq.gz",
        fq2="fastq/{sample}_2.fq.gz",
        idx="resources/genomes/concatenated",
    output:
        bam=temp("results/star/concatenated/{sample}/Aligned.out.bam"),
        unmapped1=temp("results/star/concatenated/{sample}/Unmapped.out.mate1"),
        unmapped2=temp("results/star/concatenated/{sample}/Unmapped.out.mate2")
    log:
        "logs/star/concatenated/{sample}.log",
    params:
        prefix="results/star/concatenated/{sample}/"          
    threads: 24
    conda:
        "../envs/star.yml"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.idx} \
            --readFilesIn {input.fq1} {input.fq2} --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.prefix} \
            --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic \
            --outReadsUnmapped Fastx
        """

rule extract_star_unmapped:
    input: 
        unmapped1="results/star/concatenated/{sample}/Unmapped.out.mate1",
        unmapped2="results/star/concatenated/{sample}/Unmapped.out.mate2"
    output:
        unmapped1=temp("results/unmapped/concatenated/{sample}.1.fastq.gz"),
        unmapped2=temp("results/unmapped/concatenated/{sample}.2.fastq.gz")
    threads: 8
    shell:
        """
        pigz -p {threads} -c {input.unmapped1} > {output.unmapped1}
        pigz -p {threads} -c {input.unmapped2} > {output.unmapped2}
        """

rule bowtie2:
    input:
        fq1="results/unmapped/concatenated/{sample}.1.fastq.gz",
        fq2="results/unmapped/concatenated/{sample}.2.fastq.gz",
        idx=multiext(
            "resources/genomes/concatenated",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        ref="resources/genomes/concatenated.fa", #Required for CRAM output
    output:
        cram="results/sorted/backup/{sample}_bowtie2.cram",
        fq1="results/unmapped/concatenated/{sample}_1.fq.gz",
        fq2="results/unmapped/concatenated/{sample}_2.fq.gz"
    params:
        un="results/unmapped/concatenated/{sample}_%.fq.gz",
        index="resources/genomes/concatenated"
    log:
        "logs/bowtie2/{sample}.log",    
    threads: 24
    conda:
        "../envs/bowtie2.yml"
    shell:
        """
        bowtie2 --no-unal --un-conc-gz {params.un} --sensitive-local --no-mixed \
        --no-discordant -p {threads} -x {params.index} -1 {input.fq1} \
        -2 {input.fq2} | samtools sort -@ {threads} -m 6G --reference \
        {input.ref} --write-index -o {output.cram} -
        """

rule samtools:
    input: 
        bam="results/star/concatenated/{sample}/Aligned.out.bam",
        ref="resources/genomes/concatenated.fa"
    output: 
        cram="results/sorted/backup/{sample}.cram",
        idx="results/sorted/backup/{sample}.cram.crai",
    log:
        "logs/samtools/cram_concatenated_{sample}.log",      
    threads: 12
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools sort --write-index -@ {threads} --reference {input.ref} \
            -m 6G -o {output.cram} {input.bam}
        """

rule merge_crams:
    input:
        star="results/sorted/backup/{sample}.cram",
        bowtie2="results/sorted/backup/{sample}_bowtie2.cram",
        ref="resources/genomes/concatenated.fa"
    output:
        merged="results/sorted/concatenated/{sample}.cram",
    conda:
        "../envs/samtools.yml"
    threads: 12
    shell:
        """
        samtools cat -@ 6 {input.star} {input.bowtie2} | samtools view -@ 2 \
            -m 2G -h -F 256 -T {input.ref} - | samtools sort -@ 6 -m 6G \
            --reference {input.ref} --write-index -o {output.merged} -
        """

rule separate_bams:
    input:
        merged="results/sorted/concatenated/{sample}.cram",
        ref="resources/genomes/concatenated.fa",
        regions="resources/genomes/{genome}_chrs.bed"
    output:
        "results/sorted/{genome}/{sample}.bam"
    conda:
        "../envs/samtools.yml"
    threads: 12
    shell:
        """
        samtools view --write-index -T {input.ref} -hb -@ {threads} -m 6G -P \
            --regions-file {input.regions} {input.merged} -o {output}
        """

rule deduplicate:
    input:
        bam="results/sorted/{genome}/{sample}.bam"        
    output:
        bam="results/deduplicated/{genome}/{sample}.bam",
        stats="results/deduplicated/{genome}/{sample}.stats"
    conda:
        "../envs/samtools.yml"
    threads: 12
    shell:
        """
        samtools collate -@ 12 -O {input.bam} \
            | samtools fixmate -m - - | samtools \
            sort -@ 12 -m 4G - | samtools markdup --write-index -r -@ 12 -f \
            {output.stats} - {output.bam}
        """

ruleorder: bowtie2 > samtools > merge_crams > separate_bams
