
rule phaseless_v0_impute_by_chrom:
    input:
        rules.vcf2beagle.output,
    output:
        vcf=os.path.join(OUTDIR, "phaseless_v0", "down{depth}x.{chrom}.vcf.gz"),
    log:
        os.path.join(OUTDIR, "phaseless_v0", "down{depth}x.{chrom}.vcf.gz.llog"),
    params:
        time=config["time"],
        bin=config["phaseless"]["v0"],
        C=config["phaseless"]["C"],
        chunksize=config["phaseless"]["chunksize"],
        iters=config["phaseless"]["iters"],
        out=lambda wildcards, output: output[0][:-7],
    threads: config["threads"]
    shell:
        """
        ( 
        {params.time} -v {params.bin} impute -g {input} -c {params.C} -n {threads} -S -i {params.iters} -o {params.out} && \
        bcftools index -f {output.vcf}
        ) &> {log}
        """

rule phaseless_v1_impute_by_chrom:
    input:
        rules.vcf2beagle.output,
    output:
        vcf=os.path.join(OUTDIR, "phaseless_v1", "down{depth}x.{chrom}.vcf.gz"),
    log:
        os.path.join(OUTDIR, "phaseless_v1", "down{depth}x.{chrom}.vcf.gz.llog"),
    params:
        time=config["time"],
        bin=config["phaseless"]["v1"],
        C=config["phaseless"]["C"],
        chunksize=config["phaseless"]["chunksize"],
        iters=config["phaseless"]["iters"],
        out=lambda wildcards, output: output[0][:-7],
    threads: config["threads"]
    shell:
        """
        ( 
        {params.time} -v {params.bin} impute -g {input} -c {params.C} -n {threads} -S -i {params.iters} -o {params.out} && \
        bcftools index -f {output.vcf}
        ) &> {log}
        """
