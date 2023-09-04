
rule phaseless_impute_by_chrom:
    input:
        rules.vcf2beagle.output,
    output:
        vcf=os.path.join(OUTDIR, "phaseless", "down{depth}x.{chrom}.vcf.gz"),
    log:
        os.path.join(OUTDIR, "phaseless", "down{depth}x.{chrom}.vcf.gz.llog"),
    params:
        time=config["time"],
        bin=config["phaseless"]["bin"],
        C=config["phaseless"]["C"],
        iters=config["phaseless"]["iters"],
        out=lambda wildcards, output: output[0][:-7],
    threads: config["threads"]
    shell:
        """
        {params.time} -v {params.bin} impute -g {input} -c {params.C} -n {threads} -S -i {params.iters} -o {params.out} &> {log}
        """
