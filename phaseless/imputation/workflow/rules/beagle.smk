

rule beagle41_by_chrom:
    input:
        rules.bcftools_prepare_glvcf.output.vcf,
    output:
        vcf=os.path.join(OUTDIR, "beagle4.1", "down{depth}x.{chrom}.vcf.gz"),
    log:
        os.path.join(OUTDIR, "beagle4.1", "down{depth}x.{chrom}.vcf.gz.llog"),
    params:
        time=config["time"],
        bin=config["beagle"]["version41"],
        out=lambda wildcards, output: output[0][:-7],
    threads: config["threads"]
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        {params.time} -v java -Xss5m -Xmx60g -jar {params.bin} gl={input} out={params.out} chrom={wildcards.chrom} nthreads={threads} &> {log}
        """
