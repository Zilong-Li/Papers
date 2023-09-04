
rule stitch_by_chrom:
    input:
        bamlist=rules.bamlist.output,
        posfile=rules.get_stitch_pos_file.output,
    output:
        vcf=os.path.join(OUTDIR, "stitch", "down{depth}x.{chrom}.vcf.gz"),
    log:
        os.path.join(OUTDIR, "stitch", "down{depth}x.{chrom}.vcf.gz.llog"),
    params:
        time=config["time"],
        bin=config["stitch"]["bin"],
        buffer=config["stitch"]["buffer"],
        method=config["stitch"]["method"],
        K=config["stitch"]["K"],
        nGen=config["stitch"]["nGen"],
        outdir=lambda wildcards, output: os.path.dirname(os.path.realpath(output[0])),
    threads: config["threads"]
    shell:
        """
        {params.time} -v {params.bin} \
        --bamlist={input.bamlist} \
        --posfile={input.posfile} \
        --chr={wildcards.chrom} \
        --method={params.method} \
        --K={params.K} \
        --nGen={params.nGen} \
        --nCores={threads} \
        --outputdir={params.outdir} &> {log}
        """
