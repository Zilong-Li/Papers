def get_mother_kid_bam_in_fam(wc):
    mbam = NIPTFAM[wc.fam].get("Mbam")
    kbam = NIPTFAM[wc.fam].get("Kbam")
    return {"mbam": mbam, "kbam": kbam}


rule downsample_bam_in_fam:
    input:
        unpack(get_mother_kid_bam_in_fam),
    output:
        bam=os.path.join(OUTDIR, "downsample", "{fam}_{depth}x_{ff}f_{chrom}.bam"),
        bai=os.path.join(OUTDIR, "downsample", "{fam}_{depth}x_{ff}f_{chrom}.bam.bai"),
    params:
        N="downsample_bam",
        mid=lambda wildcards: NIPTFAM[wildcards.fam]["Mid"],
        mdepth=lambda wildcards: NIPTFAM[wildcards.fam]["Mdepth"],
        kid=lambda wildcards: NIPTFAM[wildcards.fam]["Kid"],
        kdepth=lambda wildcards: NIPTFAM[wildcards.fam]["Kdepth"],
    log:
        os.path.join(OUTDIR, "downsample", "{fam}_{depth}x_{ff}f_{chrom}.bam.llog"),
    threads: 1
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        MOUT={output.bam}.{params.mid}.bam
        KOUT={output.bam}.{params.kid}.bam
        FRAC=$(echo "scale=4 ; {wildcards.depth} * (1 - {wildcards.ff}) / {params.mdepth}" | bc -l)
        samtools view -s $FRAC -o $MOUT {input.mbam} {wildcards.chrom} && samtools index $MOUT
        FRAC=$(echo "scale=4 ; {wildcards.depth} * {wildcards.ff} / {params.kdepth}" | bc -l)
        samtools view -s $FRAC -o $KOUT {input.kbam} {wildcards.chrom} && samtools index $KOUT
        samtools merge -f -c -p --no-PG -o {output.bam} $MOUT $KOUT
        samtools view -H {output.bam} | sed '/^@RG/s/SM:.[^\tCN]*/SM:{wildcards.fam}/' > {output.bam}.h
        samtools reheader {output.bam}.h {output.bam} > {output.bam}.tmp.bam
        mv {output.bam}.tmp.bam {output.bam} && samtools index {output.bam}
        """


rule bamlist:
    input:
        expand(
            rules.downsample_bam_in_fam.output.bam, fam=NIPTFAM.keys(), allow_missing=True
        ),
    output:
        bam=os.path.join(OUTDIR, "downsample", "{depth}x_{ff}f_{chrom}.bamlist"),
        ff=os.path.join(OUTDIR, "downsample", "{depth}x_{ff}f_{chrom}.fflist"),
    params:
        N="bamlist",
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        echo {input} | tr ' ' '\\n' > {output.bam}
        for i in $(cat {output.bam});do echo {wildcards.ff};done > {output.ff}
        """


