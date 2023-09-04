
rule bcftools_prepare_glvcf:
    input:
        bams=rules.bamlist.output,
        sites=rules.get_pos_from_refpanel.output.sites,
        tsv=rules.get_pos_from_refpanel.output.tsv,
    output:
        vcf=os.path.join(OUTDIR, "glvcf", "down{depth}x.{chrom}.bcf"),
        csi=os.path.join(OUTDIR, "glvcf", "down{depth}x.{chrom}.bcf.csi"),
    log:
        os.path.join(OUTDIR, "glvcf", "down{depth}x.{chrom}.bcf.llog"),
    params:
        time=config["time"],
        fasta=config["genome"]["fasta"],
        bq=config["bcftools"]["bq"],
        mq=config["bcftools"]["mq"],
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        ( \
        {params.time} -v bcftools mpileup -q {params.bq} -Q {params.mq} -f {params.fasta} \
            -I -E -A -a 'FORMAT/DP' -r {wildcards.chrom} -T {input.sites[0]} \
            -b {input.bams} -Ou | bcftools call -Aim -C alleles \
            -T {input.tsv[0]} -Ob -o {output.vcf} && bcftools index -f {output.vcf} \
        ) &> {log}
        """

rule vcf2beagle:
    input:
        rules.bcftools_prepare_glvcf.output.vcf,
    output:
        os.path.join(OUTDIR, "glvcf", "{chrom}", "down{depth}x.{chrom}.beagle.gz"),
    shell:
        """
        vcf2beagle -i {input} -o {output} -r {wildcards.chrom} -t PL
        """
