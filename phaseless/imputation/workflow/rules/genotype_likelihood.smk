
rule bcftools_prepare_glvcf:
    input:
        bams=rules.bamlist.output,
        sites=rules.get_pos_from_refpanel.output.sites,
        tsv=rules.get_pos_from_refpanel.output.tsv,
    output:
        vcf=os.path.join(OUTDIR, "glvcf", "down{depth}x.{chrom}.vcf.gz"),
        csi=os.path.join(OUTDIR, "glvcf", "down{depth}x.{chrom}.vcf.gz.csi"),
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
            -I -E -A -a 'FORMAT/DP' -r {wildcards.chrom} -T {input.sites} \
            -b {input.bams} -Ou | bcftools call -Am -C alleles \
            -T {input.tsv} -Oz -o {output.vcf} && bcftools index -f {output.vcf} \
        ) &> {log}
        """

rule vcf2beagle:
    input:
        rules.bcftools_prepare_glvcf.output.vcf,
    output:
        os.path.join(OUTDIR, "glvcf", "{chrom}.down{depth}x.beagle.gz"),
    shell:
        """
        vcf2beagle -i {input} -o {output} -r {wildcards.chrom} -t PL
        """
