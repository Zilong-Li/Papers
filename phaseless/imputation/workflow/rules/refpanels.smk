
rule get_pos_from_refpanel:
    output:
        sites=os.path.join(OUTDIR, "refpanel", "{chrom}.sites.vcf.gz"),
        tsv=os.path.join(OUTDIR, "refpanel", "{chrom}.sites.tsv.gz"),
    params:
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
    log:
        os.path.join(OUTDIR, "refpanel", "{chrom}.sites.tsv.gz.llog"),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        bcftools view -v snps -m2 -M2 -g het --threads 4 {params.vcf}| bcftools norm - -d snps | bcftools view -G -Oz -o {output.sites} --threads 4 && tabix -f {output.sites} && \
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sites} | bgzip -c > {output.tsv} && \
        tabix -s1 -b2 -e2 {output.tsv} &> {log}
        """


rule get_stitch_pos_file:
    input:
        rules.get_pos_from_refpanel.output.tsv,
    output:
        os.path.join(OUTDIR, "refpanel", "{chrom}.sites.posfile.txt"),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        zcat {input} | tr ',' '\t' > {output}
        """
