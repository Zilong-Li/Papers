

rule subset_sample_list:
    """put sample to be used in a file in case of arguments too long"""
    output:
        samples=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.kept.samples"),
    params:
        N="subset_sample_list",
        samples=get_samples_in_nipt_fam(),
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
    log:
        os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.kept.samples.llog"),
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/subset_samples.R"


rule subset_refpanel_by_region:
    input:
        rules.subset_sample_list.output.samples,
    output:
        vcf=os.path.join(
            OUTDIR_PANEL, "refsize{size}", "vcfs", "{chrom}.{start}.{end}.vcf.gz"
        ),
        csi=os.path.join(
            OUTDIR_PANEL, "refsize{size}", "vcfs", "{chrom}.{start}.{end}.vcf.gz.csi"
        ),
        hap=os.path.join(
            OUTDIR_PANEL, "refsize{size}", "vcfs", "{chrom}.{start}.{end}.vcf.hap.gz"
        ),
        leg=os.path.join(
            OUTDIR_PANEL,
            "refsize{size}",
            "vcfs",
            "{chrom}.{start}.{end}.vcf.legend.gz",
        ),
        sam=os.path.join(
            OUTDIR_PANEL, "refsize{size}", "vcfs", "{chrom}.{start}.{end}.vcf.samples"
        ),
        sites=os.path.join(
            OUTDIR_PANEL,
            "refsize{size}",
            "vcfs",
            "{chrom}.{start}.{end}.sites.vcf.gz",
        ),
        tsv=os.path.join(
            OUTDIR_PANEL, "refsize{size}", "vcfs", "{chrom}.{start}.{end}.tsv.vcf.gz"
        ),
    params:
        N="subset_refpanel_by_region",
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        start=lambda wildcards: max(
            1, int(wildcards.start) - int(config["quilt1"]["buffer"])
        ),
        end=lambda wildcards: int(wildcards.end) + int(config["quilt1"]["buffer"]),
    log:
        os.path.join(
            OUTDIR_PANEL,
            "refsize{size}",
            "vcfs",
            "{chrom}.{start}.{end}.vcf.gz.llog",
        ),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        bcftools view -v snps -m2 -M2 --samples-file {input} --threads 4 {params.vcf} {wildcards.chrom}:{params.start}-{params.end}| bcftools norm - -d snps -Ob -o {output.vcf} --threads 4 && bcftools index -f {output.vcf}
        touch -m {output.vcf}.csi
        bcftools convert --haplegendsample {params.prefix} {output.vcf}
        bcftools view -G {output.vcf} -Oz -o {output.sites} --threads 4 && tabix -f {output.sites}
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sites} | bgzip -c > {output.tsv} && tabix -s1 -b2 -e2 {output.tsv}
        """


def get_sites_subset_refpanel_by_region(wildcards):
    starts, ends = get_refpanel_chunks(wildcards.chrom)
    return expand(
        rules.subset_refpanel_by_region.output.sites,
        zip,
        start=starts,
        end=ends,
        allow_missing=True,
    )


rule concat_refpanel_sites_by_region:
    input:
        get_sites_subset_refpanel_by_region,
    output:
        sites=os.path.join(
            OUTDIR_PANEL,
            "refsize{size}",
            "vcfs",
            "{chrom}.sites.vcf.gz",
        ),
        tsv=os.path.join(OUTDIR_PANEL, "refsize{size}", "vcfs", "{chrom}.tsv.vcf.gz"),
    log:
        os.path.join(OUTDIR_PANEL, "refsize{size}", "vcfs", "{chrom}.sites.vcf.gz.llog"),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        ( \
            echo {input} | tr ' ' '\n' > {output.sites}.list && \
            bcftools concat -f {output.sites}.list -Da --threads 4 -Oz -o {output.sites} && tabix -f {output.sites} && \
            bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sites} | bgzip -c > {output.tsv} && tabix -s1 -b2 -e2 {output.tsv}
        ) & > {log}
        """

rule collect_truth_gts:
    """would be better to use sites in subrefs"""
    input:
        sites=lambda wildcards: expand(
            rules.concat_refpanel_sites_by_region.output.sites,
            size=config["refsize"],
            allow_missing=True,
        ),
    output:
        gt=os.path.join(OUTDIR_TRUTH, "truth.gts.{chrom}.txt"),
        af=os.path.join(OUTDIR_TRUTH, "af.input.panel.{chrom}.txt"),
        tmp=temp(os.path.join(OUTDIR_TRUTH, "af.input.panel.{chrom}.txt.tmp")),
        tmp2=temp(os.path.join(OUTDIR_TRUTH, "truth.gts.{chrom}.txt.tmp")),
    log:
        os.path.join(OUTDIR_PANEL, "truth.gts.{chrom}.log"),
    params:
        N="collect_truth_gts",
        samples=",".join(get_samples_in_nipt_fam()),
        truth=lambda wildcards: REFPANEL[wildcards.chrom]["truth"],
        ref=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        af=if_use_af_in_refpanel,
        ql0="%CHROM:%POS:%REF:%ALT\\n",
        ql1="%CHROM:%POS:%REF:%ALT\\t%AF\\n",
        ql2="%CHROM:%POS:%REF:%ALT[\\t%GT]\\n",
        awk="NR==FNR{a[$1]=1;} NR!=FNR{if(a[$1]){print $1,$2;}}",
        awk2="NR==FNR{a[$1]=1;} NR!=FNR{if(a[$1]){print $0;}}",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        (
        if [ -s {params.af} ];then perl -lane 'print join(":",@F[0..3])."\\t$F[4]"' {params.af} > {output.tmp}; \
        else \
            bcftools +fill-tags {input.sites[0]} -- -t AF |  bcftools query -f '{params.ql1}' > {output.tmp}; \
        fi
        awk '{params.awk}' <(bcftools query -f '{params.ql0}' {input.sites[0]}) {output.tmp} >{output.af}
        bcftools view -s {params.samples} {params.truth} | bcftools query -f '{params.ql2}' > {output.tmp2}
        HEADER="ID,"{params.samples}
        echo $HEADER | tr ',' '\t' > {output.gt}
        awk '{params.awk2}' <(bcftools query -f '{params.ql0}' {input.sites[0]}) {output.tmp2} >> {output.gt}
        ) &> {log}
        """
