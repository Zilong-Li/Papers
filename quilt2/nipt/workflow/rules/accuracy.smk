
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


rule plot_quilt_regular_by_chunk:
    input:
        truth=rules.collect_truth_gts.output.gt,
        af=rules.collect_truth_gts.output.af,
        vcf=rules.quilt_run_regular.output,
    output:
        rds=os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.{ff}f.regular.{chrom}.{start}.{end}.vcf.gz.rds",
        ),
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy_nipt.R"


rule collect_quilt_regular_accuracy:
    input:
        get_quilt_regular_accuracy_by_chunk
    output:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.{ff}f.regular.{chrom}.vcf.gz.rds.list",
        ),
    conda:
        "../envs/quilt.yaml"
    shell:
        "echo {input} | tr ' ' '\n' > {output}"


rule plot_quilt_regular:
    input:
        truth=rules.collect_truth_gts.output.gt,
        af=rules.collect_truth_gts.output.af,
        vcf=rules.quilt_ligate_regular.output.vcf,
    output:
        rds=os.path.join(
            OUTDIR_SUMMARY, "quilt.nipt.accuracy.regular.refsize{size}.{chrom}.down{depth}x.{ff}f.rds"
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY, "quilt.nipt.accuracy.regular.refsize{size}.{chrom}.down{depth}x.{ff}f.rds.llog"
        ),
    params:
        N="plot_quilt_regular",
    resources:
        slots=1,
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy_nipt.R"



rule collapse_accuracy_by_refsize:
    input:
        expand(rules.plot_quilt_regular.output.rds,
               depth=config["coverage"], ff=config["fetalfrac"]
               ,allow_missing=True),
    output:
        rds=os.path.join(
            OUTDIR_SUMMARY, "quilt.nipt.accuracy.regular.refsize{size}.{chrom}.rds"
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY, "quilt.nipt.accuracy.regular.refsize{size}.{chrom}.rds.llog"
        ),
    params:
        N="collapse_accuracy_by_refsize",
    resources:
        slots=1,
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy_collapse.R"



