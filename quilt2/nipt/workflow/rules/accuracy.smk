

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
        method=config["quilt1"]["method"],
        truth=lambda wildcards: REFPANEL[wildcards.chrom]["truth"],
    resources:
        slots=1,
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy_nipt.R"



rule collapse_regular_accuracy_by_refsize:
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



