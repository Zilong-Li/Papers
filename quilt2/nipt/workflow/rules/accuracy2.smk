

rule plot_quilt_mspbwt_by_chunk:
    input:
        truth=rules.collect_truth_gts.output.gt,
        af=rules.collect_truth_gts.output.af,
        vcf=rules.quilt_run_mspbwt.output,
    output:
        rds=os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.{ff}f.mspbwt.{chrom}.{start}.{end}.vcf.gz.rds",
        ),
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy_nipt.R"


rule collect_quilt_mspbwt_accuracy:
    input:
        get_quilt_mspbwt_accuracy_by_chunk
    output:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.{ff}f.mspbwt.{chrom}.vcf.gz.rds.list",
        ),
    conda:
        "../envs/quilt.yaml"
    shell:
        "echo {input} | tr ' ' '\n' > {output}"


rule plot_quilt_mspbwt:
    input:
        truth=rules.collect_truth_gts.output.gt,
        af=rules.collect_truth_gts.output.af,
        vcf=rules.quilt_ligate_mspbwt.output.vcf,
    output:
        rds=os.path.join(
            OUTDIR_SUMMARY, "quilt.nipt.accuracy.mspbwt.refsize{size}.{chrom}.down{depth}x.{ff}f.rds"
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY, "quilt.nipt.accuracy.mspbwt.refsize{size}.{chrom}.down{depth}x.{ff}f.rds.llog"
        ),
    params:
        N="plot_quilt_mspbwt",
        method=config["quilt2"]["method"],
        truth=lambda wildcards: REFPANEL[wildcards.chrom]["truth"],
    resources:
        slots=1,
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy_nipt.R"



rule collapse_mspbwt_accuracy_by_refsize:
    input:
        expand(rules.plot_quilt_mspbwt.output.rds,
               depth=config["coverage"], ff=config["fetalfrac"]
               ,allow_missing=True),
    output:
        rds=os.path.join(
            OUTDIR_SUMMARY, "quilt.nipt.accuracy.mspbwt.refsize{size}.{chrom}.rds"
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY, "quilt.nipt.accuracy.mspbwt.refsize{size}.{chrom}.rds.llog"
        ),
    params:
        N="collapse_accuracy_by_refsize",
    resources:
        slots=1,
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy_collapse.R"



