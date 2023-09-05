OUTDIR_TRUTH = os.path.join(OUTDIR, "truth", "")

rule collect_truth_gts:
    """would be better to use sites in subrefs"""
    input:
        sites=rules.get_pos_from_refpanel.output.sites,
    output:
        gt=os.path.join(OUTDIR_TRUTH, "truth.gts.{chrom}.txt"),
        af=os.path.join(OUTDIR_TRUTH, "af.input.panel.{chrom}.txt"),
        tmp=temp(os.path.join(OUTDIR_TRUTH, "af.input.panel.{chrom}.txt.tmp")),
        tmp2=temp(os.path.join(OUTDIR_TRUTH, "truth.gts.{chrom}.txt.tmp")),
    log:
        os.path.join(OUTDIR_TRUTH, "truth.gts.{chrom}.log"),
    params:
        samples=",".join(SAMPLES.keys()),
        truth=lambda wildcards: REFPANEL[wildcards.chrom]["truth"],
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
            bcftools +fill-tags {input.sites} -- -t AF |  bcftools query -f '{params.ql1}' > {output.tmp}; \
        fi
        awk '{params.awk}' <(bcftools query -f '{params.ql0}' {input.sites}) {output.tmp} >{output.af}
        bcftools view -s {params.samples} {params.truth} | bcftools query -f '{params.ql2}' > {output.tmp2}
        HEADER="ID,"{params.samples}
        echo $HEADER | tr ',' '\t' > {output.gt}
        awk '{params.awk2}' <(bcftools query -f '{params.ql0}' {input.sites}) {output.tmp2} >> {output.gt}
        ) &> {log}
        """


rule collect_beagle_accuracy:
    input:
        vcf=rules.beagle41_by_chrom.output.vcf,
        truth=lambda wildcards: REFPANEL[wildcards.chrom]["truth"],
        af=rules.collect_truth_gts.output.af,
    output:
        rds=os.path.join(OUTDIR, "beagle4.1", "accuracy.down{depth}x.{chrom}.rds"),
    params:
        samples=",".join(SAMPLES.keys()),
    script:
        "../scripts/accuracy_single.R"


rule collect_stitch_accuracy:
    input:
        vcf=rules.stitch_by_chrom.output.vcf,
        truth=lambda wildcards: REFPANEL[wildcards.chrom]["truth"],
        af=rules.collect_truth_gts.output.af,
    output:
        rds=os.path.join(OUTDIR, "stitch", "accuracy.down{depth}x.{chrom}.rds"),
    params:
        samples=",".join(SAMPLES.keys()),
    script:
        "../scripts/accuracy_single.R"


rule collect_phaseless_accuracy:
    input:
        vcf=rules.phaseless_impute_by_chrom.output.vcf,
        truth=lambda wildcards: REFPANEL[wildcards.chrom]["truth"],
        af=rules.collect_truth_gts.output.af,
    output:
        rds=os.path.join(OUTDIR, "phaseless", "accuracy.down{depth}x.{chrom}.rds"),
    params:
        samples=",".join(SAMPLES.keys()),
    script:
        "../scripts/accuracy_single.R"



