DIRSUPP = os.path.join(OUTDIR, "supplementary")


rule calc_mev_ci:
    input:
        vecfull=rules.run_pcaone_full.output.vec,
        vecplink2=rules.run_plink2.output.vec,
        vecpcaonea=rules.run_pcaone_arnoldi.output.vec,
        vecpcaoneh=rules.run_pcaone_alg1.output.vec,
        vecpcaonef=rules.run_pcaone_alg2.output.vec,
        vecflashpca=rules.run_flashpca.output.vec,
        vecterapca=rules.run_terapca.output.vec,
        vecpropca=rules.run_propca.output.vec,
        logplink2=rules.run_plink2.log,
        logpcaonea=rules.run_pcaone_arnoldi.log,
        logpcaoneh=rules.run_pcaone_alg1.log,
        logpcaonef=rules.run_pcaone_alg2.log,
        logflashpca=rules.run_flashpca.log,
        logterapca=rules.run_terapca.log,
        logpropca=rules.run_propca.log,
    output:
        rds=os.path.join(DIRSUPP, "{data}.acc.k{k}.rds"),
    log:
        os.path.join(DIRSUPP, "{data}.acc.k{k}.llog"),
    script:
        "../scripts/comp-mev-ci.R"


rule collapse_mev_ci:
    input:
        rds=expand(rules.calc_mev_ci.output.rds, data=DATASETS, allow_missing=True),
    output:
        rds=os.path.join(DIRSUPP, scenario, "mev.ci.k{k}.rds"),
    params:
        data=DATASETS,
    log:
        os.path.join(DIRSUPP, scenario, "mev.ci.k{k}.llog"),
    script:
        "../scripts/comp-mev-ci.R"


rule calc_binary_accuracy:
    input:
        pcaonea=rules.run_binary_pcaone_arnoldi.output.vec,
        pcaoneh=rules.run_binary_pcaone_alg1.output.vec,
        pcaonef=expand(
            rules.run_binary_pcaone_alg2.output.vec,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRSUPP, "{data}.acc.k{k}.binary.rds"),
        pdf=os.path.join(DIRSUPP, "{data}.acc.k{k}.binary.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.acc.k{k}.binary.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_acc.R"


rule calc_bgen_accuracy:
    input:
        pcaonea=rules.run_bgen_pcaone_arnoldi.output.vec,
        pcaoneh=rules.run_bgen_pcaone_alg1.output.vec,
        pcaonef=expand(
            rules.run_bgen_pcaone_alg2.output.vec,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRSUPP, "{data}.acc.k{k}.bgen.rds"),
        pdf=os.path.join(DIRSUPP, "{data}.acc.k{k}.bgen.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.acc.k{k}.bgen.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_acc.R"


rule calc_binary_log:
    input:
        pcaoneh=rules.run_binary_pcaone_alg1.log,
        pcaonef=expand(
            rules.run_binary_pcaone_alg2.log,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRSUPP, "{data}.log.k{k}.binary.rds"),
        pdf=os.path.join(DIRSUPP, "{data}.log.k{k}.binary.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.log.k{k}.binary.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_log.R"


rule calc_bgen_log:
    input:
        pcaoneh=rules.run_bgen_pcaone_alg1.log,
        pcaonef=expand(
            rules.run_bgen_pcaone_alg2.log,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRSUPP, "{data}.log.k{k}.bgen.rds"),
        pdf=os.path.join(DIRSUPP, "{data}.log.k{k}.bgen.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.log.k{k}.bgen.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_log.R"


rule plot_binary_summary:
    input:
        log=rules.calc_binary_log.output.rds,
        acc=rules.calc_binary_accuracy.output.rds,
    output:
        pdf=os.path.join(DIRSUPP, "{data}.summary.k{k}.binary.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.summary.k{k}.binary.llog"),
    wildcard_constraints:
        data=".+(?<!scrnas)",  # not end with scrna. negative lookbehind assertion. https://stackoverflow.com/questions/46856698/cant-get-this-regex-to-work-for-wildcard-constraints-in-snakemake
    script:
        "../scripts/summary_test.R"


rule plot_bgen_summary:
    input:
        log=rules.calc_bgen_log.output.rds,
        acc=rules.calc_bgen_accuracy.output.rds,
    output:
        pdf=os.path.join(DIRSUPP, "{data}.summary.k{k}.bgen.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.summary.k{k}.bgen.llog"),
    wildcard_constraints:
        data=".+(?<!scrnas)",  # not end with scrna. negative lookbehind assertion. https://stackoverflow.com/questions/46856698/cant-get-this-regex-to-work-for-wildcard-constraints-in-snakemake
    script:
        "../scripts/summary_test.R"


rule calc_stop_criteria:
    input:
        full=rules.run_pcaone_full.output.vec,
        pcaoneh=rules.run_pcaone_alg1.output.vec,
        pcaonef=expand(
            rules.run_pcaone_alg2_windows.output.vec,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRSUPP, scenario, "{data}.acc.k{k}.bfile.rds"),
    log:
        os.path.join(DIRSUPP, scenario, "{data}.acc.k{k}.bfile.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/comp-stop-criteria.R"


rule calc_bfile_accuracy:
    input:
        pcaonea=rules.run_pcaone_arnoldi.output.vec,
        pcaoneh=rules.run_pcaone_alg1.output.vec,
        pcaonef=expand(
            rules.run_pcaone_alg2_windows.output.vec,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRSUPP, "{data}.acc.k{k}.bfile.rds"),
        pdf=os.path.join(DIRSUPP, "{data}.acc.k{k}.bfile.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.acc.k{k}.bfile.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_acc.R"


rule calc_bfile_log:
    input:
        pcaoneh=rules.run_pcaone_alg1.log,
        pcaonef=expand(
            rules.run_pcaone_alg2_windows.log,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRSUPP, "{data}.log.k{k}.bfile.rds"),
        pdf=os.path.join(DIRSUPP, "{data}.log.k{k}.bfile.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.log.k{k}.bfile.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_log.R"


rule plot_bfile_summary:
    input:
        log=rules.calc_bfile_log.output.rds,
        acc=rules.calc_bfile_accuracy.output.rds,
    output:
        pdf=os.path.join(DIRSUPP, "{data}.summary.k{k}.bfile.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.summary.k{k}.bfile.llog"),
    wildcard_constraints:
        data=".+(?<!scrnas)",  # not end with scrna. negative lookbehind assertion. https://stackoverflow.com/questions/46856698/cant-get-this-regex-to-work-for-wildcard-constraints-in-snakemake
    script:
        "../scripts/summary_test.R"


rule collect_bfile_summary:
    input:
        log=expand(rules.calc_bfile_log.output.rds, data=DATASETS, allow_missing=True),
        acc=expand(
            rules.calc_bfile_accuracy.output.rds, data=DATASETS, allow_missing=True
        ),
    output:
        rds=os.path.join(DIRSUPP, scenario, "summary.k{k}.bfile.rds"),
    params:
        data=DATASETS,
    log:
        os.path.join(DIRSUPP, scenario, "summary.k{k}.bfile.llog"),
    script:
        "../scripts/comp-params.R"


rule collapse_bfile_summary:
    input:
        rds=expand(rules.collect_bfile_summary.output.rds, k=PCS),
    output:
        rds=os.path.join(DIRSUPP, scenario, "summary.bfile.rds"),
        pdf1=os.path.join(DIRSUPP, scenario, "summary.bfile.bysnps.pdf"),
        pdf2=os.path.join(DIRSUPP, scenario, "summary.bfile.bysamples.pdf"),
    params:
        k=PCS,
    log:
        os.path.join(DIRSUPP, scenario, "summary.bfile.llog"),
    script:
        "../scripts/comp-params.R"


rule collapse_binary_summary:
    input:
        log=expand(rules.calc_binary_log.output.rds, k=PCS, allow_missing=True),
        acc=expand(rules.calc_binary_accuracy.output.rds, k=PCS, allow_missing=True),
    output:
        rds=os.path.join(DIRSUPP, scenario, "summary.{data}.binary.rds"),
        pdf=os.path.join(DIRSUPP, scenario, "summary.{data}.binary.pdf"),
    params:
        k=PCS,
    log:
        os.path.join(DIRSUPP, scenario, "summary.{data}.binary.llog"),
    script:
        "../scripts/comp-params.R"


rule collapse_bgen_summary:
    input:
        log=expand(rules.calc_bgen_log.output.rds, k=PCS, allow_missing=True),
        acc=expand(rules.calc_bgen_accuracy.output.rds, k=PCS, allow_missing=True),
    output:
        rds=os.path.join(DIRSUPP, scenario, "summary.{data}.bgen.rds"),
        pdf=os.path.join(DIRSUPP, scenario, "summary.{data}.bgen.pdf"),
    params:
        k=PCS,
    log:
        os.path.join(DIRSUPP, scenario, "summary.{data}.bgen.llog"),
    script:
        "../scripts/comp-params.R"


rule get_scrnas_acc:
    input:
        pcaonea=rules.run_scrnas_pcaone_arnoldi.output.vec,
        pcaoneh=rules.run_scrnas_pcaone_alg1.output.vec,
        pcaonef=expand(
            rules.run_scrnas_pcaone_alg2.output.vec,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRSUPP, "{data}.scrnas.acc.k{k}.rds"),
        pdf=os.path.join(DIRSUPP, "{data}.scrnas.acc.k{k}.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.scrnas.acc.k{k}.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_acc.R"


rule get_scrnas_log:
    input:
        pcaoneh=rules.run_scrnas_pcaone_alg1.log,
        pcaonef=expand(
            rules.run_scrnas_pcaone_alg2.log,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRSUPP, "{data}.scrnas.log.k{k}.rds"),
        pdf=os.path.join(DIRSUPP, "{data}.scrnas.log.k{k}.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.scrnas.log.k{k}.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_log.R"


rule plot_summary_scrnas:
    input:
        log=rules.get_scrnas_log.output.rds,
        acc=rules.get_scrnas_acc.output.rds,
    output:
        pdf=os.path.join(DIRSUPP, "{data}.scrnas.summary.k{k}.pdf"),
    log:
        os.path.join(DIRSUPP, "{data}.scrnas.summary.k{k}.llog"),
    script:
        "../scripts/summary_test.R"
