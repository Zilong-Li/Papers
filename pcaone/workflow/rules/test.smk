DIRTEST = os.path.join(OUTDIR, "tests")


rule calc_test_accuracy:
    input:
        pcaonea=rules.run_pcaone_arnoldi.output.vec,
        pcaoneh=rules.run_pcaone_alg1.output.vec,
        pcaonef=expand(
            rules.run_pcaone_alg2_windows.output.vec,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRTEST, "{data}.acc.k{k}.rds"),
        pdf=os.path.join(DIRTEST, "{data}.acc.k{k}.pdf"),
    log:
        os.path.join(DIRTEST, "{data}.acc.k{k}.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_acc.R"


rule calc_test_log:
    input:
        pcaoneh=rules.run_pcaone_alg1.log,
        pcaonef=expand(
            rules.run_pcaone_alg2_windows.log,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRTEST, "{data}.log.k{k}.rds"),
        pdf=os.path.join(DIRTEST, "{data}.log.k{k}.pdf"),
    log:
        os.path.join(DIRTEST, "{data}.log.k{k}.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_log.R"


rule plot_tests_summary:
    input:
        log=rules.calc_test_log.output.rds,
        acc=rules.calc_test_accuracy.output.rds,
    output:
        pdf=os.path.join(DIRTEST, "{data}.summary.k{k}.pdf"),
    log:
        os.path.join(DIRTEST, "{data}.summary.k{k}.llog"),
    wildcard_constraints:
        data=".+(?<!scrnas)",  # not end with scrna. negative lookbehind assertion. https://stackoverflow.com/questions/46856698/cant-get-this-regex-to-work-for-wildcard-constraints-in-snakemake
    script:
        "../scripts/summary_test.R"


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
        rds=os.path.join(DIRTEST, "{data}.scrnas.acc.k{k}.rds"),
        pdf=os.path.join(DIRTEST, "{data}.scrnas.acc.k{k}.pdf"),
    log:
        os.path.join(DIRTEST, "{data}.scrnas.acc.k{k}.llog"),
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
        rds=os.path.join(DIRTEST, "{data}.scrnas.log.k{k}.rds"),
        pdf=os.path.join(DIRTEST, "{data}.scrnas.log.k{k}.pdf"),
    log:
        os.path.join(DIRTEST, "{data}.scrnas.log.k{k}.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/pcaone_log.R"


rule plot_tests_summary_scrnas:
    input:
        log=rules.get_scrnas_log.output.rds,
        acc=rules.get_scrnas_acc.output.rds,
    output:
        pdf=os.path.join(DIRTEST, "{data}.scrnas.summary.k{k}.pdf"),
    log:
        os.path.join(DIRTEST, "{data}.scrnas.summary.k{k}.llog"),
    script:
        "../scripts/summary_test.R"
