DIRSUPP = os.path.join(OUTDIR, "supplementary")


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
