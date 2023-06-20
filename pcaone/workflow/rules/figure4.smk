
DIRFIG = os.path.join(OUTDIR, "figure4")

rule get_scrnas_brain:
    input:
        pcaonea=rules.run_csv_pcaone_arnoldi.output.vec,
        logpcaonea=rules.run_csv_pcaone_arnoldi.log,
        pcaoneh=rules.run_csv_pcaone_alg1.output.vec,
        logpcaoneh=rules.run_csv_pcaone_alg1.log,
        pcaonef=expand(
            rules.run_csv_pcaone_alg2.output.vec,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
        logpcaonef=expand(
            rules.run_csv_pcaone_alg2.log,
            w=config[scenario]["windows"],
            allow_missing=True,
        ),
        onlinepca=rules.run_scrnas_online_halko.output.vec,
        logonlinepca=rules.run_scrnas_online_halko.log,
    output:
        rds=os.path.join(DIRFIG, "fig4.{data}.k{k}.rds"),
        pdf=os.path.join(DIRFIG, "fig4.{data}.k{k}.pdf"),
    log:
        os.path.join(DIRFIG, "fig4.{data}.k{k}.llog"),
    params:
        windows=config[scenario]["windows"],
    script:
        "../scripts/summary_fig4.R"
