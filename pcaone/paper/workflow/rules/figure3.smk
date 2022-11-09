DIRFIG = os.path.join(OUTDIR, "figure3")


rule fig3_table:
    input:
        pcaonea=rules.run_pcaone_arnoldi.output.vec,
        pcaoneh=rules.run_pcaone_alg1.output.vec,
        pcaonef=rules.run_pcaone_alg2.output.vec,
        terapca=rules.run_terapca.output.vec,
    output:
        rds=os.path.join(DIRFIG, "{data}.k{k}.rds"),
        pdf=os.path.join(DIRFIG, "{data}.k{k}.pdf"),
    log:
        os.path.join(DIRFIG, "{data}.k{k}.llog"),
    script:
        "../scripts/summary_fig3.R"
