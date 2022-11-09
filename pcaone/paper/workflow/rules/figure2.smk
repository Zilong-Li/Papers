DIRFIG = os.path.join(OUTDIR, "figure2")


rule fig2_dataset1_acc:
    input:
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
        rds=os.path.join(DIRFIG, "A.{data}.acc.k{k}.rds"),
        pdf=os.path.join(DIRFIG, "A.{data}.acc.k{k}.pdf"),
    log:
        os.path.join(DIRFIG, "A.{data}.acc.k{k}.llog"),
    script:
        "../scripts/accuracy_byk.R"


rule fig2_dataset2_acc:
    input:
        vecplink1=rules.run_plink1.output.vec,
        vecplink2=rules.run_plink2.output.vec,
        vecpcaonea=rules.run_pcaone_arnoldi.output.vec,
        vecpcaoneh=rules.run_pcaone_alg1.output.vec,
        vecpcaonef=rules.run_pcaone_alg2.output.vec,
        vecflashpca=rules.run_flashpca.output.vec,
        vecterapca=rules.run_terapca.output.vec,
        vecpropca=rules.run_propca.output.vec,
        logplink1=rules.run_plink1.log,
        logplink2=rules.run_plink2.log,
        logpcaonea=rules.run_pcaone_arnoldi.log,
        logpcaoneh=rules.run_pcaone_alg1.log,
        logpcaonef=rules.run_pcaone_alg2.log,
        logflashpca=rules.run_flashpca.log,
        logterapca=rules.run_terapca.log,
        logpropca=rules.run_propca.log,
    output:
        rds=os.path.join(DIRFIG, "B.{data}.acc.k{k}.rds"),
        pdf=os.path.join(DIRFIG, "B.{data}.acc.k{k}.pdf"),
    log:
        os.path.join(DIRFIG, "B.{data}.acc.k{k}.llog"),
    script:
        "../scripts/accuracy_byk.R"


rule plot_fig2A:
    input:
        expand(
            rules.fig2_dataset1_acc.output.rds,
            data=config["figure2A"]["dataset"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRFIG, "fig2A.k{k}.rds"),
        pdf=os.path.join(DIRFIG, "fig2A.k{k}.acc.pdf"),
        pdf2=os.path.join(DIRFIG, "fig2A.k{k}.time.pdf"),
    log:
        os.path.join(DIRFIG, "fig2A.k{k}.llog"),
    params:
        dataset=config["figure2A"]["dataset"],
    script:
        "../scripts/plot_fig2.R"


rule plot_fig2B:
    input:
        expand(
            rules.fig2_dataset2_acc.output.rds,
            data=config["figure2B"]["dataset"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(DIRFIG, "fig2B.k{k}.rds"),
        pdf=os.path.join(DIRFIG, "fig2B.k{k}.acc.pdf"),
        pdf2=os.path.join(DIRFIG, "fig2B.k{k}.time.pdf"),
    log:
        os.path.join(DIRFIG, "fig2B.k{k}.llog"),
    params:
        dataset=config["figure2B"]["dataset"],
    script:
        "../scripts/plot_fig2.R"


rule plot_fig2:
    input:
        A=rules.plot_fig2A.output.rds,
        B=rules.plot_fig2B.output.rds,
    output:
        pdf=os.path.join(DIRFIG, "fig2.k{k}.pdf"),
        rds=os.path.join(DIRFIG, "fig2.k{k}.rds"),
    log:
        os.path.join(DIRFIG, "fig2.k{k}.llog"),
    script:
        "../scripts/merge_fig2.R"
