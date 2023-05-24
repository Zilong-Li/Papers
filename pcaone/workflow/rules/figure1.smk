DIRFIG = os.path.join(OUTDIR, "figure1")


rule fig1_accuracy:
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
        rds=os.path.join(DIRFIG, "{data}.acc.k{k}.rds"),
        pdf=os.path.join(DIRFIG, "{data}.acc.k{k}.pdf"),
    log:
        os.path.join(DIRFIG, "{data}.acc.k{k}.llog"),
    script:
        "../scripts/accuracy_byk.R"


rule plot_fig1A:
    input:
        pcaonea=expand(rules.run_pcaone_arnoldi.log, k=PCS, allow_missing=True),
        pcaoneh=expand(rules.run_pcaone_alg1.log, k=PCS, allow_missing=True),
        pcaonef=expand(rules.run_pcaone_alg2.log, k=PCS, allow_missing=True),
        terapca=expand(rules.run_terapca.log, k=PCS, allow_missing=True),
        propca=expand(rules.run_propca.log, k=PCS, allow_missing=True),
        flashpca=expand(rules.run_flashpca.log, k=PCS, allow_missing=True),
        rds=expand(rules.fig1_accuracy.output.rds, k=PCS, allow_missing=True),
        # vecplink1=expand(rules.run_plink1.output.vec, k=PCS, allow_missing=True),
        vecfull=expand(rules.run_pcaone_full.output.vec, k=PCS, allow_missing=True),
        vecpcaoneh=expand(rules.run_pcaone_alg1.output.vec, k=PCS, allow_missing=True),
        vecpcaonef=expand(rules.run_pcaone_alg2.output.vec, k=PCS, allow_missing=True),
    output:
        rds=os.path.join(DIRFIG, "fig1A.{data}.rds"),
        pdf=os.path.join(DIRFIG, "fig1A.{data}.pdf"),
        pdf2=os.path.join(DIRFIG, "fig1A.{data}.time.pdf"),
    log:
        os.path.join(DIRFIG, "fig1A.{data}.llog"),
    params:
        pcs=PCS,
    script:
        "../scripts/plot_fig1A.R"


rule plot_fig1B:
    input:
        vecfull=expand(rules.run_pcaone_full.output.vec, k=["2"], allow_missing=True),
        vecfastpca=expand(rules.run_fastpca.output.vec, k=["2"], allow_missing=True),
        vecplink2=expand(rules.run_plink2.output.vec, k=["2"], allow_missing=True),
        vecpcaonea=expand(
            rules.run_pcaone_arnoldi.output.vec, k=["2"], allow_missing=True
        ),
        vecpcaoneh=expand(rules.run_pcaone_alg1.output.vec, k=["2"], allow_missing=True),
        vecpcaonef=expand(rules.run_pcaone_alg2.output.vec, k=["2"], allow_missing=True),
        vecflashpca=expand(rules.run_flashpca.output.vec, k=["2"], allow_missing=True),
        vecterapca=expand(rules.run_terapca.output.vec, k=["2"], allow_missing=True),
        vecpropca=expand(rules.run_propca.output.vec, k=["2"], allow_missing=True),
    output:
        pdf=os.path.join(DIRFIG, "fig1B.{data}.pdf"),
    log:
        os.path.join(DIRFIG, "fig1B.{data}.llog"),
    script:
        "../scripts/plot_fig1B.R"


rule plot_fig1:
    input:
        A=rules.plot_fig1A.output.pdf,
        B=rules.plot_fig1B.output.pdf,
    output:
        os.path.join(DIRFIG, "fig1.{data}.log"),
    shell:
        "echo make figure1 {input.A} and {input.B} done > {output}"
