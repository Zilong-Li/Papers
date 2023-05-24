import os
import sys

OUTDIR = "results"

valid_scenarios = [
    "figure1",
    "figure2",
    "figure3",
    "figure4",
    "figure2A",
    "figure2B",
    "suppl_rnaseq",
    "suppl_bgen",
    "suppl_bfile",
    "suppl_ukb",
    "test_ukb",
    "test_scrnas",
    "stop_criteria",
    "suppl_mev_ci",
]

scenario = config["scenario"]

if scenario not in valid_scenarios:
    raise RuntimeError(
        "please specify vaild scenario to run! valid scenarios are:\n"
        + str(valid_scenarios)
    )

DATASETS = config[scenario]["dataset"]
PCS = config[scenario]["pcs"]
THREADS = config[scenario]["threads"]

# global env
# Programs
CSV2BIN = config["csv2bin"]
SUMR = config["sumr"]
HALKO = config["jlhalko"]
JULIA = config["julia"]
TIME = config["time"]
PCAONE = config["pcaone"]
TERAPCA = config["terapca"]
PROPCA = config["propca"]
PLINK1 = config["plink1"]
PLINK2 = config["plink2"]
FLASHPCA = config["flashpca"]
SMARTPCA = config["smartpca"]


wildcard_constraints:
    k="\d+",
    w="\d+",


def get_all_results():
    if scenario == "figure1":
        return get_figure1_plot()
    elif scenario == "figure2":
        return get_figure2_plot()
    elif scenario == "figure3":
        return get_figure3_summary()
    elif scenario == "figure4":
        return get_figure4_plot()
    elif scenario == "figure2A":
        return get_figure2A_plot()
    elif scenario == "figure2B":
        return get_figure2B_plot()
    elif scenario == "suppl_rnaseq":
        return get_suppl_rnaseq_plot()
    elif scenario == "suppl_bgen":
        return get_suppl_bgen_plot()
    elif scenario == "test_scrnas":
        return get_test_summary_scrnas()
    elif scenario == "stop_criteria":
        return get_stop_criteira()
    elif scenario == "suppl_mev_ci":
        return get_suppl_mev_ci()
    elif scenario in ["suppl_bfile", "suppl_ukb"]:
        return rules.collapse_bfile_summary.output
    else:
        raise RuntimeError(
            "please specify vaild scenario to run! but you won't see me.\n"
        )


def get_suppl_mev_ci():
    return expand(
        rules.collapse_mev_ci.output,
        k=PCS,
    )


def get_suppl_rnaseq_plot():
    return (
        expand(
            rules.collapse_binary_summary.output,
            data=DATASETS,
        ),
        expand(rules.make_pca_plot.output, data=DATASETS),
    )


def get_suppl_bgen_plot():
    return expand(rules.collapse_bgen_summary.output, data=DATASETS,), expand(
        rules.collect_bfile_summary.output,
        k=PCS,
    )


def get_stop_criteira():
    return expand(
        rules.calc_stop_criteria.output,
        k=PCS,
        data=DATASETS,
    )


def get_test_summary_scrnas():
    return expand(
        rules.plot_summary_scrnas.output,
        k=PCS,
        data=DATASETS,
    )


def get_suppl_bfile_summary():
    return expand(
        rules.collect_bfile_summary.output,
        k=PCS,
    )


def get_test_accuracy():
    return expand(
        rules.calc_test_accuracy.output,
        k=PCS,
        data=DATASETS,
    )


def get_test_log():
    return expand(
        rules.calc_test_log.output,
        k=PCS,
        data=DATASETS,
    )


def get_figure4_plot():
    return expand(
        rules.get_scrnas_brain.output,
        k=PCS,
        data=["brain"],
    )


def get_figure3_summary():
    return expand(
        rules.fig3_table.output,
        k=PCS,
        data=DATASETS,
    )


def get_figure2_plot():
    return expand(
        rules.plot_fig2.output,
        k=PCS,
    )


def get_figure2A_plot():
    return expand(
        rules.plot_fig2A.output,
        k=PCS,
        data=config["figure2A"]["dataset"],
    )


def get_figure2B_plot():
    return expand(
        rules.plot_fig2B.output,
        k=PCS,
        data=config["figure2B"]["dataset"],
    )


def get_figure1_plot():
    return expand(
        rules.plot_fig1.output,
        k=PCS,
        data=DATASETS,
    )


def get_results_plink1():
    return expand(
        rules.run_plink1.output,
        k=[100],
        data=DATASETS,
    )


def get_results_plink2():
    return expand(
        rules.run_plink2.output,
        k=PCS,
        data=DATASETS,
    )


def get_results_flashpca():
    return expand(
        rules.run_flashpca.output,
        k=PCS,
        data=DATASETS,
    )


def get_results_pcaone_arnoldi():
    return expand(
        rules.run_pcaone_arnoldi.output,
        k=PCS,
        data=DATASETS,
    )


def get_results_pcaone_alg1():
    return expand(
        rules.run_pcaone_alg1.output,
        k=PCS,
        data=DATASETS,
    )


def get_results_pcaone_alg2():
    return expand(
        rules.run_pcaone_alg2.output,
        k=PCS,
        data=DATASETS,
    )


def get_results_pcaone_alg2_windows():
    return expand(
        rules.run_pcaone_alg2_windows.output,
        k=PCS,
        data=DATASETS,
        w=config[scenario]["windows"],
    )


def get_results_terapca():
    return expand(
        rules.run_terapca.output,
        k=PCS,
        data=DATASETS,
    )


def get_results_propca():
    return expand(
        rules.run_propca.output,
        k=PCS,
        data=DATASETS,
    )


def get_results_fastpca():
    return expand(
        rules.run_fastpca.output,
        k=PCS,
        data=DATASETS,
    )


def get_input_dataset(wildcards):
    return {
        "bed": config[wildcards.data]["bed"],
        "bim": config[wildcards.data]["bim"],
        "fam": config[wildcards.data]["fam"],
    }


rule run_pcaone_full:
    input:
        unpack(get_input_dataset),
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.full.eigvecs"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.full.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-8],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        # touch {output}
        {TIME} -v {PCAONE} --bfile {params.bfile} -k 100 -n {threads} -o {params.out} --svd 3 --verbose &> {log}
        """


rule run_pcaone_arnoldi:
    input:
        unpack(get_input_dataset),
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.a.k{k}.eigvecs"),
        val=os.path.join(OUTDIR, "{data}", "pcaone.a.k{k}.eigvals"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.a.k{k}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-8],
        pcaonea=config[scenario]["pcaonea"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --bfile {params.bfile} -k {wildcards.k} -n {threads} -o {params.out} {params.pcaonea} --verbose &> {log}
        """


rule run_pcaone_alg1:
    input:
        unpack(get_input_dataset),
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.h.k{k}.eigvecs"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.h.k{k}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-8],
        pcaoneh=config[scenario]["pcaoneh"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --bfile {params.bfile} -k {wildcards.k} -n {threads} -o {params.out} --verbose {params.pcaoneh} &> {log}
        """


rule run_pcaone_alg2:
    input:
        unpack(get_input_dataset),
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.f.k{k}.eigvecs"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.f.k{k}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-8],
        pcaonef=config[scenario]["pcaonef"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --bfile {params.bfile} -k {wildcards.k} -n {threads} -o {params.out} --verbose {params.pcaonef} &> {log}
        """


rule run_pcaone_alg2_windows:
    input:
        unpack(get_input_dataset),
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.f.w{w}.k{k}.eigvecs"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.f.w{w}.k{k}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-8],
        pcaonef=config[scenario]["pcaonef"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --bfile {params.bfile} -k {wildcards.k} -w {wildcards.w} -n {threads} -o {params.out} --verbose {params.pcaonef} &> {log}
        """


rule run_terapca:
    input:
        unpack(get_input_dataset),
    output:
        vec=os.path.join(OUTDIR, "{data}", "terapca.k{k}_singularVectors.txt"),
        val=os.path.join(OUTDIR, "{data}", "terapca.k{k}_singularValues.txt"),
    log:
        os.path.join(OUTDIR, "{data}", "terapca.k{k}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-20],
        terapca=config[scenario]["terapca"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {TERAPCA} -bfile {params.bfile} -nsv {wildcards.k} -prefix {params.out} {params.terapca} &> {log}
        """


rule run_propca:
    input:
        unpack(get_input_dataset),
    output:
        vec=os.path.join(OUTDIR, "{data}", "propca.k{k}projections.txt"),
        val=os.path.join(OUTDIR, "{data}", "propca.k{k}evals.txt"),
    log:
        os.path.join(OUTDIR, "{data}", "propca.k{k}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[1][:-9],
        propca=config[scenario]["propca"],
    threads: THREADS
    shell:
        """
        {TIME} -v {PROPCA}  -g {params.bfile}  -k {wildcards.k} -nt {threads} -o {params.out} {params.propca} &> {log}
        """


rule run_fastpca:
    input:
        unpack(get_input_dataset),
    output:
        par=os.path.join(OUTDIR, "{data}", "fastpca.k{k}.par"),
        vec=os.path.join(OUTDIR, "{data}", "fastpca.k{k}.eigvecs"),
        val=os.path.join(OUTDIR, "{data}", "fastpca.k{k}.eigvals"),
    log:
        os.path.join(OUTDIR, "{data}", "fastpca.k{k}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-20],
    threads: THREADS
    shell:
        """
        echo "genotypename: {input.bed}
        snpname: {input.bim}
        indivname: {input.fam}
        evecoutname: {output.vec}
        fastmode: YES
        numthreads: {threads}
        numoutevec: {wildcards.k}" > {output.par}
        {TIME} -v {SMARTPCA} -p {output.par} &> {log} || true
        if [ ! -f {output.vec} ]; then touch {output}; fi
        head -1 {output.vec}| perl -lane 'print join("\\n", @F[1..$#F])' >{output.val}
        """


rule run_flashpca:
    input:
        unpack(get_input_dataset),
    output:
        vec=os.path.join(OUTDIR, "{data}", "flashpca.k{k}.evecs.txt"),
        val=os.path.join(OUTDIR, "{data}", "flashpca.k{k}.evals.txt"),
        pc=os.path.join(OUTDIR, "{data}", "flashpca.k{k}.pc.txt"),
        load=os.path.join(OUTDIR, "{data}", "flashpca.k{k}.loadings.txt"),
    log:
        os.path.join(OUTDIR, "{data}", "flashpca.k{k}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        flashpca=config[scenario]["flashpca"],
    threads: 1
    shell:
        """
        {TIME} -v {FLASHPCA} --bfile {params.bfile} -d {wildcards.k} --outload {output.load} --outpc {output.pc} --outvec {output.vec} --outval {output.val} {params.flashpca} &> {log}
        """


rule run_plink1:
    input:
        unpack(get_input_dataset),
    output:
        vec=os.path.join(OUTDIR, "{data}", "plink1.k{k}.eigenvec"),
        val=os.path.join(OUTDIR, "{data}", "plink1.k{k}.eigenval"),
    log:
        os.path.join(OUTDIR, "{data}", "plink1.k{k}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-9],
    threads: THREADS
    shell:
        """
        {TIME} -v {PLINK1} --bfile {params.bfile} --pca {wildcards.k} --out {params.out} --threads {threads} &> {log}  || true
        if [ ! -f {output.vec} ]; then touch {output}; fi
        """


rule run_plink2:
    input:
        unpack(get_input_dataset),
    output:
        vec=os.path.join(OUTDIR, "{data}", "plink2.k{k}.eigenvec"),
        val=os.path.join(OUTDIR, "{data}", "plink2.k{k}.eigenval"),
    log:
        os.path.join(OUTDIR, "{data}", "plink2.k{k}.llog"),
    params:
        bfile=lambda wildcards, input: input[0][:-4],
        out=lambda wildcards, output: output[0][:-9],
    threads: THREADS
    shell:
        """
        {TIME} -v {PLINK2} --bfile {params.bfile} --pca approx {wildcards.k} --out {params.out} --threads {threads} &> {log} || true
        if [ ! -f {output.vec} ]; then touch {output}; fi
        """


rule run_binary_pcaone_arnoldi:
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.a.k{k}.binary.eigvecs"),
        val=os.path.join(OUTDIR, "{data}", "pcaone.a.k{k}.binary.eigvals"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.a.k{k}.binary.llog"),
    params:
        binary=lambda wildcards: config[wildcards.data]["binary"],
        out=lambda wildcards, output: output[0][:-8],
        pcaone=config[scenario]["pcaonea"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --binary {params.binary} -k {wildcards.k} -n {threads} -o {params.out} --verbose {params.pcaone} &> {log}
        """


rule run_binary_pcaone_alg1:
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.h.k{k}.binary.eigvecs"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.h.k{k}.binary.llog"),
    params:
        binary=lambda wildcards: config[wildcards.data]["binary"],
        out=lambda wildcards, output: output[0][:-8],
        pcaone=config[scenario]["pcaoneh"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --binary {params.binary} -k {wildcards.k} -n {threads} -o {params.out} --verbose {params.pcaone} &> {log}
        """


rule run_binary_pcaone_alg2:
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.f.k{k}.binary.eigvecs"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.f.k{k}.binary.llog"),
    params:
        binary=lambda wildcards: config[wildcards.data]["binary"],
        out=lambda wildcards, output: output[0][:-8],
        pcaone=config[scenario]["pcaonef"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --binary {params.binary} -k {wildcards.k} -n {threads} -o {params.out} --verbose {params.pcaone} &> {log}
        """


rule run_bgen_pcaone_arnoldi:
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.a.k{k}.bgen.eigvecs"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.a.k{k}.bgen.llog"),
    params:
        bgen=lambda wildcards: config[wildcards.data]["bgen"],
        out=lambda wildcards, output: output[0][:-8],
        pcaone=config[scenario]["pcaonea"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --bgen {params.bgen} -k {wildcards.k} -n {threads} -o {params.out} --verbose {params.pcaone} &> {log}
        """


rule run_bgen_pcaone_alg1:
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.h.k{k}.bgen.eigvecs"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.h.k{k}.bgen.llog"),
    params:
        bgen=lambda wildcards: config[wildcards.data]["bgen"],
        out=lambda wildcards, output: output[0][:-8],
        pcaone=config[scenario]["pcaoneh"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --bgen {params.bgen} -k {wildcards.k} -n {threads} -o {params.out} --verbose {params.pcaone} &> {log}
        """


rule run_bgen_pcaone_alg2:
    output:
        vec=os.path.join(OUTDIR, "{data}", "pcaone.f.w{w}.k{k}.bgen.eigvecs"),
    log:
        os.path.join(OUTDIR, "{data}", "pcaone.f.w{w}.k{k}.bgen.llog"),
    params:
        bgen=lambda wildcards: config[wildcards.data]["bgen"],
        out=lambda wildcards, output: output[0][:-8],
        pcaone=config[scenario]["pcaonef"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --bgen {params.bgen} -k {wildcards.k} -n {threads} -o {params.out} --verbose {params.pcaone} &> {log}
        """


rule run_scrnas_pcaone_arnoldi:
    output:
        vec=os.path.join(OUTDIR, "scrnas", "{data}", "pcaone.a.k{k}.eigvecs"),
    log:
        os.path.join(OUTDIR, "scrnas", "{data}", "pcaone.a.k{k}.llog"),
    params:
        csv=lambda wildcards: config[wildcards.data]["csv"],
        out=lambda wildcards, output: output[0][:-8],
        pcaonea=config[scenario]["pcaonea"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --csv {params.csv} -k {wildcards.k} --cpmed -n {threads} -o {params.out} --verbose {params.pcaonea} &> {log}
        """


rule run_scrnas_pcaone_alg1:
    output:
        vec=os.path.join(OUTDIR, "scrnas", "{data}", "pcaone.h.k{k}.eigvecs"),
    log:
        os.path.join(OUTDIR, "scrnas", "{data}", "pcaone.h.k{k}.llog"),
    params:
        csv=lambda wildcards: config[wildcards.data]["csv"],
        out=lambda wildcards, output: output[0][:-8],
        pcaoneh=config[scenario]["pcaoneh"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --csv {params.csv} -k {wildcards.k} --cpmed -n {threads} -o {params.out} --verbose {params.pcaoneh} &> {log}
        """


rule run_scrnas_pcaone_alg2:
    output:
        vec=os.path.join(OUTDIR, "scrnas", "{data}", "pcaone.f.w{w}.k{k}.eigvecs"),
    log:
        os.path.join(OUTDIR, "scrnas", "{data}", "pcaone.f.w{w}.k{k}.llog"),
    params:
        csv=lambda wildcards: config[wildcards.data]["csv"],
        out=lambda wildcards, output: output[0][:-8],
        pcaonef=config[scenario]["pcaonef"],
    threads: THREADS
    shell:
        """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        {TIME} -v {PCAONE} --csv {params.csv} -k {wildcards.k} --cpmed -w {wildcards.w} --no-shuffle -n {threads} -o {params.out} --verbose {params.pcaonef} &> {log}
        """


rule run_scrnas_online_sumr:
    output:
        os.path.join(
            OUTDIR,
            "scrnas",
            "{data}",
            "sumr",
            "Feature_LogCPMEDMeans.csv",
        ),
    log:
        os.path.join(OUTDIR, "scrnas", "{data}", "sumr", "summary.llog"),
    params:
        zst=lambda wildcards: config[wildcards.data]["bin"],
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    threads: THREADS
    shell:
        """
        (
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        export JULIA_NUM_THREADS={threads}
        {TIME} -v {JULIA} {SUMR} --binfile {params.zst} --outdir {params.outdir}
        ) &> {log}
        """


rule run_scrnas_online_halko:
    input:
        rules.run_scrnas_online_sumr.output,
    output:
        vec=os.path.join(
            OUTDIR,
            "scrnas",
            "{data}",
            "onlinepca.halko.k{k}.niter3",
            "Eigen_vectors.csv",
        ),
    log:
        os.path.join(OUTDIR, "scrnas", "{data}", "onlinepca.halko.k{k}.niter3", "llog"),
    params:
        zst=lambda wildcards: config[wildcards.data]["bin"],
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        sumdir=lambda wildcards, input: os.path.dirname(input[0]),
    threads: THREADS
    shell:
        """
        (
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        export JULIA_NUM_THREADS={threads}
        cper=$(cat {params.sumdir}/Sample_NoCounts.csv |R -s -e 'd=read.table("stdin");cat(median(d[,1]))')
        {TIME} -v {JULIA} {HALKO} --input {params.zst} --outdir {params.outdir} --scale log --cper $cper --dim {wildcards.k} --niter 3 --rowmeanlist {params.sumdir}/Feature_LogCPMEDMeans.csv --colsumlist {params.sumdir}/Sample_NoCounts.csv
        ) &> {log}
        """
