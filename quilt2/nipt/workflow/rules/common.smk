import os
import pandas as pd


NIPTFAM = (
    pd.read_csv(os.path.abspath(config["samples_nipt"]), sep="\t", dtype=str)
    .set_index("FamilyID")
    .to_dict(orient="index")
)


# dict : {"chr1": {"vcf": chr1.vcf.gz,"start":1, "end": 99999}, ...}
REFPANEL = (
    pd.read_csv(os.path.abspath(config["genome"]["refpanel"]), sep="\t", dtype=str)
    .set_index("chr")
    .to_dict(orient="index")
)


OUTDIR = "results"
OUTDIR_PANEL = os.path.join(OUTDIR, "refpanels", "")
OUTDIR_TRUTH = os.path.join(OUTDIR, "truth", "")
OUTDIR_QUILT1 = os.path.join(OUTDIR, "quilt1", "")
OUTDIR_QUILT2 = os.path.join(OUTDIR, "quilt2", "")
OUTDIR_SUMMARY = os.path.join(OUTDIR, "summary", "")
OUTDIR_REPORT = os.path.join(OUTDIR, "report", "")

RUN = config["RUN"]


def get_all_results():
    if RUN == "makesample":
        return makesample()
    elif RUN == "gettruth":
        return gettruth()
    elif RUN == "chunkrefpanel":
        return chunkrefpanel()
    elif RUN == "quilt1chunk":
        return quilt1chunk()
    elif RUN == "quilt1acc":
        return quilt1acc()
    elif RUN == "quilt2acc":
        return quilt2acc()
    else:
        raise RuntimeError("this is an invalid RUN!")

def get_chunks_by_chrom(chrom):
    """split chr into chunks given chunksize; return [starts],[ends]"""
    chunksize = config["chunksize"]
    starts, ends = [], []
    if REFPANEL[chrom].get("region"):
        rg = REFPANEL[chrom]["region"].split("-")
        starts = [int(rg[0])]
        ends = [int(rg[1])]
    else:
        s, e = int(REFPANEL[chrom]["start"]), int(REFPANEL[chrom]["end"])
        n = int((e - s) / chunksize) + 1
        if (n - 1) * chunksize == e - s:
            n = n - 1
        for i in range(n):
            ps = chunksize * i + s
            pe = chunksize * (i + 1) + s - 1
            pe = e if pe > e else pe
            starts.append(ps)
            ends.append(pe)
    return starts, ends


def get_glimpse_chunks_in(fn):
    starts, ends = [], []
    with open(fn) as f:
        for row in f:
            """0       chr20   chr20:82590-6074391     chr20:82590-5574162     5491573 1893"""
            tmp = row.strip().split("\t")
            rg = tmp[2].split(":")
            ps, pe = rg[1].split("-")
            starts.append(int(ps))
            ends.append(int(pe))
    return starts, ends


def get_refpanel_chunks(chrom):
    """use quilt_chunk or glimpse_chunk if defined, otherwise split the chromosome by chunks"""
    if REFPANEL[chrom].get("quilt_chunk"):
        return get_glimpse_chunks_in(REFPANEL[chrom]["quilt_chunk"])
    elif REFPANEL[chrom].get("glimpse_chunk"):
        return get_glimpse_chunks_in(REFPANEL[chrom]["glimpse_chunk"])
    else:
        return get_chunks_by_chrom(chrom)



def makesample():
    return expand(
        rules.bamlist.output,
        depth=config["coverage"],
        ff=config["fetalfrac"],
        chrom=config["chroms"],
    )


def gettruth():
    return expand(
        rules.collect_truth_gts.output,
        chrom=config["chroms"],
    )

def chunkrefpanel():
    return expand(
        rules.concat_refpanel_sites_by_region.output,
        size=config["refsize"],
        chrom=config["chroms"],
    )


def quilt1chunk():
    return expand(
        rules.collect_quilt_regular_accuracy.output,
        size=config["refsize"],
        chrom=config["chroms"],
        depth=config["coverage"],
        ff=config["fetalfrac"],
    )


def quilt1acc():
    return expand(
        rules.collapse_regular_accuracy_by_refsize.output,
        size=config["refsize"],
        chrom=config["chroms"],
    )

def quilt2acc():
    return expand(
        rules.collapse_mspbwt_accuracy_by_refsize.output,
        size=config["refsize"],
        chrom=config["chroms"],
    )

def get_samples_in_nipt_fam():
    out = []
    for fam in NIPTFAM.keys():
        out.append(NIPTFAM[fam].get("Kid"))
        out.append(NIPTFAM[fam].get("Mid"))
        out.append(NIPTFAM[fam].get("Fid"))
    return out


def if_use_quilt_map_in_refpanel(wildcards):
    if REFPANEL[wildcards.chrom].get("quilt_map"):
        # varify quilt genetic map file.
        return REFPANEL[wildcards.chrom]["quilt_map"]
    else:
        return "false"


def if_use_af_in_refpanel(wildcards):
    if REFPANEL[wildcards.chrom].get("af"):
        # varify af file. 5 columns: chr pos ref alt af
        return REFPANEL[wildcards.chrom]["af"]
    else:
        return "false"

def get_quilt_regular_output(wildcards):
    starts, ends = get_refpanel_chunks(wildcards.chrom)
    return expand(
        rules.quilt_run_regular.output, zip, start=starts, end=ends, allow_missing=True
    )

def get_quilt_mspbwt_output(wildcards):
    starts, ends = get_refpanel_chunks(wildcards.chrom)
    return expand(
        rules.quilt_run_mspbwt.output, zip, start=starts, end=ends, allow_missing=True
    )

def get_quilt_regular_accuracy_by_chunk(wildcards):
    starts, ends = get_refpanel_chunks(wildcards.chrom)
    return expand(
        rules.plot_quilt_regular_by_chunk.output, zip, start=starts, end=ends, allow_missing=True
    )

def get_quilt_mspbwt_accuracy_by_chunk(wildcards):
    starts, ends = get_refpanel_chunks(wildcards.chrom)
    return expand(
        rules.plot_quilt_mspbwt_by_chunk.output, zip, start=starts, end=ends, allow_missing=True
    )
