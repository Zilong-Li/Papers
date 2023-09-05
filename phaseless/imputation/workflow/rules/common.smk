import os
import pandas as pd


SAMPLES = (
    pd.read_csv(os.path.abspath(config["samples"]), sep="\t", dtype=str)
    .set_index("sampleid")
    .to_dict(orient="index")
)

# dict : {"chr1": {"vcf": chr1.vcf.gz,"start":1, "end": 99999}, ...}
REFPANEL = (
    pd.read_csv(os.path.abspath(config["genome"]["refpanel"]), sep="\t", dtype=str)
    .set_index("chr")
    .to_dict(orient="index")
)


OUTDIR = "results"

RUN = config["RUN"]


def get_all_results():
    if RUN == "makesample":
        return makesample()
    elif RUN == "sites":
        return run_refpanel()
    elif RUN == "beagle":
        return run_beagle()
    elif RUN == "stitch":
        return run_stitch()
    elif RUN == "phaseless":
        return run_phaseless()
    elif RUN == "all":
        return run_beagle(), run_stitch(), run_phaseless()
    else:
        raise RuntimeError("this is an invalid RUN!")


def get_regions_list_per_chrom(chrom, chunksize=5000000):
    """split chr into chunks given chunksize; return a list of '[start,end]' pairs"""
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


def makesample():
    return expand(
        rules.bamlist.output,
        depth=config["coverage"],
        chrom=config["chroms"],
    )


def run_refpanel():
    return expand(
        rules.get_pos_from_refpanel.output,
        chrom=config["chroms"],
    )

def run_beagle():
    return expand(
        rules.collect_beagle_accuracy.output,
        chrom=config["chroms"],
        depth=config["coverage"],
    )


def run_stitch():
    return expand(
        rules.collect_stitch_accuracy.output,
        chrom=config["chroms"],
        depth=config["coverage"],
    )

def run_phaseless():
    return expand(
        rules.collect_phaseless_v0_accuracy.output,
        chrom=config["chroms"],
        depth=config["coverage"],
    ), expand(
        rules.collect_phaseless_v1_accuracy.output,
        chrom=config["chroms"],
        depth=config["coverage"],
    )

def if_use_af_in_refpanel(wildcards):
    if REFPANEL[wildcards.chrom].get("af"):
        # varify af file. 5 columns: chr pos ref alt af
        return REFPANEL[wildcards.chrom]["af"]
    else:
        return "false"
