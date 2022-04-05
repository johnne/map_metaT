import pandas as pd
from collections import defaultdict


def parse_samples(f):
    df = pd.read_csv(f, sep="\t", index_col=0)
    samples = defaultdict(lambda: defaultdict(dict))
    for sample, d in df.iterrows():
        sample_id = f"{sample}_{d['unit']}"
        samples[sample_id]["fq1"] = d["fq1"]
        samples[sample_id]["fq2"] = d["fq2"]
    return samples


def parse_assemblies(f, datadir="data/"):
    df = pd.read_csv(f, sep="\t", index_col=0)
    assemblies = df.to_dict(orient="index")
    return assemblies


def clean_featurecount(sm):
    import os
    dataf = pd.DataFrame()
    for f in sm.input.tsv:
        sample = os.path.basename(f).replace(".fc.tsv", "")
        df = pd.read_csv(f, sep="\t", comment="#", usecols=[0, 1, 6])
        df.columns = ["Geneid", "Chr", sample]
        df.index = df.Chr.map(str) + ["_" + x.split("_")[-1] for x in df.Geneid]
        df.drop(["Geneid","Chr"], axis=1, inplace=True)
        dataf = pd.merge(dataf, df, left_index=True, right_index=True, how="outer")
    dataf.to_csv(sm.output.tsv, sep="\t")


def main(sm):
    toolbox = {"clean_featurecount": clean_featurecount}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)