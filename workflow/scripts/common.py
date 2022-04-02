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