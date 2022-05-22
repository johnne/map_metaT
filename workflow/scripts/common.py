import pandas as pd
import os
from collections import defaultdict


def parse_samples(f):
    df = pd.read_csv(f, sep="\t", index_col=0)
    samples = defaultdict(lambda: defaultdict(dict))
    for sample, d in df.iterrows():
        sample_id = f"{sample}_{d['unit']}"
        samples[sample_id]["fq1"] = d["fq1"]
        samples[sample_id]["fq2"] = d["fq2"]
        # If no assembly is set in the samples file, assign assembly
        # using sample
        try:
            samples[sample_id]["assembly"] = d["assembly"]
        except KeyError:
            samples[sample_id]["assembly"] = sample
    return samples


def parse_assemblies(f, assembly_dir):
    df = pd.read_csv(f, sep="\t", index_col=0, header=None)
    d = df.to_dict(orient="index")
    assemblies = {}
    for assembly in d.keys():
        if os.path.exists(f"{assembly_dir}/{assembly}/final_contigs.fa"):
            assemblies[assembly] = ""
    return assemblies


def clean_featurecount(sm):
    import os
    dataf = pd.DataFrame()
    for f in sm.input.tsv:
        sample = os.path.basename(f).replace(".fc.tsv", "")
        df = pd.read_csv(f, sep="\t", comment="#", usecols=[0, 1, 5, 6])
        df.columns = ["Geneid", "Chr", "Length", sample]
        df.index = df.Chr.map(str) + ["_" + x.split("_")[-1] for x in df.Geneid]
        df.drop(["Geneid","Chr"], axis=1, inplace=True)
        dataf = pd.merge(dataf, df, left_index=True, right_index=True, how="outer")
        try:
            dataf.drop("Length_y", axis=1, inplace=True)
        except KeyError:
            continue
        dataf.rename(columns={"Length_x": "Length"}, inplace=True)
    dataf.to_csv(sm.output.tsv, sep="\t")


def process_and_sum(q_df, annot_df):
    # Merge annotations and abundance
    # keep ORFs without annotation as "Unclassified"
    annot_q_df = pd.merge(annot_df, q_df, left_index=True, right_index=True,
                          how="right")
    annot_q_df.fillna("Unclassified", inplace=True)
    feature_cols = annot_df.columns
    annot_q_sum = annot_q_df.groupby(list(feature_cols)).sum().reset_index()
    annot_q_sum.set_index(feature_cols[0], inplace=True)
    return annot_q_sum


def sum_to_features(abundance, parsed):
    parsed_df = pd.read_csv(parsed, index_col=0, sep="\t")
    abundance_df = pd.read_csv(abundance, index_col=0, sep="\t")
    abundance_df.drop("Length", axis=1, inplace=True, errors="ignore")
    feature_sum = process_and_sum(abundance_df, parsed_df)
    return feature_sum


def count_features(sm):
    """
    Counts reads mapped to features such as KOs, PFAMs etc.

    :param sm:
    :return:
    """
    feature_sum = sum_to_features(sm.input.abund, sm.input.annot[0])
    feature_sum.to_csv(sm.output[0], sep="\t")


def get_sample_counts(samples, db):
    """
    Sub-function for iterating each sample and looking up the corresponding
    count from its assembly

    :param samples:
    :param db:
    :return:
    """
    d = {}
    annots_dict = {}
    index_col = 0
    if db=="modules":
        index_col=4
    for sample in samples.keys():
        assembly = samples[sample]["assembly"]
        f = f"results/{assembly}/{db}.parsed.counts.tsv"
        df = pd.read_csv(f, sep="\t", index_col=index_col)
        # Extract annotation columns
        annots = df.loc[:, df.dtypes == object]
        annots_dict.update(annots.to_dict(orient="index"))
        d[sample] = df.loc[:, sample]
    counts_df = pd.DataFrame(d)
    counts_df.fillna(0, inplace=True)
    counts_df = pd.merge(pd.DataFrame(annots_dict).T, counts_df, right_index=True, left_index=True)
    return counts_df


def extract_counts(sm):
    """
    This function extracts counts of features for each sample from its
    corresponding single-sample assembly.

    :param sm:
    :return:
    """
    samples = parse_samples(sm.params.sample_list)
    db = sm.wildcards.db
    counts_df = get_sample_counts(samples, db)
    counts_df.index.name = db
    counts_df.to_csv(sm.output[0], sep="\t")


def main(sm):
    toolbox = {"clean_featurecount": clean_featurecount,
               "count_features": count_features,
               "extract_counts": extract_counts}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)