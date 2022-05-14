from scripts.common import parse_assemblies, parse_samples
wildcard_constraints:
    counts_type="(counts|rpkm)",
    norm_method="(TMM|RLE)"
samples = parse_samples(config["sample_list"])
assemblies = parse_assemblies(config["assembly_list"])