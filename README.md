# Snakemake workflow: `map_metaT`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/johnne/map_metaT/workflows/Tests/badge.svg?branch=main)](https://github.com/johnne/map_metaT/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for mapping metatranscriptomic reads against a metagenomic reference.

## Overview

In short, the workflow maps paired-end input reads against assemblies and performs 
read counting and normalization both per ORF and per annotation features. 

## Usage

Required files are:

1. Assembled fasta files named `final_contigs.fa` and GFF file named
   `final_contigs.gff` under each assembly subdirectory in `assembly_dir` 
   specified in the config file. Each assembly to analyse should be specified in
   a text file given by `assembly_list` in the config file.
2. Paths to paired-end fastq files (preprocessed as necessary) specified in a 
   tab-delimited text file (given by `sample_list` in the config file).
3. An annotation directory for each assembly (given by `annotation_dir`) where
   each assembly has its own subdirectory containing tab-delimited annotation
   files in the form of `<annotation_dir>/<assembly>/<db>.parsed.tsv` where
   `db` can be for instance 'pfam', 'enzymes', 'kos', 'modules', 'pathways' etc.

### Example setup

A typical yaml-format configuration file, by default located at `config/config.yml`
may look like the following. 

```yaml
sample_list: "samples.tsv"
assembly_list: "assemblies.tsv"
assembly_dir: "data/assembly"
annotation_dir: "data/annotation"
featurecounts:
  settings: "-t CDS -g ID -M -Q 10 -B -p"
bowtie2:
  settings: "--very-sensitive"
annotation_dbs:
  - "pfam"
  - "protein-models"
  - "enzymes"
  - "kos"
  - "modules"
  - "pathways"
norm_methods:
  - "CSS"
  - "TMM"
  - "RLE"
```

An example of what the `samples.tsv` file can look like is:

| sample  | unit | fq1                            | fq2                            |
|---------|------|--------------------------------|--------------------------------|
| sample1 | 1    | /proj/data/sample1_R1.fastq.gz | /proj/data/sample1_R2.fastq.gz |
| sample2 | 1    | /proj/data/sample2_R1.fastq.gz | /proj/data/sample2_R2.fastq.gz |

The file `assemblies.tsv` lists what assemblies to analyse and only contains
each assembly name, one per row:

| assembly  | 
|-----------|
| assembly1 |
| assembly2 |

Each assembly, here `assembly1` and `assembly2` must have dedicated subdirectories
under the path given by parameter `assembly_dir`, in this case `data/assembly/`
meaning that there should exist:

```
| data/assembly/assembly1
└── final_contigs.fa
data/assembly/assembly2
└── final_contigs.fa
```

in order for the workflow to work.

Furthermore, in order for the workflow to be able to sum transcripts to the level
of annotated features, such features must be supplied for the metagenomic assembly
under the path given by `annotation_dir`. These annotation files should be in 
tab-delimited format and contain the open reading frame in the first column, 
with subsequent columns dedicated to feature IDs and/or descriptions. As an 
example, for output from running `pfam_scan` on protein sequences against the
PFAM database the output, found in file `pfam.parsed.tsv` may look like this:

| orf           | pfam    | pfam_name                     | clan   | clan_name               |
|---------------|---------|-------------------------------|--------|-------------------------|
| k141_100000_1 | PF01546 | Peptidase family M20/M25/M40  | CL0035 | Peptidase clan MH/MC/MF |
| k141_100003_2 | PF01368 | DHH                           | CL0137 | HAD superfamily |

With this setup, metatranscriptomic reads from each sample will be mapped to each
metagenomic assembly using `bowtie2`, and reads assigned to each orf counted
using `featureCounts` from the `subread` package. These read counts are then 
summed to the level of `pfam` id and normalized.

### Exmple command

To run the workflow with the above setup, use the following command: 

```bash
snakemake --use-conda --configfile config/config.yml -j 100 -rpk --profile slurm
```

The `--profile slurm` flag instructs the workflow to submit jobs to the SLURM
workload manager, using runtime and threads according to the specific rule settings.
The only thing you need to change for this to work is to update the `account:`
setting in the `config/cluster.yaml` file with your account id (_e.g._ snicXXXX-Y-ZZ).

```yaml
__default__:
  account: staff
```

### Output

Using the setup above, expected output for `assembly1` is:

```
| results/assembly1
├── count
   ├── sample1_1.fc.tsv
   ├── sample1_1.fc.tsv.summary
   ├── sample2_1.fc.tsv
   ├── sample2_1.fc.tsv.summary
├── counts.tsv
├── pfam.parsed.counts.tsv
├── pfam.parsed.CSS.tsv
├── pfam.parsed.RLE.tsv
├── pfam.parsed.TMM.tsv
├── report.html
└── rpkm.tsv
```

Under the `count/` subdirectory the raw output from `featureCounts` is located:
`count/sample1_1.fc.tsv` etc.

The file `counts.tsv` contains assigned read counts for each orf in each mapped
sample while `rpkm.tsv` has read counts normalied to reads per kilobase million.
The `report.html` file has summarized results from `multiqc` on mapped and 
assigned read metrics.

All files starting with `pfam.` contain either raw (`pfam.parsed.counts.tsv`) or
normalized read counts for annotated pfams. 

#### Collated output
Collated tables of metatranscriptomic counts and normalized values will be output to `tables/`. However, only samples with matching metagenomic assemblies listed in the `assembly_list` file, **and** with a corresponding `final_contigs.fa` under a subdirectory of `assembly_dir` will have values in these files as they only show abundances for features where a metatranscriptomic sample could be linked directly to its metagenomic assembly.

As an example, in a set up where:

- the config file contains:
```yaml
assembly_dir: "data/assemblies"
assembly_list: "assemblies.txt"
sample_list: "samples.tsv"
```
- a metatranscriptomic sample named `sample1` is specified in `samples.tsv`, and
- there is a metagenomic assembly at `data/assemblies/sample1/final_contigs.fa`

then `sample1` will also have values in collated files under `tables/`. For a metatranscriptomic sample, _e.g._ `sample2` listed in `samples.tsv` but without a corresponding metagenomic assembly, counts and normalized values for features in the other existing assemblies will be found under _e.g_ `results/sample1/kos.parsed.rpkm.tsv` _etc._
