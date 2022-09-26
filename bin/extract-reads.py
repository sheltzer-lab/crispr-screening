#!/usr/bin/env python

import re
import gzip
from typing import Iterator, List, Mapping, Tuple
from pandas import DataFrame, read_csv
import argparse


def main():
    # Parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("primer")
    parser.add_argument("reads")
    parser.add_argument("samples")
    parser.add_argument("library")
    args = parser.parse_args()

    PRIMER = args.primer
    READS = args.reads
    SAMPLES = args.samples
    LIBRARY = args.library

    # Read in samples
    samples = read_csv(SAMPLES)

    # Extract adapters
    adapters = samples["sequence"].tolist()

    # Compile the regex patterns
    patterns = [x for x in generate_regexes(adapters, PRIMER)]

    # Interate over reads and find regex matches and maintain counts
    counts: Mapping[str, Mapping[str, int]] = {}
    with gzip.open(READS, "rt") as reads:
        matches = extract_matches(reads, patterns)
        counts = count_matches(matches)

    # Convert dictionary of dictionaries to matrix of values
    df = DataFrame.from_dict(counts, orient="index")

    # Row names to column
    df.index.name = "adapter"
    df.reset_index(inplace=True)

    # Convert to flat format: adapter, sequence, count
    df = df.melt(id_vars="adapter", var_name="sequence", value_name="count")
    df["count"] = (
        df["count"].fillna(1).astype(int)
    )  # Mageck implodes when value is 0/NaN

    # Merge sgRNA library,
    df = df.merge(
        read_csv(
            LIBRARY,
            delimiter="\t",
            header=None,
            dtype=str,
            names=["gene", "id", "sequence"],
        ),
        on="sequence",
        how="inner",
    )

    # Prep for writing Mageck files
    for sample in samples["sample"].unique().tolist():
        initialAdapter = samples.loc[
            (samples["sample"] == sample) & (samples["time"] == "initial"), ["sequence"]
        ].values[0][0]
        finalAdapter = samples.loc[
            (samples["sample"] == sample) & (samples["time"] == "final"), ["sequence"]
        ].values[0][0]
        initialDf = df.loc[(df["adapter"] == initialAdapter), ["id", "gene", "count"]]
        finalDf = df.loc[df["adapter"] == finalAdapter, ["id", "gene", "count"]]
        finalDf = initialDf.merge(
            finalDf, on=["id", "gene"], suffixes=["_initial", "_final"]
        ).rename(
            columns={
                "id": "sgRNA",
                "gene": "Gene",
                "count_initial": sample + ".initial",
                "count_final": sample + ".final",
            }
        )
        finalDf.to_csv("count-" + sample + "-i.f.csv", index=False)


def generate_regexes(adapters: List[str], primer: str) -> Iterator[re.Pattern]:
    for adapter in adapters:
        yield re.compile("(" + adapter + ")" + primer + "(.{21})")


def extract_matches(
    reads: gzip.GzipFile, patterns: List[re.Pattern]
) -> Iterator[Tuple[str, str]]:
    while line := reads.readline().rstrip():
        for pattern in patterns:
            if found := pattern.search(line):
                adapter = found.group(1)
                sequence = ""

                # this 'G' character is not in the guide sequence,
                # it is likely in between sequencing primer and
                # sgRNA sequence in the vector, this can be deleted
                if found.group(2).startswith("G"):
                    sequence = found.group(2)[1:]
                else:
                    sequence = found.group(2)[:-1]

                yield (adapter, sequence)


def count_matches(
    matches: Iterator[Tuple[str, str]]
) -> Mapping[str, Mapping[str, int]]:
    counts: Mapping[str, Mapping[str, int]] = {}
    for adapter, sequence in matches:
        if adapter in counts:
            if sequence in counts[adapter]:
                counts[adapter][sequence] = counts[adapter][sequence] + 1
            else:
                counts[adapter][sequence] = 1
        else:
            counts[adapter] = {}
            counts[adapter][sequence] = 1

    return counts


if __name__ == "__main__":
    main()
