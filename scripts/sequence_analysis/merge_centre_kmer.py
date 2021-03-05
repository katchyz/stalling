#!/usr/bin/env python
""" Merge the fasta file with the centre point location before use
in the weighted kmer plots. """

import pandas as pd
from Bio import SeqIO
from pathlib import Path


def loadFastaFile(fastaPath: Path) -> pd.DataFrame:
    """ Load the fasta file containing the sequences into a DF. """
    name = []
    seq = []

    # Load the fasta into lists
    with open(fastaPath, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            name.append(record.id)
            seq.append(str(record.seq))

    # Format the list into DF
    fastaSeqs = pd.DataFrame(data=dict(id=name, seq=seq))
    fastaSeqs["len"] = fastaSeqs["seq"].apply(lambda x: len(x))
    return fastaSeqs


def loadCentrePoint(centrePath: Path) -> pd.DataFrame:
    """ Load the centre point array into memory. """
    centrePoints = pd.read_csv(centrePath)
    centrePoints.rename(columns={"tx": "id", "ss_aa": "centre"}, inplace=True)
    return centrePoints


def mergeDFs(fastaSeqs, centrePoints) -> pd.DataFrame:
    """Merge the dataframes.

    There are multiple centre points for some sequences.
    """
    fastaSeqs = fastaSeqs.merge(centrePoints, how="right", validate="one_to_many")
    return fastaSeqs


def retriveRow(row, surround=5):
    """ Retrieve the sequence about the centre point. """
    # Account for the 1-based nature of R
    centre = row["centre"] - 1
    seq = row["seq"]
    len_ = row["len"]

    start = max(centre - surround, 0)
    end = min(centre + surround + 1, len_)
    seqSlice = seq[start:end]

    nullChar = "x"

    # Test for boundary condition at end
    endOverlap = (centre + surround) - len_

    if endOverlap >= 0:
        seqSlice += nullChar * (endOverlap + 1)

    # Test for the boundary conditions at the start
    startOverlap = centre - surround
    if startOverlap <= 0:
        seqSlice = nullChar * (-startOverlap) + seqSlice

    return seqSlice


if __name__ == "__main__":
    data_dir = Path("data")
    save_dir = Path("merged")
    text_files = list(data_dir.glob("*.txt"))

    fasta_files = [data_dir / f"fasta_{i.stem}.fa" for i in text_files]
    fasta_exists = [p.exists() for p in fasta_files]
    if not all(fasta_exists):
        raise ValueError("Not matching fasta files")

    for fasta_path, centre_path in zip(fasta_files, text_files):
        fastaSeqs = loadFastaFile(fasta_path)
        centrePoints = loadCentrePoint(centre_path)
        mergedDF = mergeDFs(fastaSeqs, centrePoints)

        mergedDF["trimmed"] = mergedDF.apply(retriveRow, axis="columns", surround=10)
        save_path = save_dir / f"merged_{fasta_path.stem}.csv"

        mergedDF.to_csv(save_path)
