#!/usr/bin/env python

"""Take a set of sequences and plot a weighted heatmap."""

from dataclasses import dataclass
from pathlib import Path
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class KmerHeatmap:
    """

    Parameters
    ----------
    sequences: pd.DataFrame
        Table containing the sequences and related data
    readCol: str
        The name of the column in the sequence table containing the
        sequence string
    nKmers: int (default = 20)
        The number of kmers to include on the plot
    kmerLen: int (default = 6)
        The size of the k-mer to consider
    """

    sequences: pd.DataFrame = None

    kmerLen: int = 6
    nBins: int = 0
    nKmers: int = 20
    readCol: str = "css_seq"
    bgCol: str = None

    def __post_init__(self):
        if self.sequences is not None:
            self.sequenceLength = len(self.sequences.iloc[0][self.readCol])
            self.n_sequences = len(sequences)

        self.possibleKmers = self.sequenceLength - self.kmerLen + 1
        if self.nBins == 0:
            self.nBins = self.possibleKmers
        # if self.sequences is None:
        #   print('DEBUG: Sequences file not provided. ')

    def getKmerDF(self):
        """Create a Dataframe of the kmers with the score of the region of interest.

        We calculate the score of the region as the frequency of the kmers in the
        region of interest against the frequency in the background values.
        """

        # Create a dictionary of the kmer counts
        roiKmerDict = self.countKmers(self.sequences[self.readCol])

        # Merge into a dataframe and rename the columns
        roiKmerDF = pd.DataFrame.from_dict(data=roiKmerDict, orient="index")
        roiKmerDF.reset_index(inplace=True)
        roiKmerDF.rename(columns={"index": "kmers", 0: "roi"}, inplace=True)
        roiKmerDF["roiFreq"] = roiKmerDF["roi"] / roiKmerDF["roi"].sum()

        if self.bgCol is not None:
            # Background level
            bgKmerDict = self.countKmers(self.sequences[self.bgCol])

            # Merge into a dataframe and rename the columns
            bgKmerDF = pd.DataFrame.from_dict(data=bgKmerDict, orient="index")
            bgKmerDF.reset_index(inplace=True)
            bgKmerDF.rename(columns={"index": "kmers", 0: "bg"}, inplace=True)
            bgKmerDF["bgFreq"] = bgKmerDF["bg"] / bgKmerDF["bg"].sum()

            # Merge the dataframes together
            kmerDF = roiKmerDF.merge(bgKmerDF, how="outer", on="kmers")
            kmerDF["score"] = (kmerDF["roiFreq"] - kmerDF["bgFreq"]) * 1e3
        else:
            kmerDF = roiKmerDF
            kmerDF["score"] = kmerDF["roiFreq"] * 1e3

        kmerDF.sort_values("score", inplace=True, ascending=False)

        self.kmerScores = kmerDF
        return kmerDF

    def getHighestScoringKmers(self):
        """Get n highest scoring kmers.

        We can extend this to allow for reversing or symmetrical terms."""
        if not hasattr(self, "kmerScores"):
            self.getKmerDF()

        if self.nKmers > 0:
            highestScoringKmers = self.kmerScores.head(self.nKmers)
        else:
            self.nKmers *= -1
            highestScoringKmers = self.kmerScores.tail(self.nKmers)

        self.highestScoringKmers = highestScoringKmers.reset_index()

        return self.highestScoringKmers

    def matchATtRACTDB(self, exact=False):
        """ Match the highest scoring kmers against the ATtRACT DB"""
        ATtRACTPath = "/home/carl/projects/readRMATs/readRMATs/ATtRACT/ATtRACT_db.txt"
        ATtRACTPath = Path(ATtRACTPath)
        ATtRACTDB = pd.read_csv(ATtRACTPath, sep="\t").query(
            'Organism == "Homo_sapiens"'
        )

        tempMatches = []
        for num, row in self.highestScoringKmers.iterrows():
            kmer = row["kmers"].upper().replace("T", "U")
            # kmer = 'GGAGAAAAA'

            if exact:
                match = ATtRACTDB.query("Motif == @kmer").copy()
            else:
                match = ATtRACTDB.query(
                    f'Motif.str.contains("{kmer}") and Organism == "Homo_sapiens"',
                    engine="python",
                ).copy()

            match["kmer"] = row["kmers"]
            tempMatches.append(match)
        matches = pd.concat(tempMatches, ignore_index=True)

        matches = matches[["Gene_name", "Gene_id", "Motif", "kmer", "Pubmed"]]
        matches.to_csv("matchSites.csv")

    def getKmerLocations(self):
        """Given a list of kmers, return an array of kmer counts.

        This is calculated across the entire set of sequences of interest.
        The rows in the array correspond to different kmers, as given in the kmer
        list, the columns are positions in the sequences.
        """

        kmerLocationHeatmap = np.zeros([self.nKmers, self.nBins], dtype=np.int)

        # Create a lookup dict for the kmers
        # Much more efficient reverse lookup
        kmerLookup = OrderedDict()
        for num, row in self.highestScoringKmers.iterrows():
            kmerLookup[row.kmers] = num

        for seqNum, row in self.sequences.iterrows():
            seq = row[self.readCol]

            for kmerPos in range(self.possibleKmers):
                kmer = seq[kmerPos : kmerPos + self.kmerLen]
                if kmer in kmerLookup:
                    bin_ = self._getBinNumber(kmerPos)
                    kmerNum = kmerLookup[kmer]
                    kmerLocationHeatmap[kmerNum, bin_] += 1

        self.kmerLocationHeatmap = kmerLocationHeatmap
        return kmerLocationHeatmap

    def plotHeatmap(self, saveName, zoom=None):
        """ """
        if not hasattr(self, "kmerLocationHeatmap"):
            self.getKmerLocations()

        fig, ax = plt.subplots(figsize=(8.5, 6))

        # Pcolormesh and the +1 axis definitions are required for centring the
        # tick labels
        x = np.arange(self.nBins + 1)
        y = np.arange(self.nKmers + 1)

        im = ax.pcolormesh(
            x - 0.5,
            y - 0.5,
            # self.kmerLocationHeatmap / self.kmerLocationHeatmap.sum(),
            self.kmerLocationHeatmap / self.n_sequences,
            cmap="inferno",
        )

        # label the ticks with the highest scoring kmer names
        fontSize = 10 if self.nKmers < 22 else 6

        scored_kmers = self.highestScoringKmers["kmers"].values
        y_ticks = y[: len(scored_kmers)]

        ax.set_yticks(y_ticks)
        ax.set_yticklabels(scored_kmers, family="monospace", size=fontSize)
        # Flip the axis, remove the blank space
        ax.set_ylim(self.nKmers - 0.5, -0.5)

        # Give coordinates relative to the stall
        xticks = np.linspace(0, self.nBins - 1, 5)
        # xticks = np.linspace(0, self.possibleKmers - 1, 5)
        ax.set_xticks(xticks)
        if self.nBins == 200:
            ax.set_xticklabels(["-100", "-50", "Stall", "+50", "+100"])
        else:
            endPoint = self.nBins // 2
            midPoint = self.nBins // 4
            ax.set_xticklabels(
                [
                    f"-{endPoint}",
                    f"-{midPoint}",
                    "Stall",
                    f"+{midPoint}",
                    f"+{endPoint}",
                ]
            )

        if zoom is not None:
            ax.set_xlim(self.nBins // 2 - zoom, self.nBins // 2 + zoom)
        else:
            ax.set_xlim(-0.5, self.nBins - 0.5)

        ax.set_ylabel(f"Higest Scoring {self.kmerLen}-mers")
        ax.set_xlabel(f"Position in sequence (relative to stall)")

        plt.colorbar(im, label=f"{self.kmerLen}-mer frequency in region")
        plt.tight_layout()
        plt.savefig(saveName, dpi=300)

    def countKmers(self, sequences):
        """ Count the occurrences of each kmer in a list of sequences. """
        kmerDict = {}
        for seq in sequences:
            for i in range(self.possibleKmers):
                kmer = seq[i : i + self.kmerLen]
                if kmer not in kmerDict:
                    kmerDict[kmer] = 0
                kmerDict[kmer] += 1

        return kmerDict

    def _getBinNumber(self, kmerPos):
        """Convert the position in the array into a bin.

        For now, we assume that the sequence is 200 bases long.
        """
        if kmerPos > self.possibleKmers:
            raise ValueError("Looking for bin beyond end of sequence. ")

        return int(self.nBins * kmerPos / self.possibleKmers)


class KmerHeatmapGeneric(KmerHeatmap):
    """ Plot the kmer heatmaps, but with an arbitrary sequence length. """

    def countKmers(self, sequences):
        """ Count the occurrences of each kmer in a list of sequences. """
        kmerDict = {}
        for seq in sequences:
            for i in range(self.possibleKmers):
                kmer = seq[i : i + self.kmerLen]
                if "x" in kmer:
                    continue
                if kmer not in kmerDict:
                    kmerDict[kmer] = 0
                kmerDict[kmer] += 1

        return kmerDict


class kmerHeatmapWeighted(KmerHeatmapGeneric):
    """Score the k-mers by their location relative to the stall site.

    We do this by changing the kmer counting to increase the score of the kmers
    near the background. We might need to disable the background score as well
    """

    def __init__(self, *args, width=20, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width

    def countKmers(self, sequences):
        kmerDict = {}
        for seq in sequences:
            for i in range(self.possibleKmers):
                kmer = seq[i : i + self.kmerLen]
                # pass for null entry
                if "x" in kmer:
                    continue
                if kmer not in kmerDict:
                    kmerDict[kmer] = 0
                kmerDict[kmer] += self.weightingFunc(i)

        return kmerDict

    def weightingFunc(self, i):
        """Given a position in the sequence, return a score based on the distance
        to the stall site.
        """
        stallLocation = self.possibleKmers / 2
        distanceToStall = abs(stallLocation - i)

        # weighting = np.exp(-distanceToStall/100)
        weighting = np.exp(-distanceToStall / self.width / self.possibleKmers)
        return weighting


if __name__ == "__main__":
    widths = [1000, 1000, 1000]

    merged_data_dir = Path("merged")
    merged_paths = list(merged_data_dir.glob("merged*.csv"))
    save_dir = Path("figs")

    for merged_path in merged_paths:
        short_name = "_".join(merged_path.stem.split("_")[2:])
        for kmerLen in [1, 2, 3]:
            sequences = pd.read_csv(merged_path)
            kmers = kmerHeatmapWeighted(
                width=widths[kmerLen - 1],
                sequences=sequences,
                nKmers=21,
                kmerLen=kmerLen,
                readCol="trimmed",
                nBins=11,
            )

            save_path = save_dir / f"kmer-{kmerLen}-{short_name}.png"
            kmers.getHighestScoringKmers()
            kmers.getKmerLocations()
            kmers.plotHeatmap(saveName=save_path)
