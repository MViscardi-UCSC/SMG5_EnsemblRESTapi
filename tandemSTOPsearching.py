"""
tandemSTOPsearching.py
Marcus Viscardi      Aug 6, 2020

Going to use the CDSs and cDNA sequences pulled from Ensembl for the SMG5 homologue family to try and assess the
    prevalence of tandem stop codons (TSCs) in SMG5

pseudo-code:
    1. parse fasta files into pandas dataframes
    2. merge those dataframes by matching gene IDs (or maybe transcript IDs?)
    3. find the CDS sequence within the cDNA sequence and identify stop codon
    4. using identified index, look down-sequence from STOP for TSCs
    5. plot prevalence with frequency on the y and number of amino acids until next stop on x
        (see https://www.nature.com/articles/nature18308/figures/5)
    6. profit?
"""

import pandas as pd


def parseFastaToDataframe(fasta_filepath: str) -> pd.DataFrame:
    pass


def mergeDataframesOnID(cdna_df: pd.DataFrame, cds_df: pd.DataFrame) -> pd.DataFrame:
    pass


def matchCDStocDNA(merged_df: pd.DataFrame) -> pd.DataFrame:
    pass


def findSTOPandTSCs(merged_df: pd.DataFrame) -> pd.DataFrame:
    pass
