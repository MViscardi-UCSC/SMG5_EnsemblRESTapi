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

import regex as re
import pandas as pd

# Pandas default would cut off any columns beyond 5 so:
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)


def parseFastaToDataframe(fasta_filepath: str) -> pd.DataFrame:
    with open(fasta_filepath, 'r') as f:
        fasta_string = f.read()
        hits = re.findall('(>.+\n)(([ACTGN]+\n)+)', fasta_string, re.MULTILINE)
        output = [hits[hit_num][0:2] for hit_num in range(len(hits))]
        
        clean_output = []
        for hit in output:
            hit_id, hit_seq = hit
            hit_id = hit_id.replace('\n', '')
            hit_seq = hit_seq.replace('\n', '')
            clean_output.append([hit_id, hit_seq])
        fasta_df = pd.DataFrame(clean_output)
        # This part is specific to how ensemblJSONtoFASTA.py is prepping the fasta files
        #   It will be a source of issue if used with other fasta files
        fasta_df[['g', 't']] = fasta_df[0].str.extract(r'>(.+)\((.+)\)')
        fasta_df = fasta_df[['g', 't', 1]]
        fasta_df = fasta_df.rename(columns={'g': 'gene_id',
                                            't': 'transcript_id',
                                            1: 'sequence',
                                            })
        # print(fasta_df)
        print(fasta_df.nunique())  # There are often several transcripts per gene, this helps to show that
        #                            The presence of non-unique 'sequences' makes me think that there may be redundant
        #                            transcript IDs pointing to same sequence in Ensembl
        return fasta_df


def mergeDataframesOnID(cdna_df: pd.DataFrame, cds_df: pd.DataFrame) -> pd.DataFrame:
    # With an inner merge (accepting only overlap) it doesn't matter which of these start the merge
    merged_df = cdna_df.merge(cds_df,
                              on=['gene_id', 'transcript_id'],
                              suffixes=("_cDNA", "_CDS"))
    print(merged_df.nunique())
    return merged_df


def matchCDStocDNA(merged_df: pd.DataFrame) -> pd.DataFrame:
    # Below doesn't want to work :/
    # merged_df["3UTR"] = merged_df['sequence_cDNA'].str.extract(rf".*{merge_df['sequence_CDS']}(.+)")
    for index, row in merge_df.iterrows():
        match = re.findall(rf".*{row['sequence_CDS']}(.+)", row['sequence_cDNA'])
        if match:
            merge_df.loc[index, '3UTR'] = match[0]
    # matched_df = merge_df
    matched_df = merge_df[merge_df['3UTR'].notnull()]
    print(matched_df.nunique())  # The 10 or so more unique UTRs compared to genes probably indicates alt-spliced UTRs?
    return matched_df


def findSTOPandTSCs(matched_df: pd.DataFrame) -> pd.DataFrame:
    matched_df['STOP'] = matched_df['sequence_CDS'].str.slice(start=-3)
    print(matched_df['STOP'])
    print(matched_df.nunique())
    return matched_df


if __name__ == '__main__':
    
    cDNA = "200805_Ensembl_PTHR15696_SF1_cDNA.fasta"
    print("\ncDNA:")
    cDNA_df = parseFastaToDataframe(cDNA)
    
    CDS = "200805_Ensembl_PTHR15696_SF1_CDS.fasta"
    print("\nCDS:")
    CDS_df = parseFastaToDataframe(CDS)
    
    print("\nMerged:")
    merge_df = mergeDataframesOnID(cDNA_df, CDS_df)
    
    print("\nSequences Matched:")
    match_df = matchCDStocDNA(merge_df)
    
    print("\nFound Stops and following codons")
    findSTOPandTSCs(match_df)
