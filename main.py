from __future__ import print_function, unicode_literals, absolute_import, division

# conda
import os
import sys
import pandas as pd
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from glob import glob

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO


# pip
from tqdm import tqdm

# arguments
parser = argparse.ArgumentParser(
    description='''Extraction of sequences of ORF, Bar Code and sgRNA from Nanopore seq of ORFeome-Cas9 fusion library''', formatter_class=RawTextHelpFormatter,
    epilog='''Extract wisely''')
parser.add_argument('-i', default=os.path.join(os.getcwd()), help='directory with images. Default - WD')
parser.add_argument('-o', default=os.path.join(os.getcwd()), help='output dir. Default - WD')
parser.add_argument('--construct_fasta', default=False, type=bool, help='output dir. Default - WD')
parser.add_argument('--csv_input', default=None, type=str, help='input file with sequences, gene ID and clone ID. Default - WD'
                                                                'must have columns: DNASU Clone ID, gene ID, ref_sequence)
parser.add_argument('--fasta_out', default=None, type=str, help='fasta output dir. Default - WD')




# function to construct fasta record
def fasta_constructor(x):
    record = SeqRecord(Seq(x['ref_sequence'], IUPAC.unambiguous_dna),
                       id=x['gene ID'],
                       description=" ".join(["DNASU_CLone_ID:", x['DNASU Clone ID']]))
    return record


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
argsP = parser.parse_args()


if __name__ == "__main__":

    # construct fasta from csv?
    if argsP.construct_fasta == True:

        # import csv
        df = pd.read_csv(argsP.csv_input)

        # keep relevant cols
        df = df[['DNASU Clone ID', 'gene ID', 'ref_sequence']]

        # drop NAs
        df.dropna(inplace=True)

        # apply function to df
        records=df.apply(fasta_constructor, axis=1).tolist()

        # export
        SeqIO.write(records, argsP.fasta_out, "fasta")


    else:
        print("Fasta file is provided by the user. That means you.")



