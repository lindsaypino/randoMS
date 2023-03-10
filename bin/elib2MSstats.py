import sys
import argparse
import pandas as pd

sys.stderr.write('Imported required packages.\n')



class SQLreader:
    """Abstracts reading specific tables and columns from .elib files.
    """

    def get_peptidescores(self, elibpath):
        entry_table = pd.read_sql_table("peptidescores", f"sqlite:///{elibpath}")  # "peptidescores" database
        selected_column = ["PeptideModSeq", "SourceFile", "IsDecoy"]
        pep_table = entry_table[selected_column]
        pep_table = pep_table[pep_table['IsDecoy'] == False]
        pep_table = pep_table[['PeptideModSeq', 'SourceFile']]

        return pep_table
    
    def get_peptidetoprotein(self, elibpath):
        entry_table = pd.read_sql_table("peptidetoprotein", f"sqlite:///{elibpath}")  # "peptidescores" database
        selected_column = ["PeptideSeq", "ProteinAccession", "isDecoy"]
        selected_table = entry_table[selected_column]
        selected_table = selected_table[selected_table['isDecoy'] == False]
        selected_table = selected_table[['PeptideSeq', 'ProteinAccession']]

        return selected_table
    
    def get_peptidequants(self, elibpath):
        entry_table = pd.read_sql_table("peptidequants", f"sqlite:///{elibpath}")  # "peptidequants" database
        selected_column = ["PeptideSeq", "TotalIntensity", "SourceFile", "PrecursorCharge"]
        selected_table = entry_table[selected_column]
        
        return selected_table

    
sqlreader = SQLreader()
sys.stderr.write('Imported required functions.\n')


# usage statement and input descriptions
parser = argparse.ArgumentParser(
    description="this is a script that converts an Encyclopedia *.elib file into an MSstats-acceptable table. \
                You need to make a metadata annotation file associating your filenames to the Condition and \
                BioReplicate, but other than that it should be pretty automatic.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('elib_file', type=str,
                    help="an *.elib file from Encyclopedia's SAVE QUANT REPORT function.")


##
## MAIN
##

# parse arguments from command line
args = parser.parse_args()
elib_file = args.elib_file

## Pull peptide/precursor intensities
pepscores_df = sqlreader.get_peptidequants(elib_file)

## Pull peptide-to-protein assignments
pepprot_df = sqlreader.get_peptidetoprotein(elib_file)

## make a little template annotation file with [Run, BioReplicate, Condition]
template_df = pepscores_df[['SourceFile']].drop_duplicates()
template_df['BioReplicate'] = ''
template_df['Condition'] = ''
template_df.to_csv('./elib2msstats_template.csv', index=False)
sys.stdout.write("OK, friendo, check your working directory for a file called elib2msstats_template.csv\n \
    Fill it out, save it, and close it -- but don't change the name or location!\n")
input("Press Enter when you're finished to continue...")

## read in the completed template and merge the dataframes together
template_df = pd.read_csv('./elib2msstats_template.csv')
msstats_df = pd.merge(pd.merge(pepprot_df, pepscores_df, how='outer', on='PeptideSeq'), 
                        template_df, how='outer', on='SourceFile')

## Add a few required columns and rename header to match MSstats convention
msstats_df['IsotopeLabelType'] = 'L'
msstats_df['FragmentIon'] = 'y0'
msstats_df['ProductCharge'] = '1'
msstats_df.rename(columns={'PeptideSeq': 'PeptideSequence',
                          'ProteinAccession': 'ProteinName',
                          'SourceFile': 'Run',
                          'TotalIntensity': 'Intensity'}, inplace=True)
#print(msstats_df.head())
msstats_df.to_csv("./elib2msstats_output.csv", index=False)
sys.stdout.write("Done! Check this directory for a file called elib2msstats_output.csv")
