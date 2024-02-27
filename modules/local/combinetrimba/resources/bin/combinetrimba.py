#!/usr/bin/env python

import pandas as pd
import re
import argparse
from typing import Union, List, Tuple, Optional
import logging

# Set up logging
logger = logging.getLogger(__name__)


def read_fasta(file_path: str) -> pd.DataFrame:
    """
    Reads a FASTA file and extracts information into a pandas DataFrame.

    Parameters:
        file_path (str): The path to the FASTA file.

    Returns:
        pd.DataFrame: A DataFrame containing the extracted data with columns:
                      'Tid', 'StartCodon_pos', and 'fiveprimeutr'.
    """
    data = []
    current_id = None
    current_geneid = None
    current_genesymbol = None
    current_start_codon = None
    current_sequence = []

    # Function to add record to data
    def add_record(geneid, id, genesymbol, start_codon, sequence):
        if id is not None:
            # Join the list of sequence lines into a single string
            sequence_str = ''.join(sequence)

            # Ensure start_codon is an integer and within the sequence length
            start_pos = int(start_codon)
            # YASAR DID NOT WANT THE SPACES AROUND THE START CODON
            # if start_pos < len(sequence_str):
            #     # Insert spaces around the start codon
            #     modified_sequence = sequence_str[:start_pos]  + \
            #         sequence_str[start_pos:start_pos+3] + \
            #         sequence_str[start_pos+3:]
            # else:
            #     # If start_codon is outside the sequence, use the original sequence
            #     modified_sequence = sequence_str

            data.append([
                geneid,
                id,
                genesymbol,
                start_codon,
                sequence_str
            ])

    # Read the file
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # Add previous record to data
                add_record(
                    current_geneid,
                    current_id,
                    current_genesymbol,
                    current_start_codon,
                    ''.join(current_sequence)
                )

                # Reset current record
                parts = line[1:].split()
                current_id = parts[0]
                current_start_codon = parts[1] if len(parts) > 1 else None
                current_geneid = parts[2] if len(parts) > 2 else None
                current_genesymbol = parts[3] if len(parts) > 3 else None
                current_sequence = []
            else:
                # Append sequence
                current_sequence.append(line)

        # Add the last record to data
        add_record(
            current_geneid,
            current_id,
            current_genesymbol,
            current_start_codon,
            current_sequence
        )

    # Create DataFrame
    return pd.DataFrame(data, columns=['Gid', 'Tid', 'Symbol', 'StartCodon_pos', 'fiveprimeutr'])


def extract_meme_motifs(file_path: str) -> Union[pd.DataFrame, str]:
    """
    Extracts all motif data from a specified text file and returns it as a pandas DataFrame.

    The function reads the file, locates each motif section starting and ending with lines of dashes ('-' * 80),
    and extracts the motif name, data, adding these to a DataFrame with columns: 'Motif', 'Name', 'Sequence name', 
    'Start', 'P-value', 'Site'.

    Parameters:
        file_path (str): The path to the text file containing the motif data.

    Returns:
        pd.DataFrame: A DataFrame containing the extracted motif data.
        str: An error message if no motif sections are found.
    """
    # Define the dash pattern for easier readability and future adjustments
    dash_pattern = '-' * 80
    end_pattern = '\*' * 80

    # Read the file
    with open(file_path, 'r') as file:
        content = file.read()

    # Define a regular expression pattern to match the motif sections
    pattern = rf'{dash_pattern}\n\tMotif (.*?) MEME-(\d+) sites sorted by position p-value\n{dash_pattern}\n(.*?)\n{end_pattern}\n'

    # Find all motif sections
    matches = re.findall(pattern, content, re.DOTALL)
    if not matches:
        logger.warning("No motif sections found.")
        return None

    # Extract columns from each motif section
    data = []
    for full_motif, name, motif_data in matches:
        # Extracting just the motif name before the first space
        motif_name = full_motif.split(' ')[0]

        # Split the data into lines and skip the header lines
        lines = motif_data.strip().split('\n')[2:]  # Skipping the header lines

        for line in lines:
            if line.strip():  # skip empty lines
                if line.startswith('---'):
                    break
                # Splitting only into four parts
                parts = line.split(maxsplit=3)
                data.append([motif_name, f'MEME-{name}'] + parts)

    # Create a DataFrame
    return pd.DataFrame(data, columns=[
        'Motif', 'Motif_Name', 'Tid', 'Motif_Start', 'Motif_P-value', 'Motif_Seq'])


def combine_trimba_results(
        transcript_file: str,
        ribotish_file: str,
        meme_file: str,
        sea_file: str,
        prefix: str) -> pd.DataFrame:
    """
    Combines TRIMBA results from transcript, ribotish, and MEME motif files into a single DataFrame.

    Parameters:
        transcript_file (str): Path to the transcript file (FASTA format).
        ribotish_file (str): Path to the Ribotish result file.
        meme_file (str): Path to the MEME motif file.
        sea_file (str): Path to the SEA motif file.
        prefix (str): Prefix for the output results file.

    Returns:
        pd.DataFrame: A DataFrame containing the combined data from all three files.
    """
    t = read_fasta(transcript_file)
    # Ribotsih could be empty, if we have skipped riboSEQ analysis
    try:
        r = pd.read_csv(ribotish_file, sep='\t')
        # Only need Tid and StartCodon_pos columns from transcript file, if we have ribotish data
        rt = r.merge(t[['Tid', 'StartCodon_pos', 'fiveprimeutr']],
                     on='Tid', how='left')
    except:
        rt = t
    try:
        s = pd.read_csv(sea_file, sep='\t')
        if not s.empty and len(s) > 0:
            s.rename(columns={'seq_ID': 'Tid',
                              'Motif_Start': 'SEA_Start'}, inplace=True)
            rts = rt.merge(s, on='Tid', how='left')
        else:
            logger.warning("No SEA motif data found.")
            rts = rt
    except:
        logger.warning("No SEA motif data found.")
        rts = rt
    try:
        m = extract_meme_motifs(meme_file)
        # meme uses 1-notation for start codon, so we need to convert it to 0-notation
        if m is not None:
            m['Motif_Start'] = m['Motif_Start'].astype(int) - 1
            df = rts.merge(m, on='Tid', how='left')
        else:
            logger.warning("No MEME motif data found.")
            df = rts
    except:
        logger.warning("No MEME motif data found.")
        df = rts
    df.to_csv(f'{prefix}.tsv', index=False, header=True, sep='\t')


def main() -> None:
    """
    Main function to parse command line arguments and execute the program.
    """
    __version__ = "0.0.7"
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Combine TRIMBA results from different file sources.")

    # Add arguments
    parser.add_argument(
        "--transcript_file",
        metavar="FILE",
        type=str,
        help="Path to the transcript file (FASTA format)."
    )
    parser.add_argument(
        "--ribotish_file",
        metavar='FILE',
        type=str,
        help="Path to the Ribotish result file."
    )
    parser.add_argument(
        "--meme_file",
        metavar='FILE',
        type=str,
        help="Path to the MEME motif file."
    )
    parser.add_argument(
        "--sea_file",
        metavar='FILE',
        type=str,
        help="Path to the SEA motif file."
    )
    # prefix for filenames (usually sample name in nf-core pipeline)
    parser.add_argument(
        "--prefix",
        type=str,
        default="trimba_results",
        help="Prefix for the output results file (default: trimba_results)."
    )
    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__}')

    # Parse arguments
    args = parser.parse_args()

    # Call the combine_trimba_results function with parsed arguments
    combined_df = combine_trimba_results(
        transcript_file=args.transcript_file,
        ribotish_file=args.ribotish_file,
        meme_file=args.meme_file,
        sea_file=args.sea_file,
        prefix=args.prefix
    )
    # Print or process the combined DataFrame
    print(combined_df)


if __name__ == "__main__":
    main()
