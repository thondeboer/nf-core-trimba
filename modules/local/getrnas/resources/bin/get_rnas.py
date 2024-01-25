#!/usr/bin/env python

import argparse
import logging
import gzip
from typing import List
from pathlib import Path
import gffutils
from pyfaidx import Fasta
logger = logging.getLogger(__file__)


def get_rnas(
        db_file: str = None,
        fasta: str = None,
        gtf: str = None,
        rnas_to_filter: List[str] = None,
        fasta_out: str = None,
) -> None:
    """Get RNAs of certain type as a fasta file to use in filtering of reads.

    Args:
        db_file (str, optional): Path to the database file. Either this or genome_file and gtf_file are required. Defaults to None.
        fasta (str, optional): Path to the genome sequence in FASTA format. Either this and gtf_file or db_file are required. Defaults to None.
        gtf (str, optional): Path to the GTF file with transcript annotations. Either this and genome_file or db_file are required. Defaults to None.
        rnas_to_filter (List[str], optional): List of RNA types to filter on.
          default: 'scRNA','3prime_overlapping_ncRNA','miRNA','snRNA','macro_lncRNA','sRNA','lincRNA','Mt_rRNA','scaRNA','snoRNA',
                   'rRNA','Mt_tRNA','bidirectional_promoter_lncRNA','misc_RNA'
        fasta_out (str, optional): Path to the output FASTA file. Defaults to "[FAST_IN]_transcripts.fa".

    Returns:
        None
    """
    if rnas_to_filter is None:
        rnas_to_filter = [
            'scRNA', '3prime_overlapping_ncRNA', 'miRNA', 'snRNA', 'macro_lncRNA', 'sRNA', 'lincRNA',
            'Mt_rRNA', 'scaRNA', 'snoRNA', 'rRNA', 'Mt_tRNA', 'bidirectional_promoter_lncRNA', 'misc_RNA'
        ]
    if (db_file is None and fasta is None) or (fasta is None and gtf is None):
        raise ValueError(
            "Either db_file and fasta or fasta and gtf_file must be provided.")
    if db_file is None and fasta is not None and gtf is not None:
        if not Path(fasta).exists():
            raise FileNotFoundError(
                f"Genome file {fasta} does not exist.")
        if not Path(gtf).exists():
            raise FileNotFoundError(f"GTF file {gtf} does not exist.")
        db_file = Path(fasta).with_suffix(".db")
    if db_file is not None and Path(db_file).exists():
        db = gffutils.FeatureDB(str(db_file))
    else:
        db = gffutils.create_db(
            str(gtf),
            dbfn=str(db_file),
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            force=True,
            keep_order=True,
        )
    if fasta_out is None:
        fasta_out = Path(fasta).with_suffix(".RNAtranscripts.fa")
    genome = Fasta(fasta)
    tr_list = []
    chunk_size = 80
    for gene in db.features_of_type('gene'):
        g = gene.attributes.get('gene_name', [None])[0]
        for tr in db.children(gene, featuretype='transcript'):
            attr = tr.attributes['transcript_biotype'][0]
            if attr in rnas_to_filter:
                s = tr.sequence(genome)
                s = '\n'.join(s[i:i+chunk_size]
                              for i in range(0, len(s), chunk_size))
                tr_list.append((tr.seqid, tr.start, tr.end, tr.strand, g, gene.id, attr,
                                tr.attributes.get('transcript_id', [None])[0], s))
    if len(tr_list) > 0:
        # Save as fasta file
        with open(fasta_out, "wt") as f:
            for i, t in enumerate(tr_list):
                f.write(f">{t[7]} {t[4]} {t[5]} RNA_type: {t[6]}\n{t[8]}\n")


def setup_logging() -> logging.Logger:
    """Configure logging for the script.

    Returns:
        logging.Logger: Configured logger instance.
    """
    logging.basicConfig(
        format="%(asctime)-18s-%(levelname)-8s-|%(module)-14s| %(message)s")
    logger = logging.getLogger(__file__)
    logger.setLevel(logging.INFO)
    return logger


def comma_separated_strings(value):
    return [item.strip() for item in value.split(',')]


def main(args):
    logger = setup_logging()
    logger.info(f"Starting transcript retrieval for provided RNAs.")
    get_rnas(
        db_file=args.db_file,
        fasta=args.fasta,
        gtf=args.gtf,
        rnas_to_filter=args.rnas_to_filter,
        fasta_out=None,
    )
    logger.info(f"Finished transcript retrieval for provided RNAs.")


if __name__ == "__main__":
    __version__ = "0.0.7"
    parser = argparse.ArgumentParser(
        description="Extract transcripts for selected RNA types.")

    # Path to the database file - Either this or FASTA plus GTF file is required
    parser.add_argument(
        "--db_file",
        metavar="FILE",
        type=str,
        help=("SQLite database file with transcript annotations.")
    )

    # Genome sequence in FASTA format
    parser.add_argument(
        "--fasta",
        type=str,
        metavar="FILE",
        help="Path to the FASTA file with the genome sequence."
    )

    # GTF file with transcript annotations
    parser.add_argument(
        "--gtf",
        metavar="FILE",
        type=str,
        help="Path to the GTF file with the genome annotations."
    )
    # list of genes to include in the output
    parser.add_argument(
        "--rnas_to_filter",
        metavar="LIST",
        type=comma_separated_strings,
        help="List of RNA types to extract the transcript for (e.g. lincRNA, snoRNA etc.)."
    )

    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__}')

    args = parser.parse_args()
    main(args)
