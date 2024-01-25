#!/usr/bin/env python

import argparse
import logging
from typing import List
from pathlib import Path
import gffutils


def gff_createdb(
        gtf: str,
        db_out: str = None,
) -> None:
    """Create a SQLite database from a GTF file.

    Args:
        gtf (str): Path to the GTF file with annotations. Required.
        db_out (str, optional): Path to the output database file. Defaults to "[GTF_IN].db".

    Returns:
        None
    """
    if db_out is None:
        db_out = Path(gtf).with_suffix(".db")
    gffutils.create_db(
        str(gtf),
        dbfn=str(db_out),
        disable_infer_genes=True,
        disable_infer_transcripts=True,
        force=True,
        keep_order=True,
    )

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
    logger.info(f"Starting SQLite database creation from GTF file.")
    gff_createdb(
        gtf=args.gtf,
        db_out=args.db_out,
    )
    logger.info(f"Finished SQLite database creation from GTF file.")


if __name__ == "__main__":
    __version__ = "1.0.0"
    parser = argparse.ArgumentParser(
        description="Create a SQLite database from a GTF file")

    # GTF file with annotations
    parser.add_argument(
        "--gtf",
        metavar="FILE",
        type=str,
        help="Path to the GTF file with the genome annotations."
    )
    # list of genes to include in the output
    parser.add_argument(
        "--db_out",
        metavar="FILE",
        type=str,
        help="Path to the output database file. Defaults to [GTF_IN].db"
    )

    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__}')

    args = parser.parse_args()
    main(args)
