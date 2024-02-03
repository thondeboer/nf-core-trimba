#!/usr/bin/env python

import argparse
import logging
import gzip
from typing import List, Optional, Union
from pathlib import Path
import gffutils
from pyfaidx import Fasta
logger = logging.getLogger(__file__)


def get_transcripts(
        db_file: str,
        fasta: str,
        genes_file: List[str] = None,
        ribotish_file: Optional[Union[Path, str]] = None,
        bedgraph_file: Optional[Union[Path, str]] = None,
        p_value: float = 0.05,
        fasta_out: str = "transcripts.fa",
        padding: int = 10,
        min_size: int = 9
) -> None:
    """Get transcript sequences for a the ribotish results, limited by an optional list of genes as a FASTA file with codon_start position in header.

    Args:
        db_file (str, required): Path to the database file. Either this or genome_file and gtf_file are required. Defaults to None.
        genome_file (Optional[Union[Path, str]], required): Path to the genome sequence in FASTA format. Either this and gtf_file or db_file are required. Defaults to None.
        genes_file (Optional[Union[Path, str]], optional): Path to the gene names to include in the output. Defaults to None.
        ribotish_file (str, optional): Path to the ribotish results file. Defaults to None.
        fasta_file (str, optional): Path to the output FASTA file. Defaults to "transcripts.fa.gz".
        padding (int, optional): The amount of basepairs AFTER the start codon to include in the transcript for motif finding. Defaults to 10.
        min_size (int, optional): Minimum bp length required for the transcript, does not include padding. Defaults to 9.

    Returns:
        None

    Note:
        The output FASTA file contains the transcript sequences with codon_start position in the header.
        The codon_start position is 1-based.
        The padding parameter is used to include additional basepairs after the start codon for motif finding.
    """
    genes = set()
    ribotish = set()
    if genes_file is not None and Path(genes_file).exists():
        genes = set([g.lower().strip()
                    for g in open(genes_file, 'r').readlines() if g.strip() != ''])

    if ribotish_file is not None and Path(ribotish_file).exists():
        ribotish = [r.strip().split('\t')
                    for r in open(ribotish_file, 'r').readlines()]
        ribotish_header = ribotish[0]
        ribotish = ribotish[1:]
        # FIlter on p-value if provided
        if p_value is not None:
            ribotish = [r for r in ribotish if (
                (r[12] != 'None') and (float(r[12]) <= p_value))]
        logger.warning(f"Found {len(ribotish)} significant ribotish results.")
        # Ribotish header
        # 0 Gid
        # 1 Tid
        # 2 Symbol
        # 3 GeneType
        # 4 GenomePos
        # 5 StartCodon
        # 6 Start
        # 7 Stop
        # 8 TisType
        # 9 TISGroup
        # 10 TISCounts
        # 11 TISPvalue
        # 12 RiboPvalue
        # 13 RiboPStatus
        # 14 FisherPvalue
        # 15 TISQvalue
        # 16 FrameQvalue
        # 17 FisherQvalue
        # 18 AALen
        # 19 Seq
        rtranscripts = set([r[1] for r in ribotish])
        rgenes = set([r[2].lower() for r in ribotish])
        logger.warning(f"Found {len(rtranscripts)} transcripts in the ribotish file."
                       f" Found {len(rgenes)} genes in the ribotish file.")
    else:
        rtranscripts = set()
        rgenes = set()
    if not genes and not ribotish:
        raise FileNotFoundError(
            "Either genes_file or ribotish must be provided.")
    elif not genes and rgenes:
        genes = rgenes
    elif len(genes) > 0 and len(rgenes) > 0:
        logger.warning(len(genes))
        # Create intersection of genes and ribotish
        genes = genes.intersection(rgenes)
        if not genes:
            logger.warning(
                "No genes found in the ribotish file and the provided gene list, after p-value filtering and intersecting the two.")
            return None
    if db_file is not None and Path(db_file).exists():
        db = gffutils.FeatureDB(str(db_file))
    else:
        raise FileNotFoundError("Database file is required.")
    if fasta is not None and Path(fasta).exists():
        genome = Fasta(fasta)
    else:
        raise FileNotFoundError("Genome file is required.")
    tr_list = []
    for gene in db.features_of_type('gene'):
        g = gene.attributes.get('gene_name', [None])[0]
        logger.debug(f"Processing gene {g}.")
        if g is not None and g.lower() in genes:
            for tr in db.children(gene, featuretype='transcript'):
                if rtranscripts and tr.id not in rtranscripts:
                    continue
                logger.debug(f"Processing transcript {tr.id} for gene {g}.")
                s = ""
                start_codon = 0
                has_utr5 = False
                utr5_size = db.children_bp(
                    tr, child_featuretype='five_prime_utr')
                for f in db.children(tr):
                    if f.featuretype == 'five_prime_utr':
                        has_utr5 = True
                    if f.featuretype == 'exon':
                        s += f.sequence(genome)
                    if f.featuretype == 'start_codon':
                        start_codon = utr5_size
                if has_utr5:
                    tr_list.append((
                        tr.seqid,
                        tr.start,
                        tr.end,
                        tr.strand,
                        g,
                        gene.id,
                        tr.attributes.get('transcript_id', [None])[0],
                        start_codon,
                        s))
    if len(tr_list) > 0:
        # Save as fasta file
        with open(fasta_out, "wt") as f:
            for i, t in enumerate(tr_list):
                if t[7]+padding-1 > min_size:
                    s = t[8][:t[7]+padding-1]
                    seq = '\n'.join(s[i:i+80] for i in range(0, len(s), 80))
                    f.write(f">{t[6]} {t[7]} {t[5]} {t[4]}\n{seq}\n")
        # Filter bedgraph file on transcripts and save
        filtered_transcripts = [t[6] for t in tr_list]
        bedgraph_out = Path(fasta_out).with_suffix(".filtered.bedgraph")
        ribotish_out = Path(fasta_out).with_suffix(".filtered.ribotish.txt")
        with open(bedgraph_file, 'rt') as fin, open(bedgraph_out, "wt") as fout:
            for line in fin:
                if line.split('\t')[0] in filtered_transcripts:
                    fout.write(line)
        with open(ribotish_out, "wt") as fout:
            fout.write('\t'.join(ribotish_header)+'\n')
            for line in ribotish:
                if line[1] in filtered_transcripts:
                    fout.write('\t'.join(line)+'\n')


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


def main(args):
    logger = setup_logging()
    logger.info(f"Starting transcript retrieval for provided genes.")
    get_transcripts(
        db_file=args.db_file,
        fasta=args.fasta,
        genes_file=args.genes_file,
        ribotish_file=args.ribotish_file,
        bedgraph_file=args.bedgraph,
        p_value=args.p_value,
        fasta_out=f"{args.prefix}.fa",
        padding=args.padding,
        min_size=args.min_size
    )
    logger.info(f"Finished transcript retrieval for provided genes.")


if __name__ == "__main__":
    __version__ = "0.0.7"
    parser = argparse.ArgumentParser(
        description="Extract transcripts for selected genes.")

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

    # list of genes to include in the output
    parser.add_argument(
        "--genes_file",
        metavar="FILE",
        type=str,
        help="File containing the gene (symbol) names, one gene per line."
    )

    # Ribotish file with ORF finding results
    parser.add_argument(
        "--ribotish_file",
        metavar="FILE",
        type=str,
        help="Path to file with the Ribotish results."
    )

    # Ribotish file with ORF finding results
    parser.add_argument(
        "--bedgraph",
        metavar="FILE",
        type=str,
        help="Path to file with the transcriptome coverage data."
    )

    # prefix for filenames (usually sample name in nf-core pipeline)
    parser.add_argument(
        "--prefix",
        type=str,
        default="transcripts",
        help="Prefix for the output FASTA file (default: transcripts)"
    )

    # The p-value cutoff for significant ribotish results
    parser.add_argument(
        "--p_value",
        type=float,
        default=0.05,
        metavar="N",
        help="The p-value cutoff for significant Ribotish results (Default: 0.05)."
    )

    # The amount of basepairs AFTER the start codon to include in the transcript for motif finding
    parser.add_argument(
        "--padding",
        type=int,
        default=10,
        metavar="N",
        help="The amount of basepairs AFTER the start codon to include in the transcript for motif finding. (default: 10)."
    )

    #
    parser.add_argument(
        "--min_size",
        type=int,
        default=9,
        metavar="N",
        help="Minimum bp length required for the transcript, does not include padding. Defaults to 9."
    )

    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__}')

    args = parser.parse_args()
    main(args)
