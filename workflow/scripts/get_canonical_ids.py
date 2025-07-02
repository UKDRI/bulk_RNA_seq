"""Script to get the canonical transcript of each gene to get the mapping"""

import argparse
import csv
import logging
import re
from typing import Dict

if snakemake is not None:
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(msg)s",
    )
else:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(msg)s",
    )

logger = logging.getLogger()


def parse_args() -> argparse.Namespace:
    """Parse commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gtf-file",
        type=str,
        required=True,
        help="Path to the Ensembl GTF file to retrieve the ids from",
    )
    parser.add_argument(
        "--mapping-file",
        type=str,
        required=True,
        help="Path to the mapping file to create",
    )
    parser.add_argument("--version", action="version", version="0.1.0")
    return parser.parse_args()


def parse_gtf_file(path: str) -> Dict[str, str]:
    """Parse GTF file to get the gene to canonical transcript mapping

    Args:
        path: Path to the GTF file

    Returns:
        Dictionary where the key is the gene ID and the value is the transcript ID

    Raises:
        ValueError when no mapping was found in the GTF file
    """
    logger.info(path)
    gene_pattern = re.compile(r'gene_id "([^"]+)"')
    transcript_pattern = re.compile(r'transcript_id "([^"]+)"')
    mapping_ids = {}
    with open(path, mode="r", encoding="utf-8", newline="") as gtf_file:
        csv_reader = csv.reader(gtf_file, delimiter="\t")
        for row in csv_reader:
            if len(row) == 9 and row[2] == "transcript":
                attributes = row[8]
                if "Ensembl_canonical" in attributes:
                    gene_match = gene_pattern.search(attributes)
                    transcript_match = transcript_pattern.search(attributes)
                    if gene_match and transcript_match:
                        gene_id = gene_match.group(1)
                        mapping_ids[gene_id] = transcript_match.group(1)
    if len(mapping_ids) == 0:
        logger.error("Did not find canonical transcripts in: '%s'", path)
        raise ValueError(f"Could not find canonical transcripts in {path}")

    return mapping_ids


def write_mapping_id_file(mapping_ids: Dict[str, str], output_file: str):
    """Write the mapping data to a TSV file

    Args:
        mapping_ids: Dictionary with the gene_id as key and the transcript_id as value.
        output_file: Path to the output file
    """
    with open(output_file, mode="w", encoding="utf8", newline="") as mapping_file:
        csv_writer = csv.writer(mapping_file, delimiter="\t")
        csv_writer.writerows(mapping_ids.items())


# Input arguments
if snakemake is not None:
    mapping_data = parse_gtf_file(snakemake.input[0])
    write_mapping_id_file(mapping_data, snakemake.output[0])
elif __name__ == "__main__":
    args = parse_args()
    mapping_data = parse_gtf_file(args.gtf_file)
    write_mapping_id_file(mapping_data, args.mapping_file)
