#!/usr/bin/env python3
"""
pileup3seq.py

Convert samtools pileup files to LASER .seq format.

Python 3 compatible replacement for the original LASER pileup2seq.py script.
"""

import argparse
import gzip
import logging
import math
import os
import re
import sys
from collections import OrderedDict
from time import ctime


def banner():
    msg = """
==================================================================
====         pileup3seq: Convert pileup to seq format         ====
====        Python 3 compatible version for PANDA              ====
====        Original LASER pileup2seq by Zhan and Wang         ====
==================================================================
"""
    print(msg)


def open_text(filename):
    """
    Open plain text or gzip-compressed files in text mode.
    """
    with open(filename, "rb") as handle:
        magic = handle.read(2)

    if magic == b"\x1f\x8b":
        return gzip.open(filename, "rt")

    return open(filename, "r")


def merge_regions(regions):
    """
    Merge overlapping BED regions by chromosome.
    """
    for chrom in regions:
        regions[chrom].sort(key=lambda x: x[0])

        merged = []

        for beg, end in regions[chrom]:
            if not merged:
                merged.append([beg, end])
                continue

            last_end = merged[-1][1]

            if beg < last_end:
                if end > last_end:
                    merged[-1][1] = end
            else:
                merged.append([beg, end])

        regions[chrom] = merged

    return regions


class BedFile:
    def __init__(self):
        self.data = {}

    def open(self, filename, trim_chr_prefix=False):
        for line in open_text(filename):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip().split()
            chrom, beg, end = fields[:3]

            if chrom.lower().startswith("chr") and trim_chr_prefix:
                chrom = chrom[3:]

            beg = int(beg)
            end = int(end)

            self.data.setdefault(chrom, []).append([beg, end])

        self.data = merge_regions(self.data)

        n_regions = sum(len(self.data[chrom]) for chrom in self.data)
        return len(self.data), n_regions

    def contain(self, chrom, pos):
        chrom = chrom.replace("chr", "")

        if chrom not in self.data:
            return False

        try:
            pos = int(pos)
        except ValueError:
            return False

        regions = self.data[chrom]

        lo = 0
        hi = len(regions)

        while lo < hi:
            mid = (lo + hi) // 2
            beg, end = regions[mid]

            if pos < beg:
                hi = mid
            elif beg <= pos < end:
                return True
            else:
                lo = mid + 1

        return False


class IndexedFasta:
    """
    Minimal FASTA reader using a .fai index.

    This is only used if -f reference.fa is supplied.
    """

    def __init__(self):
        self.handle = None
        self.index = {}

    def open(self, fasta):
        fai = fasta + ".fai"

        if not os.path.exists(fasta):
            raise FileNotFoundError(f"Reference FASTA does not exist: {fasta}")

        if not os.path.exists(fai):
            raise FileNotFoundError(
                f"Reference FASTA index does not exist: {fai}. "
                f"Create it with: samtools faidx {fasta}"
            )

        self.handle = open(fasta, "r")

        for line in open_text(fai):
            fields = line.strip().split()
            chrom = fields[0].replace("chr", "")
            self.index[chrom] = [int(x) for x in fields[1:]]

        return len(self.index)

    def get_base_1based(self, chrom, pos):
        chrom = chrom.replace("chr", "")
        pos = int(pos) - 1

        if chrom not in self.index:
            return "N"

        chrom_len, file_offset, line_bases, line_width = self.index[chrom]

        if pos < 0 or pos >= chrom_len:
            return "N"

        line_number, line_offset = divmod(pos, line_bases)
        offset = file_offset + line_number * line_width + line_offset

        self.handle.seek(offset)
        return self.handle.read(1)


def check_reference_alleles(reference_file, site_file, logger):
    fasta = IndexedFasta()
    n_chroms = fasta.open(reference_file)

    logger.info(f"Loaded reference file: {reference_file} ({n_chroms} chromosomes)")

    mismatches = 0

    for line in open_text(site_file):
        if not line.strip():
            continue

        fields = line.strip().split()

        if fields[:5] in (
            ["CHR", "POS", "ID", "REF", "ALT"],
            ["CHROM", "POS", "ID", "REF", "ALT"],
        ):
            continue

        chrom, pos, rsid, ref, alt = fields[:5]
        true_ref = fasta.get_base_1based(chrom, pos)

        if true_ref.lower() != ref.lower():
            mismatches += 1
            logger.warning(
                "Mismatched reference allele at %s:%s. "
                "[%s in %s vs. %s in %s]",
                chrom,
                pos,
                ref,
                site_file,
                true_ref,
                reference_file,
            )

    if mismatches:
        logger.error(
            "Detected %d mismatched reference alleles. Please fix and rerun.",
            mismatches,
        )
        sys.exit(1)

    logger.info("No mismatched reference alleles detected.")


def strip_indels(reads):
    """
    Remove insertion/deletion annotations from a samtools pileup read-bases string.

    Pileup indels look like:
      +2AG
      -12ACGTACGTACGT

    The original Python 2 script removed the +/- and length digits, but did not
    fully remove the inserted/deleted sequence. This version removes the full
    pileup indel annotation.
    """
    output = []
    i = 0

    while i < len(reads):
        char = reads[i]

        if char in "+-":
            i += 1

            length_digits = []

            while i < len(reads) and reads[i].isdigit():
                length_digits.append(reads[i])
                i += 1

            if length_digits:
                indel_len = int("".join(length_digits))
                i += indel_len

            continue

        output.append(char)
        i += 1

    return "".join(output)


def count_ref_alt(ref, reads):
    """
    Return ref and alt counts from a samtools pileup read-bases string.
    """
    if len(ref) != 1:
        return 0, 0

    # Remove start-of-read marker and following mapping-quality character.
    reads = re.sub(r"\^.", "", reads)

    # Remove end-of-read markers.
    reads = reads.replace("$", "")

    # Remove insertion/deletion annotations.
    reads = strip_indels(reads)

    # Remove deletion placeholder.
    reads = reads.replace("*", "")

    ref_count = 0
    alt_count = 0

    bases = set("acgtnACGTN")

    for base in reads:
        if base in ".,":  # reference base on forward/reverse strand
            ref_count += 1
        elif base in bases:
            alt_count += 1
        else:
            print(f"Unrecognized base in pileup read string: {reads}", file=sys.stderr)
            break

    return ref_count, alt_count


_PHRED_TO_ERROR = None


def calculate_mean_quality(qual_string):
    """
    Calculate mean base quality in the same style as the original script.
    """
    global _PHRED_TO_ERROR

    if _PHRED_TO_ERROR is None:
        _PHRED_TO_ERROR = {
            chr(i): 10.0 ** (-(i - 33.0) / 10.0)
            for i in range(33, 256)
        }

    total_error = 0.0
    total_count = 0

    for q in qual_string:
        total_error += _PHRED_TO_ERROR.get(q, 1.0)
        total_count += 1

    if total_count == 0:
        return 0

    return round(-10.0 * math.log10(total_error / total_count))


def load_site_file(site_file, logger):
    map_content = [
        line.strip().split()
        for line in open_text(site_file)
        if line.strip()
    ]

    if len(map_content) < 1:
        print("Site file is too short or incorrect.", file=sys.stderr)
        sys.exit(1)

    expected_headers = (
        ["CHR", "POS", "ID", "REF", "ALT"],
        ["CHROM", "POS", "ID", "REF", "ALT"],
    )

    if map_content[0][:5] in expected_headers:
        map_content = map_content[1:]
    else:
        print(
            "Site file does not have a standard header. Continuing anyway.",
            file=sys.stderr,
        )
        print(map_content[0], file=sys.stderr)

    positions = [
        f"{fields[0].replace('chr', '')}:{fields[1]}"
        for fields in map_content
    ]

    rsids = [fields[2] for fields in map_content]

    col_dict = OrderedDict(zip(positions, rsids))

    if len(positions) != len(col_dict):
        duplicated = len(positions) - len(col_dict)
        logger.warning(
            "Total [ %d ] duplicated markers are found in site file [ %s ]",
            duplicated,
            site_file,
        )

        seen = set()
        dup_pos = set()

        for pos in positions:
            if pos in seen:
                dup_pos.add(pos)
            else:
                seen.add(pos)

        if dup_pos:
            logger.warning("Duplicated site: %s", ",".join(sorted(dup_pos)))

    return map_content, col_dict


def load_id_file(id_file, logger):
    pileup_id = {}

    for line in open_text(id_file):
        fields = line.strip().split()

        if not fields:
            continue

        pileup_key = fields[0]

        if "." in pileup_key:
            logger.error(
                "First column in id file contains '.' as in %s. "
                "Only the part before the first dot will be used.",
                pileup_key,
            )
            pileup_key = pileup_key.split(".")[0]

        values = fields[1:]

        if len(values) > 2:
            values = values[:2]
        elif len(values) == 1:
            values = [values[0], values[0]]
        elif len(values) == 2:
            pass
        else:
            logger.error("Wrong line in idFile: %s", line.strip())
            continue

        if pileup_key in pileup_id:
            logger.info("Skipped duplicated id in idFile: %s", pileup_key)
            continue

        pileup_id[pileup_key] = values

    logger.info("%d sample IDs loaded", len(pileup_id))

    return pileup_id


def parse_pileup_file(filename, col_dict, logger):
    result = {}

    for line in open_text(filename):
        fields = line.strip().split()

        if len(fields) == 4:
            # Truncated pileup line.
            chrom, pos, ref, depth = fields
            ref_count, alt_count, qual = 0, 0, 0

        elif len(fields) == 6:
            chrom, pos, ref, depth, reads, quals = fields

            try:
                ref_count, alt_count = count_ref_alt(ref, reads)
                qual = calculate_mean_quality(quals)
            except Exception as exc:
                logger.warning(
                    "Cannot parse pileup data in file %s line: %s",
                    filename,
                    line.strip(),
                )
                raise exc

        else:
            logger.warning(
                "File [ %s ] is empty or not valid. Problem line: [ %s ]",
                filename,
                line.strip(),
            )
            break

        key = f"{chrom.replace('chr', '')}:{pos}"

        if key not in col_dict:
            continue

        result[key] = (ref_count, alt_count, qual)

    return result


def main():
    banner()

    parser = argparse.ArgumentParser(
        description="Convert samtools pileup files to LASER .seq format."
    )

    parser.add_argument(
        "-b",
        dest="bed_file",
        default=None,
        help="BED file of regions to exclude.",
    )

    parser.add_argument(
        "-i",
        dest="id_file",
        default=None,
        help="Optional ID file: pileup_id GWAS_ID.",
    )

    parser.add_argument(
        "-f",
        dest="reference_file",
        default=None,
        help="Optional reference FASTA for checking reference alleles.",
    )

    parser.add_argument(
        "-m",
        dest="site_file",
        required=True,
        help="Site file with columns: CHR POS ID REF ALT.",
    )

    parser.add_argument(
        "-o",
        dest="out_prefix",
        required=True,
        help="Output prefix.",
    )

    parser.add_argument(
        "pileup_files",
        nargs="+",
        help="Input pileup files.",
    )

    args = parser.parse_args()

    logger = logging.getLogger("pileup2seq")
    logger.setLevel(logging.DEBUG)

    file_handler = logging.FileHandler(args.out_prefix + ".log")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(stream_handler)

    logger.info("Started time: %s", ctime())

    bed_file = BedFile()

    if args.bed_file is not None:
        loaded = bed_file.open(args.bed_file, trim_chr_prefix=True)
        logger.info("Loaded BED lines and unique regions: %s", str(loaded))
    else:
        logger.info("Skip loading BED file. All loci will be processed.")

    map_content, col_dict = load_site_file(args.site_file, logger)

    if args.reference_file is not None:
        check_reference_alleles(args.reference_file, args.site_file, logger)
    else:
        logger.info("Skip loading reference file. Reference alleles will not be checked.")

    exclude_pos = [
        bed_file.contain(fields[0].replace("chr", ""), fields[1])
        for fields in map_content
    ]

    positions = [
        f"{fields[0].replace('chr', '')}:{fields[1]}"
        for fields in map_content
    ]

    exclude_pos = {
        positions[idx]
        for idx, should_exclude in enumerate(exclude_pos)
        if should_exclude
    }

    logger.info("Excluding %d on-target markers", len(exclude_pos))
    logger.info("%d markers loaded", len(col_dict))
    logger.info("%d pileup files", len(args.pileup_files))

    pileup_id = None

    if args.id_file is not None:
        pileup_id = load_id_file(args.id_file, logger)

    seq_filename = args.out_prefix + ".seq"

    with open(seq_filename, "w") as seq_file:
        for pileup_file in args.pileup_files:
            result = parse_pileup_file(pileup_file, col_dict, logger)

            if pileup_id is not None:
                key = os.path.basename(pileup_file).split(".")[0]
                gwas_id = pileup_id.get(key, [key, key])

                if key not in pileup_id:
                    logger.info(
                        "Not updating ID: id file does not have correct entry for %s",
                        key,
                    )

                seq_file.write("\t".join(gwas_id))

            else:
                gwas_id = os.path.basename(pileup_file).replace(".pileup", "")
                seq_file.write("\t".join([gwas_id, gwas_id]))

            total_depth = 0
            site_has_pileup = 0
            site_no_pileup = 0

            for key in col_dict:
                if key in exclude_pos:
                    seq_file.write("\t0 0 0")
                    logger.debug("Exclude\t%s", key)
                    continue

                if key in result:
                    ref_count, alt_count, qual = result[key]
                    depth = ref_count + alt_count
                    seq_file.write(f"\t{depth} {ref_count} {qual}")
                    site_has_pileup += 1
                    total_depth += depth
                else:
                    seq_file.write("\t0 0 0")
                    site_no_pileup += 1

            seq_file.write("\n")

            avg_depth = total_depth / (1e-10 + len(col_dict))
            pct_site_pileup = (
                100.0 * site_has_pileup / (1e-10 + site_has_pileup + site_no_pileup)
            )

            logger.info(
                "%s, avgDepth=%.4f, ptgSiteHasPileup=%.4f%%",
                pileup_file,
                avg_depth,
                pct_site_pileup,
            )

    logger.info("Sequence file is generated: %s", seq_filename)

    site_filename = args.out_prefix + ".site"

    with open(site_filename, "w") as site_file:
        site_file.write("\t".join(["CHR", "POS", "ID", "REF", "ALT"]))
        site_file.write("\n")

        for fields in map_content:
            chrom, pos, rsid, ref, alt = fields[:5]
            site_file.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\n")

    logger.info("Site file is generated: %s", site_filename)
    logger.info("Pileup2seq finished at %s", ctime())


if __name__ == "__main__":
    main()