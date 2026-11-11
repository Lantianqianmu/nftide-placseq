#!/usr/bin/env python3
"""Validate/normalize BED peaks for FitHiChIP without third-party packages."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile


def prepare_peaks(input_file, chrom_sizes, output_file, report_file):
    lengths = {}
    with open(chrom_sizes, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 2 or int(fields[1]) <= 0 or fields[0] in lengths:
                raise ValueError("Invalid/duplicate chromosome size on line %d" % number)
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise ValueError("The chromosome-size file is empty")

    source = Path(input_file)
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
        handle.seek(0)
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    opener = gzip.open if compressed else open
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    stats = dict(input_intervals=0, retained_intervals=0, clipped_intervals=0,
                 skipped_nonreference=0, skipped_outside_reference=0)
    unknown = set()
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent,
                                         prefix=".peaks-", delete=False) as target:
            temporary = target.name
            with opener(source, "rt", encoding="utf-8-sig") as handle:
                for number, line in enumerate(handle, 1):
                    line = line.strip()
                    if not line or line.startswith("#") or line.split()[0] in ("track", "browser"):
                        continue
                    fields = line.split("\t") if "\t" in line else line.split()
                    if len(fields) < 3:
                        raise ValueError("Peak line %d has fewer than three BED columns" % number)
                    try:
                        start, end = int(fields[1]), int(fields[2])
                    except ValueError:
                        raise ValueError("Peak line %d requires integer BED start/end coordinates" % number)
                    if end <= start:
                        raise ValueError("Peak line %d has end <= start" % number)
                    stats["input_intervals"] += 1
                    if fields[0] not in lengths:
                        stats["skipped_nonreference"] += 1
                        unknown.add(fields[0])
                        continue
                    clipped_start = max(0, start)
                    clipped_end = min(end, lengths[fields[0]])
                    if clipped_end <= clipped_start:
                        stats["skipped_outside_reference"] += 1
                        continue
                    if (start, end) != (clipped_start, clipped_end):
                        stats["clipped_intervals"] += 1
                        # Preserve the summit's genomic position when clipping narrowPeak.
                        if len(fields) == 10:
                            try:
                                summit = int(fields[9])
                            except ValueError:
                                raise ValueError("Invalid narrowPeak summit on line %d" % number)
                            if summit >= 0:
                                summit += start - clipped_start
                                fields[9] = str(summit if 0 <= summit < clipped_end - clipped_start else -1)
                    fields[1:3] = [str(clipped_start), str(clipped_end)]
                    target.write("\t".join(fields) + "\n")
                    stats["retained_intervals"] += 1
        if not stats["retained_intervals"]:
            raise ValueError("No usable peak intervals overlap the chromosome-size file: %s" % source)
        os.replace(temporary, output)
        temporary = None
    finally:
        if temporary is not None:
            os.unlink(temporary)
    stats.update(input_sha256=digest.hexdigest(),
                 output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                 compressed_input=compressed, skipped_chromosomes=sorted(unknown))
    Path(report_file).write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    try:
        stats = prepare_peaks(args.input, args.chrom_sizes, args.output, args.report)
    except (OSError, ValueError, EOFError) as error:
        parser.exit(1, "Peak preparation failed: %s\n" % error)
    print("Retained %(retained_intervals)d/%(input_intervals)d peak intervals" % stats)


if __name__ == "__main__":
    main()
