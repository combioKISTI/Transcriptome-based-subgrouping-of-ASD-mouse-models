#!/usr/bin/env python3
"""
01_id_mapping -- Ensembl gene ID to HGNC symbol (human only)

For the cross-species step the manuscript recommends BioMart release 102,
but for fast preprocessing this script uses the BioTools API
(https://biotools.fr/human/ensembl_symbol_converter/) to resolve large lists
of human Ensembl gene IDs to HGNC symbols.

Mouse data are typically already in symbol form, so this step is unnecessary
for the mouse pipeline.

Usage:
    python ensembl_to_symbol.py <input_ids.txt> <output.tsv>

Input  : one ENSG ID per line (any version suffix such as .14_1 is stripped)
Output : input_id<TAB>ensembl_clean<TAB>symbol
"""
import sys
import re
import json
import time
from typing import List, Dict
import requests

API_URL = "https://biotools.fr/human/ensembl_symbol_converter/"


def extract_ENSG(s: str):
    """Return the canonical 'ENSG' + 11-digit ID, or None if not present."""
    if not s:
        return None
    m = re.search(r"ENSG\d{11}", s.strip().upper())
    return m.group(0) if m else None


def chunked(items: List[str], n: int):
    for i in range(0, len(items), n):
        yield items[i:i + n]


def post_convert(ids: List[str], max_retries: int = 3, timeout: int = 30) -> Dict[str, str]:
    """POST a batch of IDs to the BioTools converter and return {ENSG: SYMBOL}."""
    if not ids:
        return {}
    payload = {"api": 1, "ids": json.dumps(ids)}
    for attempt in range(1, max_retries + 1):
        try:
            r = requests.post(API_URL, data=payload, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except Exception:
            if attempt == max_retries:
                raise
            time.sleep(1.5 * attempt)
    return {}


def main():
    if len(sys.argv) < 3:
        sys.stderr.write(__doc__)
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    with open(in_path, encoding="utf-8") as f:
        raw_ids = [line.strip() for line in f if line.strip()]
    if not raw_ids:
        sys.exit("No valid IDs in input file.")

    clean_ids = [extract_ENSG(x) for x in raw_ids]
    uniq_clean = sorted({c for c in clean_ids if c})

    mapping: Dict[str, str] = {}
    for batch in chunked(uniq_clean, 500):
        mapping.update(post_convert(batch))

    with open(out_path, "w", encoding="utf-8") as out:
        out.write("input_id\tensembl_clean\tsymbol\n")
        for raw, clean in zip(raw_ids, clean_ids):
            sym = mapping.get(clean, "NA") if clean else "NA"
            out.write(f"{raw}\t{clean or 'NA'}\t{sym}\n")

    total = len(raw_ids)
    cleaned = sum(1 for x in clean_ids if x)
    mapped = sum(1 for x in clean_ids if x and mapping.get(x))
    sys.stderr.write(
        f"[OK] total={total} cleaned={cleaned} mapped={mapped} -> {out_path}\n"
    )


if __name__ == "__main__":
    main()
