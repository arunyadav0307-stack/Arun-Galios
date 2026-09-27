#!/usr/bin/env python3
"""
P0 scaffolding utility.

Extracts the three primary-source PDFs from `Galios Hull.zip` to page-marked
plain text under `sources/`, so that provenance lookups (SOURCE_INVENTORY.md,
discrepancy protocol §7.3, flags F2/F3) can cite page numbers durably.

The prior extraction lived in /tmp/hull_txt/ (outside the repository, not
committed, not persisted). This script regenerates it from the ZIP, which is
the ONLY primary source.

NOT a research computation: no mathematics is evaluated here. This is
document extraction for provenance only.

CAVEAT (recorded in sources/README.md): PDF text extraction mangles much of
the mathematics (lost superscripts, garbled symbols). The output is fit for
*locating* statements by page, NOT for reading formulas verbatim. Any formula
used in the paper must be re-read from the rendered PDF.

Usage:  python3 scripts/extract_sources.py
"""

from __future__ import annotations

import hashlib
import re
import sys
import zipfile
from pathlib import Path

import pymupdf  # PyMuPDF

REPO = Path(__file__).resolve().parent.parent
ZIP = REPO / "Galios Hull.zip"
OUTDIR = REPO / "sources"

# Basename inside the ZIP -> output filename (P1 / P2 / P3 per BLUEPRINT.md §0)
MAPPING = {
    "Galois hulls of constacyclic codes over finite fields.pdf": (
        "P1_galois_hulls_constacyclic.txt",
        "P1",
    ),
    "Galios Hull over Affin Algebra.pdf": (
        "P2_affine_algebra.txt",
        "P2",
    ),
    "Average_dimensions_of_Galois_hulls_of_constacyclic.pdf": (
        "P3_average_dimensions.txt",
        "P3",
    ),
}

# Page counts asserted by BLUEPRINT.md §0 (17 / 27 / 36). Verified, not assumed.
EXPECTED_PAGES = {"P1": 17, "P2": 27, "P3": 36}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if not ZIP.is_file():
        print(f"ERROR: primary source ZIP not found: {ZIP}", file=sys.stderr)
        return 2

    OUTDIR.mkdir(parents=True, exist_ok=True)

    ziphash = sha256_of(ZIP)
    print(f"ZIP      : {ZIP.name}")
    print(f"ZIP sha256: {ziphash}\n")

    ok = True
    with zipfile.ZipFile(ZIP) as zf:
        members = {Path(n).name: n for n in zf.namelist()}
        for src_name, (out_name, tag) in MAPPING.items():
            if src_name not in members:
                print(f"ERROR: missing member in ZIP: {src_name}", file=sys.stderr)
                ok = False
                continue

            data = zf.read(members[src_name])
            doc = pymupdf.open(stream=data, filetype="pdf")
            npages = doc.page_count

            parts = []
            for i in range(npages):
                text = doc[i].get_text("text")
                parts.append(f"<<<PAGE {i + 1}>>>\n{text}")
            doc.close()

            out_path = OUTDIR / out_name
            body = "\n".join(parts)
            out_path.write_text(body, encoding="utf-8")

            expected = EXPECTED_PAGES[tag]
            status = "OK" if npages == expected else "MISMATCH"
            if npages != expected:
                ok = False

            print(
                f"{tag}  {out_name}\n"
                f"    pages      : {npages} (BLUEPRINT.md §0 asserts {expected}) [{status}]\n"
                f"    chars      : {len(body):,}\n"
                f"    sha256     : {sha256_of(out_path)}"
            )

    print(f"\nRESULT   : {'ALL OK' if ok else 'FAILURES PRESENT'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
