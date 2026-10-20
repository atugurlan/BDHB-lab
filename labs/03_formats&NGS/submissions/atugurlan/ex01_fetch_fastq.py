"""
Exercise 03 — FASTQ Download (student-owned)

Objective:
- Choose a TP53-related accession (e.g., SRR..., ERR...) and DOWNLOAD a FASTQ file.
- Save it in data/work/<handle>/lab03/your_reads.fastq.gz

Minimum requirements:
- The script must accept an accession (e.g., from command line args).
- The script downloads at least one FASTQ (one file is enough for the exercise).
- The script prints to stdout the path of the downloaded file.

Recommended:
- Support .fastq or .fastq.gz.

NOTE:
- The chosen library does not matter (requests/urllib/etc.), but avoid heavy packages.
"""

import argparse
from pathlib import Path
import requests


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    ap.add_argument("--handle", required=True)
    args = ap.parse_args()

    response = requests.get(
        "https://www.ebi.ac.uk/ena/portal/api/filereport",
        params={
            "accession": args.accession,
            "result": "read_run",
            "fields": "fastq_ftp",
            "format": "json"
        },
        timeout=30
    )
    records = response.json()

    if not records or not records[0]["fastq_ftp"]:
        raise ValueError("No FASTQ file is available for this accession.")

    fastq_link = records[0]["fastq_ftp"].split(";")[0]
    url = "https://" + fastq_link.removeprefix("ftp://")

    out_path = Path(f"data/work/{args.handle}/lab03/your_reads.fastq.gz")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(url, timeout=60)
    response.raise_for_status()
    out_path.write_bytes(response.content)

    print("Downloaded:", out_path)


if __name__ == "__main__":
    main()
