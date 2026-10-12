"""Download the course-specific negative dataset, then prepare its outputs.

Requires internet access and Python's standard library. Existing raw inputs
are backed up before replacement. Search pagination supplies the authoritative
accession set; FASTA stream exports are reconciled to that same set.
"""

import io
import json
import re
import shutil
import tempfile
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import filter_negatives as negative


def search_url():
    return "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({
        "query": negative.QUERY, "format": "tsv",
        "fields": negative.EXPORT_FIELDS, "size": 500,
    })


def open_url(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers={
        "User-Agent": "LB2-Group9-negative-collection/1.0",
        "Accept-Encoding": "identity",
    }), timeout=60)


def collect(root=negative.ROOT):
    import csv
    raw = root / "data/raw"
    raw.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="negative_download_", dir=raw) as temp:
        stage = Path(temp)
        url = search_url()
        rows, fields = [], None
        expected_count = release = release_date = None
        page = 0
        while url:
            with open_url(url) as response:
                count = int(response.headers["X-Total-Results"])
                current_release = response.headers.get("X-UniProt-Release")
                if not current_release:
                    raise ValueError("UniProt release header is missing")
                if expected_count is None:
                    expected_count, release = count, current_release
                    release_date = response.headers.get("X-UniProt-Release-Date")
                elif (count, current_release) != (expected_count, release):
                    raise ValueError("UniProt search changed during collection; retry")
                reader = csv.DictReader(io.StringIO(response.read().decode("utf-8")), delimiter="\t")
                if fields is None:
                    fields = reader.fieldnames
                elif fields != reader.fieldnames:
                    raise ValueError("TSV columns changed during collection")
                rows.extend(reader)
                link = response.headers.get("Link")
            url = None
            if link:
                match = re.search(r'<([^>]+)>;\s*rel="next"', link)
                if not match or not match.group(1).startswith("https://rest.uniprot.org/"):
                    raise ValueError("Unexpected UniProt pagination link")
                url = match.group(1)
            page += 1
            if page % 5 == 0 or not url:
                print(f"Downloaded {len(rows):,} / {expected_count:,} annotations", flush=True)
        if len(rows) != expected_count:
            raise ValueError("Incomplete annotation download")
        accessions = [row["Entry"] for row in rows]
        if len(set(accessions)) != len(accessions):
            raise ValueError("Duplicate search accession")
        fasta_url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode({
            "query": negative.QUERY, "format": "fasta",
        })
        print("Downloading canonical protein sequences", flush=True)
        stream_path = stage / "stream.fasta"
        with open_url(fasta_url) as response, stream_path.open("wb") as handle:
            if response.headers.get("X-UniProt-Release") != release:
                raise ValueError("UniProt release changed before sequence download")
            shutil.copyfileobj(response, handle)
        streamed = negative.read_fasta(stream_path)
        if not set(accessions).issubset(streamed):
            raise ValueError("Some search accessions have no downloaded sequence")
        extra = sorted(set(streamed) - set(accessions))
        records = {accession: streamed[accession] for accession in accessions}
        negative.validate_candidates(fields, rows, records)
        negative.positive_overlap(root, accessions)
        negative.write_tsv(stage / "negative_candidates.tsv", fields, rows)
        negative.write_fasta(stage / "negative_candidates.fasta", records, accessions)
        manifest = {
            "collected_utc": datetime.now(timezone.utc).isoformat(),
            "uniprot_release": release, "uniprot_release_date": release_date,
            "query": negative.QUERY, "export_fields": negative.EXPORT_FIELDS,
            "locations": negative.LOCATIONS, "search_result_count": len(rows),
            "accession_set_source": "UniProtKB search endpoint, complete pagination",
            "accession_set_sha256": negative.accession_digest(accessions),
            "sha256": {name: negative.sha256(stage / name) for name in
                       ["negative_candidates.tsv", "negative_candidates.fasta"]},
            "stream_accessions_outside_search": extra,
        }
        (stage / "negative_collection.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        names = ["negative_candidates.tsv", "negative_candidates.fasta", "negative_collection.json"]
        existing = [raw / name for name in names if (raw / name).exists()]
        if existing:
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            backup = root / "data/archive" / (stamp + "_negative_raw")
            backup.mkdir(parents=True)
            for path in existing:
                shutil.copy2(path, backup / path.name)
        for name in names:
            shutil.copy2(stage / name, raw / name)
    negative.build(root)


if __name__ == "__main__":
    try:
        collect()
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f"Collection stopped: {error}")
