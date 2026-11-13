"""Download and verify a dated snapshot of ChEMBL EGFR IC50 activities."""
import csv
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urljoin
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[2]
API = "https://www.ebi.ac.uk/chembl/api/data/"
FIELDS = (
    "activity_id molecule_chembl_id canonical_smiles standard_type standard_relation "
    "standard_value standard_units assay_chembl_id target_chembl_id document_chembl_id "
    "assay_type assay_description target_organism target_pref_name pchembl_value "
    "data_validity_comment potential_duplicate activity_comment type relation value units"
).split()
QUERY = API + "activity.json?" + urlencode({
    "target_chembl_id": "CHEMBL203", "standard_type": "IC50", "order_by": "activity_id",
    "limit": 500, "only": ",".join(FIELDS),
})


def get_json(url):
    for attempt in range(5):
        try:
            with urlopen(url, timeout=120) as response:
                return json.load(response)
        except HTTPError as exc:
            if exc.code not in {408, 429, 500, 502, 503, 504} or attempt == 4:
                raise
        except (URLError, TimeoutError, ConnectionError) as exc:
            if attempt == 4:
                raise
            print(f"Retrying request after {exc}", file=sys.stderr, flush=True)
        time.sleep(2 ** attempt)


def download(folder):
    url = QUERY
    visited = set()
    total = None
    csv_path = folder / "egfr_ic50.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, extrasaction="ignore")
        writer.writeheader()
        count = 0
        while url:
            if url in visited:
                raise ValueError("Repeated pagination link")
            visited.add(url)
            page = get_json(url)
            meta = page["page_meta"]
            if total is None:
                total = meta["total_count"]
            if meta["offset"] != count or meta["total_count"] != total:
                raise ValueError("Page offset or reported total changed during retrieval")
            for row in page["activities"]:
                if row["target_chembl_id"] != "CHEMBL203" or row["standard_type"] != "IC50":
                    raise ValueError("Received activity outside requested query")
                writer.writerow(row)
                count += 1
            print(f"Page {len(visited)}: {count}/{total}", flush=True)
            url = urljoin(url, meta["next"]) if meta["next"] else None

    ids, compounds, assays = set(), set(), set()
    with csv_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            aid = row["activity_id"]
            if not aid or aid in ids:
                raise ValueError(f"Missing or repeated activity ID: {aid!r}")
            ids.add(aid)
            compounds.add(row["molecule_chembl_id"])
            assays.add(row["assay_chembl_id"])
    if not total or len(ids) != total or count != total:
        raise ValueError(f"Incomplete snapshot: {len(ids)} distinct activities vs {total} expected")
    return {
        "query": QUERY, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "pages": len(visited), "reported_total": total, "verified_unique_activity_ids": len(ids),
        "distinct_compounds_with_id": len(compounds - {""}), "distinct_assays_with_id": len(assays - {""}),
        "sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
    }


def main():
    raw = ROOT / "data" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=raw) as temp:
        staging = Path(temp) / "snapshot"
        staging.mkdir()
        manifest = download(staging)
        (staging / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        dest = raw / datetime.now(timezone.utc).strftime("egfr_ic50_%Y%m%dT%H%M%S%fZ")
        staging.rename(dest)
    print(dest)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
