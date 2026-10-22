#!/usr/bin/env python3
"""Download exact public GEO files with bounded transfers and SHA256 receipts."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import urllib.request

import pandas as pd
from human_config import load_config, data_path, resource


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def fetch(url, path, expected_size=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if expected_size and path.stat().st_size != expected_size:
            raise ValueError(f"Existing cache size mismatch; preserve and inspect: {path}")
        return {"url": url, "path": str(path), "bytes": path.stat().st_size,
                "sha256": digest(path), "status": "REUSED_PRESENT_SIZE_HASH_RECORDED"}
    partial = path.with_name(path.name + ".part")
    if partial.exists():
        raise RuntimeError(f"Partial transfer exists; inspect it before a bounded retry: {partial}")
    request = urllib.request.Request(url, headers={"User-Agent": "IVDD-public-reproducibility/1"})
    started = time.time()
    checksum = hashlib.sha256()
    with urllib.request.urlopen(request, timeout=90) as response, partial.open("wb") as handle:
        advertised = response.headers.get("Content-Length")
        for block in iter(lambda: response.read(1024 * 1024), b""):
            handle.write(block)
            checksum.update(block)
    size = partial.stat().st_size
    expected = expected_size or (int(advertised) if advertised else None)
    if expected and size != expected:
        raise ValueError(f"Incomplete transfer: {url}, {size} versus {expected}")
    partial.replace(path)
    return {"url": url, "path": str(path), "bytes": size, "sha256": checksum.hexdigest(),
            "status": "DOWNLOADED", "seconds": round(time.time() - started, 2)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--project-root")
    parser.add_argument("--dataset", choices=["tu", "discovery", "all"], default="all")
    args = parser.parse_args()
    config = load_config(args.config, args.project_root)
    records = []
    try:
        if args.dataset in ["tu", "all"]:
            manifest = pd.read_csv(resource(config, "tu_downloads", "resources/metadata/TU_DOWNLOAD_URLS_AND_SIZES.tsv"), sep="\t")
            destination = data_path(config, config["tu"]["raw_dir"])
            for row in manifest.to_dict("records"):
                for kind in ["counts", "cellname"]:
                    receipt = fetch(row[kind + "_url"], destination / row[kind + "_file"], int(row[kind + "_size_bytes"]))
                    receipt.update(dataset="GSE165722", gsm=row["gsm"], tissue="NP")
                    records.append(receipt)
                    print(json.dumps(receipt), flush=True)
        if args.dataset in ["discovery", "all"]:
            manifest = pd.read_csv(resource(config, "discovery_samples", "resources/metadata/DISCOVERY_SAMPLE_MANIFEST.tsv"), sep="\t")
            manifest = manifest[(manifest.geo_series == "GSE230808") & (manifest.tissue == "NP")]
            destination = data_path(config, config["discovery"]["phase2_workspace"]) / "02_geo_downloads" / "GSE230808"
            for row in manifest.to_dict("records"):
                url = row["matrix_url"]
                for suffix in ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]:
                    current = url.replace("matrix.mtx.gz", suffix)
                    receipt = fetch(current, destination / current.rsplit("/", 1)[-1])
                    receipt.update(dataset="GSE230808", gsm=row["gsm"], donor=row["donor_id"], tissue="NP")
                    records.append(receipt)
                    print(json.dumps(receipt), flush=True)
    finally:
        audit = config["_root"] / "00_audit"
        audit.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(records).to_csv(audit / "HUMAN_DOWNLOAD_RECEIPTS.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
