#!/usr/bin/env python3
"""Fetch real public inputs and record actual hashes; no synthetic substitutes."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import stat
import time
import urllib.request
import zipfile
from project_common import (pipeline, module_config, load, data_path, repo_path,
                            writable, resolved_b_config, sha, dump, now)


def fetch(url, path, expected_bytes=None, expected_sha256=None):
    if path.exists():
        status, elapsed = "REUSED_PRESENT_INPUT", None
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_name(path.name + ".part")
        if partial.exists():
            raise RuntimeError("Partial input exists; inspect before another bounded transfer: " + str(partial))
        started = time.time()
        request = urllib.request.Request(url, headers={"User-Agent": "IVDD-own-code-reproducibility/1"})
        with urllib.request.urlopen(request, timeout=90) as response, partial.open("wb") as stream:
            advertised = response.headers.get("Content-Length")
            for block in iter(lambda: response.read(1024 * 1024), b""):
                stream.write(block)
        if advertised and partial.stat().st_size != int(advertised):
            raise ValueError("Incomplete Content-Length for " + url)
        if expected_bytes and partial.stat().st_size != int(expected_bytes):
            raise ValueError("Frozen byte count differs for " + url)
        if expected_sha256 and sha(partial) != expected_sha256:
            raise ValueError("Frozen SHA256 differs for " + url)
        partial.replace(path)
        status, elapsed = "DOWNLOADED", time.time() - started
    digest = sha(path)
    if expected_bytes and path.stat().st_size != int(expected_bytes):
        raise ValueError("Existing input byte count differs: " + str(path))
    if expected_sha256 and digest != expected_sha256:
        raise ValueError("Existing input SHA256 differs: " + str(path))
    return {"url": url, "path": str(path), "bytes": path.stat().st_size, "sha256": digest,
        "expected_sha256": expected_sha256, "checksum_basis": "FROZEN_EXPECTED_SHA256" if expected_sha256 else "ACTUAL_HASH_RECORDED_NO_PRETENDED_FROZEN_SHA",
        "status": status, "seconds": elapsed}


def verify_source(source, lock):
    records = []
    for relative, expected in lock["git_blob_sha1"].items():
        path = source / relative
        raw = path.read_bytes()
        actual = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        if actual != expected:
            raise ValueError("Official source differs from pinned Git blob: " + relative)
        records.append({"relative_path": relative, "git_blob_sha1": actual, "sha256": sha(path)})
    return records


def install_source(archive, destination, commit):
    """Extract the pinned archive under a new private runtime root, never repo code."""
    prefix = "SigXTalk-" + commit
    staging = destination.with_name(destination.name + ".extracting")
    if staging.exists():
        raise RuntimeError("Preserve and inspect an earlier extraction: " + str(staging))
    staging.mkdir(parents=True)
    with zipfile.ZipFile(archive) as handle:
        for item in handle.infolist():
            parts = Path(item.filename).parts
            if not parts or parts[0] != prefix or ".." in parts or stat.S_ISLNK(item.external_attr >> 16):
                raise ValueError("Unexpected or unsafe pinned archive member")
            relative = Path(*parts[1:])
            if str(relative) == ".":
                continue
            target = (staging / relative).resolve()
            if not target.is_relative_to(staging.resolve()):
                raise ValueError("Archive path escapes its private staging directory")
            if item.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with handle.open(item) as source, target.open("wb") as out:
                    shutil.copyfileobj(source, out)
    staging.replace(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--project-root")
    parser.add_argument("--dataset", choices=["all", "human", "spatial", "official"], default="all")
    parser.add_argument("--include-official-source", action="store_true")
    parser.add_argument("--include-official-demo", action="store_true")
    parser.add_argument("--from-cached-inputs", action="store_true", help="Record skipped downloads; relevant modules verify only their needed actual caches")
    args = parser.parse_args()
    cfg = pipeline(args.config, args.project_root)
    root = cfg["_root"]
    records = []
    status = {"started_utc": now(), "records": records, "status": "RUNNING"}
    output = root / "00_audit/PUBLIC_DOWNLOAD_RECEIPTS.json"
    if args.from_cached_inputs:
        dump({"status": "SKIPPED_CACHED_MODE_DOWNLOADS_NOT_RUN", "finished_utc": now(),
            "raw_network_transfers": 0, "cache_validation": "actual selected inputs are verified by their relevant modules; unused four libraries, Tu raw files and spatial tar are not required"}, root / "00_audit/PUBLIC_DOWNLOAD_CACHE_SKIP.json")
        print("Cached mode: public downloads skipped; stage-specific actual input checks remain required.")
        return
    def download(url, path, **kwargs):
        if not path.exists():
            writable(root, path)
        receipt = fetch(url, path, **kwargs)
        records.append(receipt)
        dump(status, output)
        print(json.dumps(receipt), flush=True)
    try:
        if args.dataset in ["all", "human"]:
            raw = load(module_config(cfg, "raw_selection")).get("raw_selection", {})
            manifest = repo_path(raw.get("sample_manifest", "resources/metadata/SIGXTALK_GSE230808_SOURCE_MANIFEST.tsv"))
            raw_dir = data_path(root, raw.get("raw_dir", "data_cache/GSE230808"))
            with manifest.open(encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            if len(rows) != 16:
                raise ValueError("Expected all sixteen public disease AF/NP libraries")
            for row in rows:
                for suffix in ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]:
                    url = row["matrix_url"].replace("matrix.mtx.gz", suffix)
                    name = url.rsplit("/", 1)[-1]
                    candidates = [p for p in [raw_dir / name, raw_dir / row["gsm"] / name] if p.exists()]
                    if len(candidates) > 1:
                        raise ValueError("Ambiguous raw GEO cache layout")
                    download(url, candidates[0] if candidates else raw_dir / row["gsm"] / name)
            human = load(module_config(cfg, "human"))
            tu_manifest = repo_path(human.get("resources", {}).get("tu_downloads", "resources/metadata/TU_DOWNLOAD_URLS_AND_SIZES.tsv"))
            with tu_manifest.open(encoding="utf-8") as stream:
                tu_rows = list(csv.DictReader(stream, delimiter="\t"))
            if len(tu_rows) != 8:
                raise ValueError("Expected eight Tu released matrices and cell-name companions")
            tu_dir = data_path(root, human["tu"]["raw_dir"])
            for row in tu_rows:
                for kind in ["counts", "cellname"]:
                    download(row[kind + "_url"], tu_dir / row[kind + "_file"], expected_bytes=int(row[kind + "_size_bytes"]))
        if args.dataset in ["all", "spatial"]:
            spatial = load(module_config(cfg, "spatial"))
            download(spatial["raw_url"], data_path(root, spatial["raw_tar"]),
                expected_bytes=spatial["raw_bytes"], expected_sha256=spatial["raw_sha256"])
        downloads = cfg.get("downloads", {})
        if args.dataset in ["all", "official"]:
            b_path, b = resolved_b_config(cfg)
            if args.include_official_source or downloads.get("include_official_source", False):
                lock = load(repo_path("scripts/sigxtalk_upstream_lock.json"))
                source = data_path(root, b["official_source"])
                if not source.exists():
                    writable(root, source)
                    archive = root / "data_cache" / ("SigXTalk_" + lock["commit"] + ".zip")
                    download("https://codeload.github.com/LithiumHou/SigXTalk/zip/" + lock["commit"], archive)
                    install_source(archive, source, lock["commit"])
                proof = {"repository": lock["repository"], "commit": lock["commit"],
                    "critical_files": verify_source(source, lock), "not_vendored_in_own_code_repository": True}
                dump(proof, root / "00_audit/OFFICIAL_SOURCE_DOWNLOAD_PROOF.json")
            if args.include_official_demo or downloads.get("include_official_demo", False):
                download(downloads.get("demo_url", "https://zenodo.org/api/records/15531285/files/SigXTalk_demo_data.rds/content"),
                    data_path(root, b["pbmc_rds"]), expected_sha256=downloads.get("demo_sha256"))
        status["status"] = "DONE_ACTUAL_PUBLIC_INPUTS"
    except Exception as error:
        status["status"], status["error"] = "PARTIAL_FAILED", str(error)
        raise
    finally:
        status["ended_utc"] = now()
        dump(status, output)


if __name__ == "__main__":
    main()
