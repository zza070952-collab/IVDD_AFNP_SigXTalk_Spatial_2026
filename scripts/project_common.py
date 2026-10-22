"""Portable paths, hashes and checked subprocesses for the own-code workflow."""
from __future__ import annotations
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[1]


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def dump(value, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(value, indent=2, allow_nan=False) + "\n"
    if not path.exists() or path.read_text(encoding="utf-8") != raw:
        path.write_text(raw, encoding="utf-8")


def load(path):
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        return json.loads(text)
    import yaml
    return yaml.safe_load(text)


def repo_path(value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPO / path).resolve()


def data_path(root, value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def writable(root, path):
    path = path.resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("Runtime output must remain in the new project root: " + str(path))
    return path


def pipeline(path, project_root=None):
    cfg = load(path)
    cfg["_root"] = data_path(REPO, project_root or cfg.get("project_root", "analysis"))
    cfg["_config"] = Path(path).resolve()
    return cfg


def module_config(cfg, name):
    return repo_path(cfg.get("modules", {}).get(name + "_config", {
        "human": "config/human.example.yaml", "raw_selection": "config/raw_selection.pipeline.example.json",
        "sigxtalk": "configs/sigxtalk.config.example.json", "spatial": "configs/spatial.example.yaml",
        "integration": "config/integration.example.yaml"}[name]))


def science_python(cfg):
    return cfg.get("executables", {}).get("science_python", sys.executable)


def execute(cfg, name, argv):
    root = cfg["_root"]
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    logfile = root / "logs" / ("project_" + name + "_" + stamp + ".log")
    logfile.parent.mkdir(parents=True, exist_ok=True)
    started, clock = now(), time.time()
    env = os.environ.copy()
    addons = cfg.get("executables", {}).get("python_addons")
    if addons:
        env["PYTHONPATH"] = str(data_path(root, addons)) + os.pathsep + env.get("PYTHONPATH", "")
    with logfile.open("w", encoding="utf-8") as stream:
        stream.write("# actual argv=" + json.dumps(list(map(str, argv))) + "\n")
        stream.flush()
        result = subprocess.run(list(map(str, argv)), cwd=REPO, env=env, stdout=stream, stderr=subprocess.STDOUT)
    receipt = {"name": name, "started_utc": started, "ended_utc": now(),
        "argv": list(map(str, argv)), "exit_code": result.returncode,
        "seconds": time.time() - clock, "log_sha256": sha(logfile), "log": str(logfile)}
    dump(receipt, root / "00_audit" / ("EXECUTION_" + name + ".json"))
    dump(receipt, logfile.with_suffix(".execution.json"))
    if result.returncode:
        raise RuntimeError(f"{name} exited {result.returncode}; inspect {logfile}")
    return receipt


def require_files(paths):
    records = []
    for path in map(Path, paths):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("Required measured output absent/empty: " + str(path))
        records.append({"path": str(path), "bytes": path.stat().st_size, "sha256": sha(path)})
    return records


def resolved_b_config(cfg):
    """Preserve an already-correct frozen B config byte-for-byte on cache reuse."""
    original = module_config(cfg, "sigxtalk")
    content = load(original)
    old_root = data_path(original.parent, content["project_root"])
    if old_root == cfg["_root"]:
        return original, content
    content["project_root"] = str(cfg["_root"])
    runtime = cfg["_root"] / "00_audit/runtime_configs/sigxtalk.resolved.json"
    dump(content, runtime)
    return runtime, content
