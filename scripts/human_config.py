"""Portable configuration helpers for the human-data module."""
from pathlib import Path
import json

REPO = Path(__file__).resolve().parents[1]


def load_config(path, project_root=None):
    path = Path(path).resolve()
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        config = json.loads(text)
    else:
        import yaml
        config = yaml.safe_load(text)
    root = Path(project_root or config.get("project_root", "."))
    config["_root"] = (root if root.is_absolute() else REPO / root).resolve()
    config["_repo"] = REPO
    return config


def data_path(config, value):
    path = Path(value)
    return path if path.is_absolute() else config["_root"] / path


def resource(config, key, default):
    path = Path(config.get("resources", {}).get(key, default))
    return path if path.is_absolute() else config["_repo"] / path

