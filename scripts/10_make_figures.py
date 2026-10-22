#!/usr/bin/env python3
"""Regenerate four main figures only from completed actual module outputs."""
import argparse
from project_common import *
from run_project import spatial_context, verify_completed_depth, verify_b_completion


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--project-root")
    parser.add_argument("--from-raw", action="store_true", help="Use the new raw-mode spatial baseline cache")
    args = parser.parse_args()
    cfg = pipeline(args.config, args.project_root)
    root = cfg["_root"]
    values, baseline = spatial_context(cfg, args.from_raw)
    if verify_completed_depth(values, baseline) is None or verify_b_completion(cfg) is None:
        raise RuntimeError("Figure export requires completed new B and the single fixed-depth C check")
    b_path, _ = resolved_b_config(cfg)
    steps = [execute(cfg, "figures12", [science_python(cfg), repo_path("scripts/10_make_human_figures.py"),
        "--config", module_config(cfg, "human"), "--project-root", root])]
    steps.append(execute(cfg, "figure3", [science_python(cfg), repo_path("scripts/10_plot_new_sigxtalk_figure3.py"),
        "--config", b_path, "--out", root / "05_figures"]))
    steps.append(execute(cfg, "figure4", [science_python(cfg), repo_path("scripts/plot_figure4.py"),
        "--spatial-input", baseline, "--depth-input", values["depth_output"],
        "--roi-json", values["roi_config"], "--out", root / "05_figures"]))
    stems = [root / "05_figures/Figure1_Figure2/Figure1_Adult_Cohorts_and_Identity",
        root / "05_figures/Figure1_Figure2/Figure2_Fixed_NP_Programs_CrossCohort",
        root / "05_figures/Figure3_New_Official_SigXTalk", root / "05_figures/Figure4_Adult_Mouse_Spatial_Context"]
    output = require_files([str(stem) + "." + ext for stem in stems for ext in ["png", "pdf", "svg"]])
    dump({"status": "DONE_FOUR_MAIN_FIGURES_FROM_ACTUAL_OUTPUTS", "finished_utc": now(),
        "main_figure_count": 4, "executions": steps, "outputs": output,
        "visual_QA": "REQUIRES_ACTUAL_IMAGE_INSPECTION_NOT_AUTOMATICALLY_PASSED",
        "new_experiment_figure_created": False}, root / "05_figures/FIGURE_EXPORT_STATUS.json")


if __name__ == "__main__":
    main()
