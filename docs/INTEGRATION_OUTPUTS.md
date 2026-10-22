The integration generator uses measured source tables and actual new formal
execution receipts. It produces candidate evidence, negative evidence, a
conditional future-experiment decision and manuscript text. It does not train
a communication model, run another spatial sampling analysis, invent a clinical
test, or calculate a combined mechanism score.

```bash
python scripts/09_integrate_evidence.py --project-root "$PROJECT" \
  --config config/integration.example.yaml \
  --tu-candidate-table "$PROJECT/01_human_programs/TU_CANDIDATE_GENE_BY_DONOR_IDENTITY.tsv" \
  --spatial-input "$BASELINE_SPATIAL" --sigxtalk-config "$ACTUAL_RUN_CONFIG"
```

All example paths are relative to the analysis project root. Human candidate
input must retain each of the ten prespecified genes across all actual donors,
with unavailable features explicit. Spatial candidate rows must cover each
actual disc/region. None of these missing observations is replaced by positive
expression. Independent within-donor identity direction is stated separately
from unresolved clinical grade.

Default mode requires NEW_R_PYTHON_R_COMPLETION.json, the formal final status,
HGNN training status and raw GSE230808 input provenance. The source is required
to be formal disease counts, with no old pathways/PRS used. Prediction and
input-provenance hashes must agree across receipts; all-tested and positive PRS
row counts must reproduce the actual complete tables. A smoke run cannot pass
the formal gate. The ten-candidate source retains every tested PA/PRS row and
NA. Target absence, PA failure, PRS nonreturn, relative filtering and absence
of source-supported receptor context have separate reasons. Missing/unrun
quantities remain NA rather than being represented as zero tests.

`--allow-incomplete-b` permits an explicitly PARTIAL pending draft from genuine
completed human/spatial tables. It does not consume partially written or smoke
SigXTalk tables as formal results and cannot overwrite a completed integration
with a pending draft. Its new formal pathway counts remain unassessed. Strict
mode should be used by the final top-level workflow.

The main candidate matrix contains fixed GEP2 membership/rank, per-donor NP
expression/detection and within-identity direction, original AF ligand LR
complexes, explicit complex component roles, receiver Receptor–SSC–target keys,
original PA/PRS and applicable official metrics, per-disc mouse detectability,
depth attenuation, evidence level and counterevidence. Representative original
rows are selected only for compact display by decreasing official PRS, then PA,
then LR probability and stable original keys; this display order is not a new
score. Complete candidate axes and joined rows remain available separately.

FN1 is evaluated as a sender-ligand candidate and ITGB1 as an integrin-complex
component; neither is required to be a selected downstream target. PLOD2 and
all prespecified comparison genes retain independent target testability. A
passing FN1/ITGB1-containing receptor-key context with PLOD2 can retain the
factorial perturbation as a hypothesis. If that context is absent, the decision
does not present it as a supported causal chain and lists only genuine passing
alternative target contexts, without selecting one through a custom score.

Terminology is fixed throughout generated text:

| Term | Scope |
|---|---|
| GEP2 | Frozen human NP-associated program, not a disease grade |
| CellChat LR probability | Original inferred AF-to-NP ligand–complex context |
| HGNN PA | Original receiver-axis activation prediction |
| official PRS | Original returned official pathway relevance quantity |
| SSC | Database-listed intermediate; source column TF is retained for audit |
| receptor-key context | LR/component linkage to receiver model; not ligand-specific induction |
| donor / disc nested in capture | Disclosed descriptive units; cells/bins are not independent patients/animals |

Methods and Results use actual A/B/C table values and recorded run parameters.
R versions are read from the formal R_SESSION receipt when present; missing
version provenance is not filled from memory. Only existing verified figure
legends are incorporated, so run figure generation before the final integration.
No Figure5 experimental data or old THBS1 validation is presented as new work.

Code Availability defaults to BLOCKED_UNVERIFIED. To enable an actual public
statement, supply `--github-proof` pointing to independently verified JSON with
`verified_public: true`, `repo_url`, a full forty-character `head_sha`, and
`verified_utc`. Preparing a ZIP or knowing a suggested URL is insufficient.
Data Availability distinguishes public reused GEO data, the local review
artifact and any still-unassigned permanent derived-data deposition.

Dependencies: Python, numpy, pandas, PyYAML. Five tests exercise strict rejection,
partial/full document plumbing in disposable clearly marked synthetic structure
fixtures, target-specific negative/NA reasons, and public-release proof gating.
Technical fixture output is deleted automatically and never used as research
evidence. Real-run execution logs and artifact hashes must accompany final
manuscript claims. Source code preparation alone is not proof of formal analysis.
