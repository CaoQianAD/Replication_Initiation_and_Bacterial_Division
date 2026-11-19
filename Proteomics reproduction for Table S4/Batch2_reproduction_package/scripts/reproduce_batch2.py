"""Rebuild the manuscript's Batch2 protein tables from an original DIA-NN report."""

import argparse
import hashlib
import json
import os
import platform
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow
import pyarrow.parquet as pq


GROUP = "Protein.Group"
MEASURES = ("MaxLFQ_linear", "MaxLFQ_log2", "PPM")
EXPERIMENTS = ("dnaA", "hns")
GROUP_COLUMNS = [GROUP, "Protein.Ids", "Genes", "Protein.Names"]
TARGET = {"ECOLI", "ALLOWED_CUSTOM"}
ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def write(frame, path, index=False, missing="NA"):
    frame.to_csv(path, index=index, encoding="utf-8-sig", na_rep=missing)


def stable_join(values):
    result = {}
    for value in values:
        if pd.isna(value):
            continue
        for token in re.split(r"\s*;\s*|\s*//\s*", str(value)):
            token = token.strip()
            if token and token.upper() != "NA":
                result[token] = None
    return ";".join(result) if result else pd.NA


def group_category(accessions, source):
    members = set()
    for token in re.split(r"[;,]", accessions):
        token = token.strip()
        if token:
            members.add(token.split("|")[1] if token.count("|") >= 2 else token)
    categories = {source.get(member, "UNKNOWN") for member in members}
    if "UNKNOWN" in categories:
        return "UNKNOWN"
    if "DECOY" in categories:
        return "DECOY"
    if "CONTAMINANT" in categories:
        return "MIXED_TARGET_CONTAMINANT" if categories & TARGET else "CONTAMINANT"
    return "TARGET_ELIGIBLE" if categories and categories <= TARGET else "UNKNOWN"


def classify_batch2(report, accessions):
    # Stage 2 derives identity from all report rows, before run and q-value selection.
    info = pq.read_table(report, columns=GROUP_COLUMNS).to_pandas().drop_duplicates()
    info = info.groupby(GROUP, sort=True)[GROUP_COLUMNS[1:]].agg(stable_join).reset_index()
    ids = info["Protein.Ids"].fillna(info[GROUP])
    info["Final_Source_Class"] = [group_category(value, accessions) for value in ids]
    return info


def prepare_quantification(report, manifest):
    selected = manifest[
        manifest["Run_Type"].eq("BIOLOGICAL") & manifest["Frozen_Run_Decision"].eq("Include")
    ]
    if not selected["Run_ID"].is_unique or not selected["Sample"].is_unique:
        raise ValueError("This Batch2 workflow expects one unique retained run per sample")
    run_to_sample = selected.set_index("Run_ID")["Sample"]
    rows = pq.read_table(report, columns=["Run", GROUP, "PG.MaxLFQ", "Global.PG.Q.Value"]).to_pandas()
    if missing := set(run_to_sample.index) - set(rows["Run"].unique()):
        raise ValueError(f"Runs absent from the DIA-NN report: {sorted(missing)}")
    rows = rows[rows["Run"].isin(run_to_sample.index)].copy()
    intensity = pd.to_numeric(rows["PG.MaxLFQ"], errors="coerce").astype(np.float64)
    qvalue = pd.to_numeric(rows["Global.PG.Q.Value"], errors="coerce")
    rows["PG.MaxLFQ"] = intensity
    valid = rows[qvalue.lt(.01) & np.isfinite(intensity) & intensity.gt(0)]
    quantity = valid.groupby(["Run", GROUP], as_index=False)["PG.MaxLFQ"].max()
    quantity["Sample"] = quantity["Run"].map(run_to_sample)
    return quantity, selected["Sample"].tolist()


def quantify(quantity, groups, samples, audit_dir):
    info = groups.set_index(GROUP)
    if not info.index.is_unique:
        raise ValueError("Protein group classification has duplicate IDs")
    missing = set(quantity[GROUP]) - set(info.index)
    if missing:
        raise ValueError(f"{len(missing)} quantified groups lack classification")
    selected = quantity[quantity[GROUP].isin(info.index[info["Final_Source_Class"].eq("TARGET_ELIGIBLE")])].copy()
    selected["MaxLFQ_log2"] = np.log2(selected["PG.MaxLFQ"])
    # The historical numerical path computed 2**log2(MaxLFQ), also for single runs.
    selected["MaxLFQ_linear"] = np.exp2(selected["MaxLFQ_log2"])
    denominators = selected.groupby("Sample")["MaxLFQ_linear"].sum().reindex(samples)
    if denominators.isna().any() or not denominators.gt(0).all():
        raise ValueError("A retained sample has no positive target denominator")
    selected["PPM"] = selected["MaxLFQ_linear"] / selected["Sample"].map(denominators) * 1_000_000
    ppm_sums = selected.groupby("Sample")["PPM"].sum().reindex(samples)
    if not np.allclose(ppm_sums, 1e6, rtol=0, atol=1e-6):
        raise ValueError("PPM normalization check failed")
    write(pd.DataFrame({"target_MaxLFQ_sum": denominators, "PPM_sum": ppm_sums}),
          audit_dir / "ppm_denominators.csv", index=True)
    write(selected, audit_dir / "target_run_level.csv")
    matrices = {key: selected.pivot(index=GROUP, columns="Sample", values=key).reindex(columns=samples)
                for key in MEASURES}
    return matrices


def sample_design(metadata, paper, retained):
    names = paper["Sample"].tolist()
    if paper["Sample"].duplicated().any() or paper["paper_label"].duplicated().any():
        raise ValueError("The manuscript sample mapping contains duplicates")
    if not set(names) <= set(retained):
        raise ValueError("A manuscript sample is absent from the retained Batch2 runs")
    for key in ["Group", "Biological_replicates", "Time_h", "Phase"]:
        mapped = metadata.set_index("Sample").loc[names, key].tolist()
        if mapped != paper[key].tolist():
            raise ValueError(f"Manuscript sample metadata differ for {key}")
    if paper["Group"].value_counts().to_dict() != {"dnaA_Upshift": 36, "hns_Upshift": 32}:
        raise ValueError("Expected 36 dnaA and 32 hns manuscript samples")


def stage2_filter(matrices, info, metadata, paper, experiment, mode_dir):
    # ABC columns used by the original Stage 2 are all in Batch2.
    design = paper[paper["Group"].eq(f"{experiment}_Upshift")].set_index("Sample")
    ppm = matrices["PPM"]
    linear = matrices["MaxLFQ_linear"]
    complete = linear[design.index].notna().all(axis=1)
    triggers = {}
    for rep in "ABC":
        names = design.index[design["Biological_replicates"].str.endswith(f"replicate{rep}")]
        data = ppm[names]
        triggers[rep] = data.max(axis=1).lt(10) & data.min(axis=1).lt(1)
    low = pd.DataFrame(triggers).any(axis=1)
    keep = complete & ~low
    audit = pd.DataFrame({"stage2_complete": complete, "stage2_low_expression": low, "stage2_keep": keep})
    write(audit, mode_dir / f"{experiment}_stage2_filter_audit.csv", index=True)
    # R selects Batch2 A/B/C from the experiment-wide matrices.
    names = metadata.loc[metadata["Group"].str.startswith(experiment) & metadata["Sample"].isin(ppm.columns), "Sample"]
    for key in ("MaxLFQ_log2", "PPM"):
        matrix = matrices[key].loc[keep, names].copy()
        matrix.insert(0, "gene", info["Genes"].reindex(matrix.index))
        write(matrix, mode_dir / f"{experiment}_{key}_filtered.csv", index=True)
    return int(keep.sum())


def run_r_model(rscript, r_library, mode_dir, out_root, metadata, experiment):
    model_dir = mode_dir / experiment
    model_dir.mkdir(exist_ok=True)
    env = {**os.environ}
    temp = str((out_root / "tmp").resolve())
    env.update({"TEMP": temp, "TMP": temp, "TMPDIR": temp})
    if r_library:
        env["BATCH2_R_LIBRARY"] = r_library.resolve().as_posix()
    command = [
        rscript, "--vanilla", Path(__file__).with_name("run_ppm_model.R").as_posix(), experiment,
        (mode_dir / f"{experiment}_MaxLFQ_log2_filtered.csv").as_posix(),
        (mode_dir / f"{experiment}_PPM_filtered.csv").as_posix(), metadata.as_posix(), model_dir.as_posix(),
    ]
    log = out_root / "logs" / f"{mode_dir.name}_{experiment}_R.log"
    with log.open("w", encoding="utf-8") as stream:
        result = subprocess.run(command, env=env, cwd=out_root, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{experiment} R model failed. See {log}")
    path = model_dir / "06_PPM" / "04_response" / f"{experiment}_PPM_track_formal_response_all.csv"
    model = pd.read_csv(path).set_index(GROUP)
    if not model.index.is_unique:
        raise ValueError(f"{experiment} model has duplicate protein groups")
    return model


def make_exports(matrices, groups, paper, selection, division, products, models):
    labels = paper["paper_label"].tolist()
    names = paper["Sample"].tolist()
    keep = matrices["PPM"][names].notna().any(axis=1)
    raw = {metric: matrices[metric].loc[keep, names] for metric in MEASURES}
    info = groups.set_index(GROUP).loc[raw["PPM"].index]
    division_genes = set(division.loc[division["Analysis_Used"].str.casefold().eq("true"), "gene"].str.strip().str.lower())
    product_map = products.assign(gene=products["gene"].str.strip().str.lower()).groupby("gene")["product"].agg(stable_join)
    parsed = [str(value).split(";") if pd.notna(value) else [] for value in info["Genes"]]
    result = pd.DataFrame(index=info.index)
    result["gene"] = info["Genes"]
    result["Is_division_genes"] = [int(any(g.strip().lower() in division_genes for g in genes)) for genes in parsed]
    dcw = selection["Is_dcw_genes"].replace({"True": "1", "False": "0"})
    if not dcw.dropna().isin(["0", "1"]).all():
        raise ValueError("Table S4 dcw flags must be 0 or 1")
    result["Is_dcw_genes"] = dcw.reindex(info.index).fillna("")
    for exp in EXPERIMENTS:
        result[f"Retained after preprocessing ({exp} upshift)"] = np.where(info.index.isin(models[exp].index), "Yes", "No")
    result["product"] = [stable_join(product_map.get(g.strip().lower(), pd.NA) for g in genes) for genes in parsed]
    result = result.join(raw["PPM"].set_axis(labels, axis=1))
    result[["Protein.Names", "Protein.Ids"]] = info[["Protein.Names", "Protein.Ids"]]
    s4 = result.reindex(selection.index)
    return raw, result, s4


def compare_matrix(new, old, atol=1e-12, rtol=1e-9):
    common = new.index.intersection(old.index)
    keys = new.columns.intersection(old.columns)
    a, b = new.loc[common, keys].to_numpy(dtype=float), old.loc[common, keys].to_numpy(dtype=float)
    finite = np.isfinite(a) & np.isfinite(b)
    mismatch = ~np.isclose(a, b, atol=atol, rtol=rtol, equal_nan=True)
    relative = np.abs(a[finite] - b[finite]) / np.maximum(np.abs(b[finite]), 1e-100)
    return {"new_rows": len(new), "reference_rows": len(old), "common_rows": len(common),
            "common_samples": len(keys), "new_only": new.index.difference(old.index).tolist(),
            "reference_only": old.index.difference(new.index).tolist(),
            "missing_mismatches": int(np.count_nonzero(np.isnan(a) != np.isnan(b))),
            "cell_mismatches": int(np.count_nonzero(mismatch)), "finite_cells": int(finite.sum()),
            "max_relative_error": float(relative.max()) if relative.size else None}


def compare_history(raw, full, s4, models, mode_dir, reference, paper):
    result = {"matrices": {}, "stage2": {}, "retention": {}}
    for metric in MEASURES:
        old = pd.read_csv(reference / f"{metric}.csv").set_index(GROUP)
        old = old.loc[old.notna().any(axis=1), raw[metric].columns]
        result["matrices"][metric] = compare_matrix(raw[metric], old)
    for exp in EXPERIMENTS:
        for metric in ("MaxLFQ_log2", "PPM"):
            new = pd.read_csv(mode_dir / f"{exp}_{metric}_filtered.csv").set_index(GROUP).drop(columns="gene")
            old = pd.read_csv(reference / f"{exp}_{metric}_filtered.csv").set_index(GROUP).drop(columns="gene")
            result["stage2"][f"{exp}_{metric}"] = compare_matrix(new[old.columns], old)
        old_model = pd.read_csv(reference / f"{exp}_formal_response_all.csv").set_index(GROUP)
        new_model = models[exp]
        result["retention"][exp] = {
            "new_count": len(new_model), "reference_count": len(old_model),
            "new_only": new_model.index.difference(old_model.index).tolist(),
            "reference_only": old_model.index.difference(new_model.index).tolist(),
            "limma_numeric_check": compare_matrix(
                new_model.select_dtypes(include="number"), old_model.select_dtypes(include="number"),
                atol=1e-9, rtol=1e-6,
            ),
            "significant_set_matches": set(new_model.index[new_model["Significant"]]) ==
                                       set(old_model.index[old_model["Significant"]]),
        }
    labels = paper["paper_label"].tolist()
    old_s4 = pd.read_csv(reference / "TableS4_Protein.csv").set_index(GROUP)
    s4_values = compare_matrix(s4[labels], old_s4[labels])
    checks = {}
    for column in s4.columns.difference(labels):
        left, right = s4[column].astype("string"), old_s4[column].astype("string")
        if column in ("Is_division_genes", "Is_dcw_genes"):
            left = left.replace({"True": "1", "False": "0", "1.0": "1", "0.0": "0"})
            right = right.replace({"True": "1", "False": "0", "1.0": "1", "0.0": "0"})
        different = ~left.fillna("").eq(right.fillna(""))
        checks[column] = left.index[different].tolist()
    result["table_s4"] = {"values": s4_values, "annotation_or_flag_mismatches": checks,
                          "row_order_matches": s4.index.equals(old_s4.index),
                          "columns_match": s4.columns.equals(old_s4.columns)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--inputs", type=Path, default=ROOT / "inputs")
    parser.add_argument("--output", type=Path, default=ROOT)
    parser.add_argument("--reference", type=Path, help="Optional historical CSV comparison directory")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--r-library", type=Path)
    parser.add_argument("--classification", choices=["historical", "batch2", "both"], default="batch2")
    args = parser.parse_args()
    root, inputs, report = args.output.resolve(), args.inputs.resolve(), args.report.resolve()
    for sub in ("intermediate", "outputs", "logs", "tmp"):
        (root / sub).mkdir(parents=True, exist_ok=True)

    def get(name):
        return pd.read_csv(inputs / f"{name}.csv", dtype=str, encoding="utf-8-sig")

    papers = get("paper_samples").assign(paper_order=lambda d: d["paper_order"].astype(int)).sort_values("paper_order")
    meta = get("batch2_sample_metadata")
    freeze = get("batch2_run_manifest")
    selection = get("table_s4_selection").set_index(GROUP)
    if not selection.index.is_unique:
        raise ValueError("Table S4 selection has duplicate protein groups")
    division = get("division_gene_annotations")
    products = get("product_annotations")
    provenance = {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__,
                  "pyarrow": pyarrow.__version__, "report_sha256": sha256(report),
                  "inputs_sha256": {path.name: sha256(path) for path in sorted(inputs.glob("*.csv"))}}
    quantity, retained = prepare_quantification(report, freeze)
    sample_design(meta, papers, retained)
    print(f"Rebuilding {len(retained)} retained Batch2 samples and selecting {len(papers)} manuscript samples.", flush=True)
    modes = ("historical", "batch2") if args.classification == "both" else (args.classification,)
    summary = {"provenance": provenance, "modes": {}}
    group_tables = {}
    for mode in modes:
        mode_dir = root / "intermediate" / mode
        mode_dir.mkdir(exist_ok=True)
        if mode == "historical":
            groups = get("historical_protein_group_classification")
        else:
            accession_classes = get("fasta_accession_sources").set_index("accession")["source_class"].to_dict()
            groups = classify_batch2(report, accession_classes)
        group_tables[mode] = groups
        write(groups, mode_dir / "protein_group_classification.csv")
        if groups["Final_Source_Class"].eq("UNKNOWN").any():
            raise ValueError(f"{mode} has unknown groups. See {mode_dir}")
        matrices = quantify(quantity, groups, retained, mode_dir)
        info = groups.set_index(GROUP)
        counts, models = {}, {}
        for exp in EXPERIMENTS:
            counts[exp] = stage2_filter(matrices, info, meta, papers, exp, mode_dir)
            models[exp] = run_r_model(args.rscript, args.r_library, mode_dir, root, inputs / "batch2_sample_metadata.csv", exp)
        raw, full, s4 = make_exports(matrices, groups, papers, selection, division, products, models)
        out_dir = root / "outputs" / mode
        out_dir.mkdir(exist_ok=True)
        for metric in MEASURES[:2]:
            write(raw[metric].set_axis(papers["paper_label"], axis=1), out_dir / f"Batch2_{metric}_full.csv", index=True)
        write(full, out_dir / "Batch2_PPM_full.csv", index=True)
        write(s4, out_dir / "TableS4_Protein_reproduced.csv", index=True)
        metrics = {"batch2_protein_groups": len(matrices["PPM"]), "paper_protein_groups": len(full),
                   "paper_samples": len(papers), "stage2_counts": counts,
                   "retained_models": {exp: len(models[exp]) for exp in EXPERIMENTS}}
        if args.reference:
            metrics["historical_comparison"] = compare_history(raw, full, s4, models, mode_dir, args.reference.resolve(), papers)
        summary["modes"][mode] = metrics
        print(mode, {k: v for k, v in metrics.items() if k != "historical_comparison"}, flush=True)
    if len(group_tables) == 2:
        a = group_tables["historical"].set_index(GROUP)
        b = group_tables["batch2"].set_index(GROUP)
        common = a.index.intersection(b.index)
        diffs = a.loc[common, "Final_Source_Class"].compare(b.loc[common, "Final_Source_Class"])
        summary["classification_comparison"] = {
            "shared_protein_groups": len(common), "historical_only_groups": len(a.index.difference(b.index)),
            "batch2_only_groups": len(b.index.difference(a.index)),
            "changed_source_class": diffs.index.tolist(),
        }
        full_a = pd.read_csv(root / "outputs" / "historical" / "Batch2_PPM_full.csv").set_index(GROUP)
        full_b = pd.read_csv(root / "outputs" / "batch2" / "Batch2_PPM_full.csv").set_index(GROUP)
        if full_a.index.equals(full_b.index):
            differences = full_a.compare(full_b, result_names=("historical", "batch2"))
            write(differences, root / "intermediate" / "classification_output_differences.csv", index=True)
            summary["classification_comparison"]["full_ppm_table_difference_columns"] = list(
                differences.columns.get_level_values(0).unique()
            )
            summary["classification_comparison"]["full_ppm_table_difference_groups"] = differences.index.tolist()
        write(diffs, root / "intermediate" / "changed_source_classes.csv", index=True)
    (root / "outputs" / "validation.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
