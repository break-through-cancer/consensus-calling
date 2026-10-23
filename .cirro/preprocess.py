import json
import re
from cirro.helpers.preprocess_dataset import PreprocessDataset


# ── Expected filenames ────────────────────────────────────────────────────────
# Callers are identified by the caller name in the filename. Each caller
# pipeline publishes its VCF under a caller-named file:
#
#   MuTect1      <sample>.mutect1.pass.vcf               (plain, not bgzipped)
#   DeepSomatic  <sample>.deepsomatic.vcf.gz
#   Strelka      <tumor>_vs_<normal>.strelka.somatic_snvs.vcf.gz
#                <tumor>_vs_<normal>.strelka.somatic_indels.vcf.gz
#                (nf-core/sarek writes SNVs and indels separately; both are needed)
#   VarScan      <sample>.varscan.vcf[.gz]               (assumed; adjust if different)
#
# <sample> is the tumor sample ID; it is the key used to group the callers.
# For Strelka it is everything before _vs_.

# (caller, file kind) -> filename pattern. The first group is the sample ID.
# "vcf" is a complete VCF for the caller; Strelka contributes two parts that the
# workflow concatenates.
FILE_PATTERNS = [
    ("mutect1",     "vcf",    re.compile(r"^(.+)\.mutect1\.pass\.vcf$", re.IGNORECASE)),
    ("deepsomatic", "vcf",    re.compile(r"^(.+)\.deepsomatic\.vcf\.gz$", re.IGNORECASE)),
    ("strelka",     "snvs",   re.compile(r"^(.+?)_vs_.+\.strelka\.somatic_snvs\.vcf\.gz$", re.IGNORECASE)),
    ("strelka",     "indels", re.compile(r"^(.+?)_vs_.+\.strelka\.somatic_indels\.vcf\.gz$", re.IGNORECASE)),
    ("varscan",     "vcf",    re.compile(r"^(.+)\.varscan\.vcf(\.gz)?$", re.IGNORECASE)),
]

# Every sample must have all of these to be used.
REQUIRED_CALLERS = ["mutect1", "deepsomatic", "strelka", "varscan"]

# Parts each caller must have.
REQUIRED_KINDS = {
    "mutect1":     {"vcf"},
    "deepsomatic": {"vcf"},
    "strelka":     {"snvs", "indels"},
    "varscan":     {"vcf"},
}

# Callers whose indel calls count towards the indel consensus.
INDEL_CALLERS = ["mutect1", "deepsomatic", "strelka", "varscan"]


def classify_file(filename: str):
    """
    Return (sample_id, caller, kind) if the filename matches a known caller
    pattern, otherwise (None, None, None).
    """
    basename = filename.split("/")[-1]
    for caller, kind, pattern in FILE_PATTERNS:
        m = pattern.match(basename)
        if m:
            return m.group(1), caller, kind
    return None, None, None


def discover_sample_groups(all_files: list) -> dict:
    """
    Group files into {sample_id: {caller: {kind: path}}}.
    """
    groups = {}
    for f in all_files:
        sample_id, caller, kind = classify_file(f)
        if sample_id is None:
            continue
        parts = groups.setdefault(sample_id, {}).setdefault(caller, {})
        if kind in parts:
            print(
                f"[WARN] Duplicate {caller} {kind} file for sample {sample_id} — "
                f"keeping first, ignoring {f}"
            )
        else:
            parts[kind] = f
    return groups


def missing_for_sample(callers: dict) -> list:
    """Return the list of missing '<caller>' / '<caller>:<kind>' items."""
    missing = []
    for caller in REQUIRED_CALLERS:
        have = set(callers.get(caller, {}))
        need = REQUIRED_KINDS[caller]
        if not have:
            missing.append(caller)
        else:
            missing.extend(f"{caller}:{kind}" for kind in sorted(need - have))
    return missing


def validate_groups(groups: dict) -> list:
    """
    Keep samples that have every required caller (and every part of it).
    Warn and skip incomplete samples rather than crashing.
    Returns a list of sample dicts for the workflow.
    """
    valid = []
    for sample_id, callers in sorted(groups.items()):
        missing = missing_for_sample(callers)
        if missing:
            found = ", ".join(sorted(callers)) or "none"
            print(
                f"[WARN] Skipping sample '{sample_id}' — "
                f"missing: {', '.join(missing)} (found callers: {found})"
            )
            continue
        valid.append({
            "sample_id": sample_id,
            "vcf_mutect1":         callers["mutect1"]["vcf"],
            "vcf_deepsomatic":     callers["deepsomatic"]["vcf"],
            "vcf_strelka_snvs":    callers["strelka"]["snvs"],
            "vcf_strelka_indels":  callers["strelka"]["indels"],
            "vcf_varscan":         callers["varscan"]["vcf"],
        })
        print(
            f"[OK] Sample '{sample_id}' — "
            + ", ".join(
                f"{c}={'/'.join(p.split('/')[-1] for p in (callers[c][k] for k in sorted(callers[c])))}"
                for c in REQUIRED_CALLERS
            )
        )
    return valid


def main():
    ds = PreprocessDataset.from_running()

    print("[DEBUG] ds.files columns:", list(ds.files.columns))
    print(ds.files.head().to_string())

    # Pick the column that contains paths / filenames
    file_col = None
    for col in ["s3_path", "file", "path", "uri", "name", "filename"]:
        if col in ds.files.columns:
            file_col = col
            break

    if file_col is None:
        raise ValueError(f"Could not find file path column in ds.files. Columns: {list(ds.files.columns)}")

    all_files = ds.files[file_col].astype(str).tolist()

    print(f"[INFO] Total files in dataset: {len(all_files)}")
    print(f"[INFO] VCF files found: {sum(1 for f in all_files if f.endswith(('.vcf', '.vcf.gz')))}")

    # Discover and validate sample groups
    groups = discover_sample_groups(all_files)
    print(f"\n[INFO] Sample groups detected: {len(groups)}")

    valid_samples = validate_groups(groups)
    print(f"[INFO] Valid samples (all {len(REQUIRED_CALLERS)} callers present): {len(valid_samples)}")

    if not valid_samples:
        raise ValueError(
            "No valid sample groups found. Every sample needs: "
            "<sample>.mutect1.pass.vcf, <sample>.deepsomatic.vcf.gz, "
            "<tumor>_vs_<normal>.strelka.somatic_snvs.vcf.gz and .somatic_indels.vcf.gz, "
            "and <sample>.varscan.vcf[.gz]"
        )

    # Emit to Nextflow params
    ds.add_param("samples", json.dumps(valid_samples))
    ds.add_param("n_samples", len(valid_samples))
    ds.add_param("n_callers", len(REQUIRED_CALLERS))
    ds.add_param("indel_callers", ",".join(INDEL_CALLERS))
    ds.add_param("indel_min_callers", 2)

    print("\n[INFO] Sample manifest:")
    print(json.dumps(valid_samples, indent=2))


if __name__ == "__main__":
    main()

# import json
# from cirro.helpers.preprocess_dataset import PreprocessDataset


# def validate_vcf(vcf_path):
#     if not vcf_path:
#         raise ValueError("Empty VCF path provided")

#     if not str(vcf_path).endswith(".vcf.gz"):
#         raise ValueError(f"Specified file is not a .vcf.gz file: {vcf_path}")

#     tbi_path = f"{vcf_path}.tbi"
#     return {"vcf": vcf_path, "tbi": tbi_path}


# def main():
#     ds = PreprocessDataset.from_running()

#     user_vcfs = [
#         ds.params.get("vcf1"),
#         ds.params.get("vcf2"),
#         ds.params.get("vcf3"),
#         ds.params.get("vcf4"),
#         ds.params.get("vcf5"),
#         ds.params.get("vcf6"),
#         ds.params.get("vcf7"),
#     ]

#     user_vcfs = [v for v in user_vcfs if v not in (None, "", "null")]

#     if len(user_vcfs) < 3:
#         raise ValueError(
#             f"Consensus calling requires at least 3 VCFs, but only {len(user_vcfs)} were specified."
#         )

#     if len(user_vcfs) > 7:
#         raise ValueError(
#             f"This workflow supports at most 7 VCFs, but {len(user_vcfs)} were specified."
#         )

#     vcf_pairs = [validate_vcf(vcf) for vcf in user_vcfs]

#     print(f"Using {len(vcf_pairs)} user-specified VCFs")

#     ds.add_param("n_callers", len(vcf_pairs))
#     ds.add_param("snv_min_callers", len(vcf_pairs) - 1)
#     ds.add_param("indel_min_callers", 2)

#     print("\nSelected VCF pairs:")
#     print(json.dumps(vcf_pairs, indent=2))

#     print("\nFinal parameters:")
#     print(json.dumps(ds.params, indent=2, default=str))


# if __name__ == "__main__":
#     main()