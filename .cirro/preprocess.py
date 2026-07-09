import json
import re
from cirro.helpers.preprocess_dataset import PreprocessDataset


# ── Caller suffix patterns ────────────────────────────────────────────────────
# Strelka and Freebayes: ID_vs_normalID.strelka.vcf.gz
#                        ID_vs_normalID.freebayes.filtered.vcf.gz
# Mutect2:               ID.merged.vcf.gz
#
# For paired callers we capture everything before _vs_ as the tumor sample ID.
# For Mutect2 we capture everything before .merged.
# All three yield the same tumor ID so grouping works cleanly.

CALLER_PATTERNS = [
    # More specific first: freebayes.filtered before a hypothetical freebayes-only
    ("freebayes", re.compile(r"^(.+?)_vs_.+\.freebayes\.filtered\.vcf\.gz$")),
    ("strelka",   re.compile(r"^(.+?)_vs_.+\.strelka\.vcf\.gz$")),
    ("mutect2",   re.compile(r"^(.+)\.merged\.vcf\.gz$")),
]

REQUIRED_CALLERS = {"freebayes", "strelka", "mutect2"}


def classify_file(filename: str):
    """
    Return (sample_id, caller_name) if the filename matches a known caller
    pattern, otherwise return (None, None).

    Sample ID is always the tumor ID portion:
      - Paired callers (Strelka, Freebayes): everything before _vs_
      - Mutect2: everything before .merged
    """
    basename = filename.split("/")[-1]
    for caller_name, pattern in CALLER_PATTERNS:
        m = pattern.match(basename)
        if m:
            return m.group(1), caller_name
    return None, None


def discover_sample_groups(all_files: list) -> dict:
    """
    Scan all files in the dataset, extract sample ID + caller for each
    matching VCF, and group into {sample_id: {caller: path}}.
    """
    groups = {}
    for f in all_files:
        sample_id, caller = classify_file(f)
        if sample_id is None:
            continue
        if sample_id not in groups:
            groups[sample_id] = {}
        if caller in groups[sample_id]:
            print(
                f"[WARN] Duplicate {caller} VCF for sample {sample_id} — "
                f"keeping first, ignoring {f}"
            )
        else:
            groups[sample_id][caller] = f
    return groups


def validate_groups(groups: dict) -> list:
    """
    Filter to samples that have all three required callers.
    Warn and skip incomplete samples rather than crashing.
    Returns a list of validated sample dicts.
    """
    valid = []
    for sample_id, callers in sorted(groups.items()):
        missing = REQUIRED_CALLERS - set(callers.keys())
        if missing:
            print(
                f"[WARN] Skipping sample '{sample_id}' — "
                f"missing callers: {', '.join(sorted(missing))}"
            )
            continue
        valid.append({
            "sample_id": sample_id,
            "vcf_mutect2":   callers["mutect2"],
            "vcf_strelka":   callers["strelka"],
            "vcf_freebayes": callers["freebayes"],
        })
        print(
            f"[OK] Sample '{sample_id}' — "
            f"mutect2={callers['mutect2'].split('/')[-1]}, "
            f"strelka={callers['strelka'].split('/')[-1]}, "
            f"freebayes={callers['freebayes'].split('/')[-1]}"
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
    print(f"[INFO] VCF files found: {sum(1 for f in all_files if f.endswith('.vcf.gz'))}")

    # Discover and validate sample groups
    groups = discover_sample_groups(all_files)
    print(f"\n[INFO] Sample groups detected: {len(groups)}")

    valid_samples = validate_groups(groups)
    print(f"[INFO] Valid samples (all 3 callers present): {len(valid_samples)}")

    if not valid_samples:
        raise ValueError(
            "No valid sample groups found. Every sample needs "
            ".merged.vcf.gz, .strelka.vcf.gz, and .freebayes.filtered.vcf.gz"
        )

    # Emit to Nextflow params
    ds.add_param("samples", json.dumps(valid_samples))
    ds.add_param("n_samples", len(valid_samples))
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