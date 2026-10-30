import json
import os
import re
import subprocess
from cirro.helpers.preprocess_dataset import PreprocessDataset
# ── Expected filenames ────────────────────────────────────────────────────────
# Callers are identified by the caller name in the filename:
#
#   MuTect       <sample>.mutect[1|2].<anything>.vcf[.gz]
#   DeepSomatic  <sample>.deepsomatic.vcf.gz
#   Strelka      <tumor>_vs_<normal>.strelka.somatic_snvs.vcf.gz
#                <tumor>_vs_<normal>.strelka.somatic_indels.vcf.gz
#                (nf-core/sarek writes SNVs and indels separately; both are needed)
#
# <sample> is the tumor sample ID; it is the key used to group the callers.
# For Strelka it is everything before _vs_.

# (caller, file kind) -> filename pattern. The first group is the sample ID.
FILE_PATTERNS = [
    ("mutect",      "vcf",    re.compile(r"^(.+?)\.mutect[12]?\..*vcf(\.gz)?$", re.IGNORECASE)),
    ("deepsomatic", "vcf",    re.compile(r"^(.+)\.deepsomatic\.vcf\.gz$", re.IGNORECASE)),
    ("strelka",     "snvs",   re.compile(r"^(.+?)_vs_.+\.strelka\.somatic_snvs\.vcf\.gz$", re.IGNORECASE)),
    ("strelka",     "indels", re.compile(r"^(.+?)_vs_.+\.strelka\.somatic_indels\.vcf\.gz$", re.IGNORECASE)),
]

# Callers that can be chosen in the UI, and the ones that only call SNVs.
ALL_CALLERS = ["mutect", "deepsomatic", "strelka"]
SNV_ONLY_CALLERS = ["mutect"]

# Parts each caller must have.
REQUIRED_KINDS = {
    "mutect":      {"vcf"},
    "deepsomatic": {"vcf"},
    "strelka":     {"snvs", "indels"},
}

# Set to False once input discovery is working, to keep logs short.
DEBUG_INPUTS = True


def _dump(title, text):
    print(f"\n===== {title} =====\n{text}", flush=True)


def _find_s3_paths(obj):
    """Yield every s3:// string found anywhere in a nested dict/list."""
    if isinstance(obj, str) and obj.startswith("s3://"):
        yield obj
    elif isinstance(obj, dict):
        for v in obj.values():
            yield from _find_s3_paths(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _find_s3_paths(v)


def debug_inputs(ds, max_lines=200):
    """Print everything the preprocess step can see, to diagnose missing inputs."""
    _dump("PARAMS", json.dumps(ds.params, indent=2, default=str))
    _dump("METADATA", json.dumps(ds.metadata, indent=2, default=str))
    samplesheet = getattr(ds, "samplesheet", None)
    _dump("SAMPLESHEET", samplesheet.to_string() if samplesheet is not None else "None")
    _dump("FILES", ds.files.to_string() if not ds.files.empty else "(empty)")
    _dump("LOCAL WORKING DIR", "\n".join(sorted(os.listdir("."))))

    for path in sorted(set(_find_s3_paths(ds.metadata))):
        try:
            result = subprocess.run(
                ["aws", "s3", "ls", "--recursive", path],
                capture_output=True, text=True,
            )
            lines = (result.stdout or result.stderr or "(empty)").splitlines()
        except FileNotFoundError:
            lines = ["aws CLI not available in this container"]
        if len(lines) > max_lines:
            lines = lines[:max_lines] + [f"... ({len(lines) - max_lines} more lines truncated)"]
        _dump(f"S3 LISTING: {path}", "\n".join(lines))


def as_bool(value, default=True):
    """Cirro may hand booleans over as real bools or as 'true'/'false' strings."""
    if value is None or value == "":
        return default
    if isinstance(value, str):
        return value.strip().lower() in ("true", "1", "yes", "on")
    return bool(value)


def selected_callers(params) -> list:
    """Callers ticked in the UI (use_<caller> checkboxes), in canonical order."""
    return [c for c in ALL_CALLERS if as_bool(params.get(f"use_{c}"))]


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


def complete_callers(callers: dict, wanted: list) -> list:
    """Wanted callers that have every part they need (Strelka needs both SNVs and indels)."""
    return [c for c in wanted if REQUIRED_KINDS[c] <= set(callers.get(c, {}))]


def validate_groups(groups: dict, wanted: list) -> list:
    """
    Keep only samples where every selected caller is complete. Returns a list of
    sample dicts for the workflow.
    """
    keys = {
        ("mutect", "vcf"): "vcf_mutect",
        ("deepsomatic", "vcf"): "vcf_deepsomatic",
        ("strelka", "snvs"): "vcf_strelka_snvs",
        ("strelka", "indels"): "vcf_strelka_indels",
    }
    valid = []
    for sample_id, callers in sorted(groups.items()):
        have = complete_callers(callers, wanted)
        missing = [c for c in wanted if c not in have]
        if missing:
            print(f"[WARN] Skipping sample '{sample_id}' - missing/incomplete: {', '.join(missing)}")
            continue
        entry = {"sample_id": sample_id, "callers": have}
        for c in have:
            for kind, path in callers[c].items():
                entry[keys[(c, kind)]] = path
        valid.append(entry)
        print(f"[OK] Sample '{sample_id}' - " + ", ".join(
            f"{c}={'/'.join(callers[c][k].split('/')[-1] for k in sorted(callers[c]))}"
            for c in have))
    return valid


def main():
    ds = PreprocessDataset.from_running()

    wanted = selected_callers(ds.params)
    if not wanted:
        raise ValueError("No callers selected; tick at least one of MuTect, DeepSomatic, Strelka")
    indel_callers = [c for c in wanted if c not in SNV_ONLY_CALLERS]
    print(f"[INFO] Callers selected: {', '.join(wanted)}")
    print(f"[INFO] Indel callers: {', '.join(indel_callers) or '<none - SNV-only callers selected>'}")

    if DEBUG_INPUTS:
        debug_inputs(ds)

    if ds.files.empty:
        raise ValueError(
            "No files found in the input dataset(s). Check which datasets were "
            "selected when launching the run and whether their files are registered in Cirro."
        )

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

    valid_samples = validate_groups(groups, wanted)
    print(f"[INFO] Samples to process: {len(valid_samples)}")

    if not valid_samples:
        raise ValueError(
            f"No sample has every selected caller ({', '.join(wanted)}). Expected files such as: "
            "<sample>.mutect*.vcf[.gz], <sample>.deepsomatic.vcf.gz, "
            "and <tumor>_vs_<normal>.strelka.somatic_snvs.vcf.gz / .somatic_indels.vcf.gz"
        )

    # Emit to Nextflow params
    ds.add_param("samples", json.dumps(valid_samples))
    ds.add_param("n_samples", len(valid_samples))
    ds.add_param("callers", ",".join(wanted))

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