import json
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset


def extract_vcfs(ds):
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    vcf_files = df[df["file"].str.endswith(".vcf.gz")]["file"].tolist()
    tbi_files = df[df["file"].str.endswith(".vcf.gz.tbi")]["file"].tolist()

    if not vcf_files:
        raise ValueError("No VCF files found in dataset")

    pairs = []
    for vcf in sorted(vcf_files):
        tbi = vcf + ".tbi"
        if tbi not in tbi_files:
            raise ValueError(f"Missing TBI index for {vcf}")
        pairs.append({"vcf": vcf, "tbi": tbi})

    return pairs


def main():
    ds = PreprocessDataset.from_running()

    vcf_pairs = extract_vcfs(ds)
    n_callers = len(vcf_pairs)

    if n_callers < 3:
        raise ValueError(
            f"Consensus calling requires at least 3 VCFs, but only {n_callers} were found."
        )

    if n_callers > 7:
        raise ValueError(
            f"This workflow supports at most 7 VCFs, but {n_callers} were found."
        )

    print(f"Found {n_callers} VCFs")

    for i, pair in enumerate(vcf_pairs):
        ds.add_param(f"caller_vcf_{i}", pair["vcf"])
        ds.add_param(f"caller_tbi_{i}", pair["tbi"])

    ds.add_param("n_callers", n_callers)
    ds.add_param("snv_min_callers", n_callers - 1)
    ds.add_param("indel_min_callers", 2)

    print("\nFinal parameters:")
    print(json.dumps(ds.params, indent=2, default=str))


if __name__ == "__main__":
    main()