import json
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset


def extract_vcfs(ds):
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    # Keep VCFs and their TBI indices
    vcf_files = df[df["file"].str.endswith(".vcf.gz")]["file"].tolist()
    tbi_files = df[df["file"].str.endswith(".vcf.gz.tbi")]["file"].tolist()

    if not vcf_files:
        raise ValueError("No VCF files found in dataset")

    # Match VCFs to their TBIs by basename
    pairs = []
    for vcf in sorted(vcf_files):
        tbi = vcf + ".tbi"
        if tbi not in tbi_files:
            raise ValueError(f"Missing TBI index for {vcf}")
        pairs.append({"vcf": vcf, "tbi": tbi})

    return pairs


def main():
    ds = PreprocessDataset.from_running()

    print("=== ds.files preview ===")
    print(ds.files.to_string(index=False))
    print(f"\n=== shape: {ds.files.shape} ===")
    print(f"\n=== unique samples: {ds.files['sample'].unique().tolist()} ===")
    print(f"\n=== unique datasets: {ds.files['dataset'].unique().tolist()} ===")
    # vcf_pairs = extract_vcfs(ds)
    # n_callers = len(vcf_pairs)

    # # Fail fast — consensus is meaningless with only one caller
    # if n_callers < 2:
    #     raise ValueError(
    #         f"Consensus calling requires at least 2 VCFs, only {n_callers} found. "
    #         "Please provide VCFs from multiple callers."
    #     )
    # print(f"Found {n_callers} VCFs")

    # # Store each VCF/TBI with a predictable index-based name
    # for i, pair in enumerate(vcf_pairs):
    #     ds.add_param(f"caller_vcf_{i}", pair["vcf"])
    #     ds.add_param(f"caller_tbi_{i}", pair["tbi"])

    # # Store n so Nextflow knows how many callers there are
    # ds.add_param("n_callers", n_callers)

    # # n-1 threshold: for SNVs accept if called by at least n-1 callers
    # # for indels accept if called by both (min 2)
    # ds.add_param("snv_min_callers", max(1, n_callers - 1))
    # ds.add_param("indel_min_callers", min(2, n_callers))a

    # print("\nFinal parameters:")
    # print(json.dumps(ds.params, indent=2, default=str))

if __name__ == "__main__":
    main()