import json
from cirro.helpers.preprocess_dataset import PreprocessDataset


def validate_vcf(vcf_path):
    if not vcf_path:
        raise ValueError("Empty VCF path provided")

    if not str(vcf_path).endswith(".vcf.gz"):
        raise ValueError(f"Specified file is not a .vcf.gz file: {vcf_path}")

    tbi_path = f"{vcf_path}.tbi"
    return {"vcf": vcf_path, "tbi": tbi_path}


def main():
    ds = PreprocessDataset.from_running()

    user_vcfs = [
        ds.params.get("vcf1"),
        ds.params.get("vcf2"),
        ds.params.get("vcf3"),
        ds.params.get("vcf4"),
        ds.params.get("vcf5"),
        ds.params.get("vcf6"),
        ds.params.get("vcf7"),
    ]

    user_vcfs = [v for v in user_vcfs if v not in (None, "", "null")]

    if len(user_vcfs) < 3:
        raise ValueError(
            f"Consensus calling requires at least 3 VCFs, but only {len(user_vcfs)} were specified."
        )

    if len(user_vcfs) > 7:
        raise ValueError(
            f"This workflow supports at most 7 VCFs, but {len(user_vcfs)} were specified."
        )

    vcf_pairs = [validate_vcf(vcf) for vcf in user_vcfs]

    print(f"Using {len(vcf_pairs)} user-specified VCFs")

    ds.add_param("n_callers", len(vcf_pairs))
    ds.add_param("snv_min_callers", len(vcf_pairs) - 1)
    ds.add_param("indel_min_callers", 2)

    print("\nSelected VCF pairs:")
    print(json.dumps(vcf_pairs, indent=2))

    print("\nFinal parameters:")
    print(json.dumps(ds.params, indent=2, default=str))


if __name__ == "__main__":
    main()