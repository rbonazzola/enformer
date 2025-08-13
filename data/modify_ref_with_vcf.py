import os
import gzip
import subprocess
import urllib.request
import tempfile

def download_geuvadis_sample_list(output_path: str = "datasets/E-GEUV-1.sdrf.txt") -> str:
    """
    Downloads the GEUVADIS SDRF file from ArrayExpress and returns its path.
    Skips download if the file already exists.
    """
    url = "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt"

    if os.path.exists(output_path):
        print(f"📎 GEUVADIS SDRF already exists at {output_path}, skipping download.")
    else:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        print(f"⬇️ Downloading GEUVADIS SDRF from {url}")
        urllib.request.urlretrieve(url, output_path)
        print(f"💾 Saved to {output_path}")

    return output_path



def extract_geuvadis_ids(sdrf_path: str, european_only=True) -> list[str]:
    """
    Extracts GEUVADIS sample IDs from the SDRF file and returns them as a list.
    """

    # european_pops = {"CEU", "FIN", "GBR", "TSI"}
    european_pops = {"Utah", "Finnish", "British", "Tuscan"}

    print(f"📋 Extracting {'European ' if european_only else ''}GEUVADIS sample IDs from {sdrf_path}")
    sample_ids = set()

    with open(sdrf_path) as f:
        header = f.readline().strip().split("\t")
        pop_idx = header.index("Characteristics[ancestry category]")
        sample_idx = header.index("Source Name")

        for line in f:
            fields = line.strip().split("\t")
            population = fields[pop_idx].replace("population: ", "").strip()
            sample = fields[sample_idx].strip()
            if population in european_pops:
                sample_ids.add(sample)

    print(f"✅ Extracted {len(sample_ids)} unique GEUVADIS sample IDs")
    return sorted(sample_ids)


def extract_geuvadis_european_ids(sdrf_path: str) -> list[str]:
    """
    Extracts European GEUVADIS sample IDs (CEU, FIN, GBR, TSI) from the SDRF file.
    """
    print(f"📋 Extracting European GEUVADIS sample IDs from {sdrf_path}")
    european_pops = {"CEU", "FIN", "GBR", "TSI"}
    sample_ids = set()

    with open(sdrf_path) as f:
        header = f.readline().strip().split("\t")
        pop_idx = header.index("Characteristics[population]")
        sample_idx = header.index("Source Name")

        for line in f:
            fields = line.strip().split("\t")
            population = fields[pop_idx].replace("population: ", "").strip()
            sample = fields[sample_idx].strip()
            if population in european_pops:
                sample_ids.add(sample)

    print(f"✅ Found {len(sample_ids)} European GEUVADIS sample IDs")
    return sorted(sample_ids)



def filter_vcf_by_samples(input_vcf: str, output_vcf: str, sample_ids: list[str]) -> None:
    """
    Uses bcftools to subset a VCF file based on a list of sample IDs, passed via stdin.
    """
    print(f"🧬 Filtering {os.path.basename(input_vcf)} → {os.path.basename(output_vcf)}")
    proc = subprocess.Popen([
        "bcftools", "view",
        "-S", "-",  # read sample IDs from stdin
        "-Oz", "-o", output_vcf,
        input_vcf
    ], stdin=subprocess.PIPE, text=True)
    
    proc.communicate("\n".join(sample_ids))
    if proc.returncode != 0:
        raise RuntimeError(f"bcftools view failed on {input_vcf}")

    subprocess.run(["bcftools", "index", output_vcf], check=True)


def summarize_vcf(input_vcf: str, geuvadis_sample_ids: list[str]) -> None:
    """
    Summarizes sample and variant statistics from a VCF file, including overlap with GEUVADIS sample IDs.
    """    

    # Get all sample IDs from the VCF header
    with gzip.open(input_vcf, 'rt') as f:
        for line in f:
            if line.startswith("#CHROM"):
                vcf_sample_ids = line.strip().split()[9:]
                break

    # Compute overlaps
    matched_ids = [sid for sid in geuvadis_sample_ids if sid in vcf_sample_ids]

    print("👥 Sample statistics:")
    print(f"   • GEUVADIS sample IDs       : {len(geuvadis_sample_ids)}")
    print(f"   • Samples in original VCF   : {len(vcf_sample_ids)}")
    print(f"   • GEUVADIS samples in VCF   : {len(matched_ids)}")

    # Run bcftools stats
    stats_file = input_vcf + ".stats.txt"
    with open(stats_file, "w") as fout:
        subprocess.run(["bcftools", "stats", "-s", "-", input_vcf], stdout=fout, check=True)

    # Parse bcftools stats
    n_variants = 0
    n_snps = 0
    n_indels = 0
    n_biallelic_snps = 0

    with open(stats_file) as f:
        for line in f:
            if line.startswith("SN\t0\tnumber of records:"):
                n_variants = int(line.strip().split("\t")[-1])
            elif line.startswith("SN\t0\tnumber of SNPs:"):
                n_snps = int(line.strip().split("\t")[-1])
            elif line.startswith("SN\t0\tnumber of indels:"):
                n_indels = int(line.strip().split("\t")[-1])
            elif line.startswith("SN\t0\tnumber of multiallelic SNP sites:"):
                multi = int(line.strip().split("\t")[-1])
                n_biallelic_snps = n_snps - multi

    print("📊 Variant statistics:")
    print(f"   • Total variants      : {n_variants:,}")
    print(f"   • SNPs                : {n_snps:,}")
    print(f"   • Indels              : {n_indels:,}")
    print(f"   • Biallelic SNPs only : {n_biallelic_snps:,}")


def get_sample_ids_from_vcf(vcf_path: str) -> set[str]:

    """
    Extracts the list of sample IDs from the header of a .vcf.gz file.
    """
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith("#CHROM"):
                fields = line.strip().split()
                return set(fields[9:])  # samples start at column 10
    return set()


def process_all_vcfs(input_dir: str, output_dir: str, sample_ids: list[str]) -> None:

    """
    Filters all VCF files in the input directory using the GEUVADIS sample IDs,
    and writes the filtered VCFs to the output directory.
    """
    os.makedirs(output_dir, exist_ok=True)
    vcf_files = sorted(f for f in os.listdir(input_dir) if f.endswith(".vcf.gz") and f.startswith("ALL.chr"))

    for fname in vcf_files:
        input_vcf = os.path.join(input_dir, fname)
        chrom = fname.split(".")[1]
        output_vcf = os.path.join(output_dir, f"{chrom}_GEUVADIS.vcf.gz")
        vcf_samples = get_sample_ids_from_vcf(input_vcf)
        common_ids = [sid for sid in sample_ids if sid in vcf_samples]

        if not common_ids:
            print(f"⚠️ Skipping {fname} — no overlapping GEUVADIS samples")
            continue

        summarize_vcf(input_vcf, common_ids)
        filter_vcf_by_samples(input_vcf, output_vcf, common_ids)

    print("🏁 All VCFs have been filtered by GEUVADIS samples.")


def main():
    
    vcf_dir = "datasets/vcfs"  # Folder with all the original 1000G VCFs
    out_dir = "transforms/vcfs_geuvadis"  # Output folder

    # Step 1: Download GEUVADIS SDRF and extract sample IDs
    sdrf_path = download_geuvadis_sample_list()
    sample_ids = extract_geuvadis_ids(sdrf_path)

    # Step 2: Process all VCFs using the sample IDs
    process_all_vcfs(vcf_dir, out_dir, sample_ids)


if __name__ == "__main__":
    main()