# Oxford-Nanopore-Mapping-and-Variant-Calling-Pipeline-
Modular pipeline for variant calling from Oxford Nanopore long-read sequencing data, with support for diploids and ployploids
______________________________________________________________________________________________________________________________

## Overview

This pipeline provides separate shell scripts for each analysis step, designed for integration with YAD GUI interfaces or command-line use. It supports both standard diploid analysis and advanced polyploid genotyping.

## Quick Start

### For Diploid Samples (2n) using 16 threads:
```bash
# Standard workflow
./01_qc.sh reads.fastq.gz qc_output 16
./02_align.sh reads.fastq.gz reference.fasta aligned.bam 16
./03_process_bam.sh aligned.bam 16
./04a_call_variants_clair3_updated.sh aligned.sorted.bam reference.fasta output model 16
./05_filter_variants_updated.sh output/merge_output.vcf.gz filtered.vcf.gz 20 5
./06_variant_stats_v2.sh filtered.vcf.gz reference.fasta stats
```

### For Polyploid Samples (3n, 4n, etc.):
```bash
# Use integrated Clair3 + FreeBayes pipeline
./01_qc.sh reads.fastq.gz qc_output 16
./02_align.sh reads.fastq.gz reference.fasta aligned.bam 16
./03_process_bam.sh aligned.bam 16
./04a_call_variants_clair3_triploid_v2.sh aligned.sorted.bam reference.fasta output model 16 3
./06_variant_stats_v2.sh output/03_final_genotypes/final_ploidy3.vcf.gz reference.fasta stats
```

## Pipeline Components

### Core Scripts (Required)

| Script | Function | Input | Output |
|--------|----------|-------|--------|
| **01_qc.sh** | Quality control (NanoPlot) | FASTQ | HTML reports, statistics |
| **02_align.sh** | Alignment (minimap2) | FASTQ + Reference | BAM file |
| **03_process_bam.sh** | Sort, index BAM | BAM | Sorted BAM + index |
| **04a_call_variants_clair3_updated.sh** | Variant calling - Diploid | Sorted BAM | VCF (diploid genotypes) |
| **04a_call_variants_clair3_triploid_v2.sh** | Variant calling - Any ploidy | Sorted BAM | VCF (correct ploidy genotypes) |
| **04c_call_variants_longshot_updated.sh** | Fast SNV calling | Sorted BAM | VCF (SNVs only) |
| **05_filter_variants_updated.sh** | Filter by quality/depth | VCF | Filtered VCF |
| **06_variant_stats_v2.sh** | Comprehensive statistics | VCF | Statistics files |

### Utility Scripts (Optional)

| Script | Function | When to Use |
|--------|----------|-------------|
| **compare_vcfs.sh** | Compare 2 VCFs | Evaluate caller concordance |
| **compare_three_vcfs.sh** | Compare 3 VCFs | Compare Clair3, Longshot, Pepper |
| **estimate_triploid_genotypes.sh** | Estimate triploid GTs from AF | Quick triploid analysis without FreeBayes |

### Example Workflows (Optional)

| File | Purpose | Users |
|------|---------|-------|
| **config.example.sh** | Parameter template | Command-line users, testing |
| **run_complete_pipeline.sh** | Automated workflow | Batch processing, HPC |

**Note for GUI Users:** These are optional - the YAD interface will handle parameter collection and script execution.

## System Requirements

- **OS**: Ubuntu 20.04 or later (tested on Ubuntu 20.04/22.04)
- **CPU**: Minimum 16 cores (more cores = faster)
- **RAM**: Minimum 256 GB (some steps need substantial memory)
- **Storage**: ~500 GB free space per sample (depends on coverage)
- **Conda**: Miniforge or Miniconda/Anaconda

## Installation

### 1. Install Miniforge (Recommended)

```bash
# Download Miniforge
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# Install
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3

# Initialize
~/miniforge3/bin/conda init bash
source ~/.bashrc

# Configure
conda config --set auto_activate_base false
conda config --set channel_priority strict
```

### 2. Create Pipeline Environments

```bash
# Main pipeline environment
mamba create -n nanopore_pipeline -c bioconda -c conda-forge \
    minimap2 samtools bcftools nanoplot \
    longshot rtg-tools bedtools freebayes \
    python=3.9 -y

# Clair3 environment
mamba create -n clair3 -c bioconda clair3 -y

# Optional: Pepper environment (if using Pepper caller)
mamba create -n pepper python=3.8 -y
conda activate pepper
pip install pepper-deepvariant
conda deactivate
```

### 3. Make Scripts Executable

```bash
cd nanopore_pipeline
chmod +x *.sh
```

## Detailed Usage

### Script 1: Quality Control

```bash
./01_qc.sh <input_fastq> <output_dir> <threads>
```

**Parameters:**
- `input_fastq`: FASTQ or FASTQ.GZ file
- `output_dir`: Directory for QC reports
- `threads`: Number of CPU threads (recommend: 8-16)

**Output:**
- `NanoPlot-report.html` - Interactive quality report
- `NanoStats.txt` - Summary statistics
- Various plots (PNG/HTML)

**Example:**
```bash
./01_qc.sh sample.fastq.gz qc_results 16
```

---

### Script 2: Read Alignment

```bash
./02_align.sh <input_fastq> <reference_fasta> <output_bam> <threads>
```

**Parameters:**
- `input_fastq`: FASTQ or FASTQ.GZ file
- `reference_fasta`: Reference genome FASTA
- `output_bam`: Output BAM file path
- `threads`: Number of CPU threads (recommend: 16-32)

**Output:**
- Unsorted BAM file
- Alignment log with statistics

**Example:**
```bash
./02_align.sh sample.fastq.gz reference.fasta aligned.bam 32
```

**Notes:**
- Uses minimap2 with `map-ont` preset for Oxford Nanopore
- Includes MD tags for variant calling
- Uses =/X CIGAR operators for match/mismatch

---

### Script 3: BAM Processing

```bash
./03_process_bam.sh <input_bam> <threads>
```

**Parameters:**
- `input_bam`: Unsorted BAM from alignment
- `threads`: Number of CPU threads

**Output:**
- `<input>.sorted.bam` - Sorted BAM
- `<input>.sorted.bam.bai` - BAM index
- `<input>.sorted.stats.txt` - Alignment statistics

**Example:**
```bash
./03_process_bam.sh aligned.bam 16
```

---

### Script 4a: Variant Calling (Diploid)

```bash
./04a_call_variants_clair3_updated.sh <bam> <reference> <output_dir> <model> <threads>
```

**Parameters:**
- `bam`: Sorted and indexed BAM file
- `reference`: Reference genome FASTA
- `output_dir`: Output directory
- `model`: Clair3 model name (see below)
- `threads`: Number of CPU threads (recommend: 16-32)

**Clair3 Models:**
- `r1041_e82_400bps_sup_v500` - R10.4.1 chemistry, super accuracy (recommended)
- `r1041_e82_400bps_sup_v410` - R10.4.1 chemistry, older
- `r941_prom_sup_g5014` - R9.4.1 chemistry, super accuracy
- `r941_prom_hac_g360+g422` - R9.4.1 chemistry, high accuracy

**Output:**
- `merge_output.vcf.gz` - Final variant calls
- `pileup.vcf.gz` - Pileup-based calls (QC only)
- `full_alignment.vcf.gz` - Full alignment calls (QC only)

**Example:**
```bash
./04a_call_variants_clair3_updated.sh \
    aligned.sorted.bam \
    reference.fasta \
    clair3_output \
    r1041_e82_400bps_sup_v500 \
    32
```

**Use this for:** Standard diploid samples (animals, most plants)

---

### Script 4a: Variant Calling (Polyploid - Integrated)

NOTE: Clair3 reports diploid genotypes (0/1, 1/1, etc) and not polyploid (for example in triploids 0/0/1, 0/1/1, etc).  This can affect the allele frequencies and heterozygosity estimates. The same may be true for Pepper. The script provided here uses FreeBayes after variants are called with Clair3 to re-analyze allele frequencies to estimate triploid genotypes based on read counts. Other tools such as Octopus may be useful for polyploids, but have not been tested. 

```bash
./04a_call_variants_clair3_triploid_v2.sh <bam> <reference> <output_dir> <model> <threads> <ploidy> [min_qual] [min_depth]
```

**Required Parameters:**
- `bam`: Sorted and indexed BAM file
- `reference`: Reference genome FASTA
- `output_dir`: Output directory
- `model`: Clair3 model name
- `threads`: Number of CPU threads
- `ploidy`: Sample ploidy (2, 3, 4, 6, etc.)

**Optional Parameters:**
- `min_qual`: Minimum QUAL score (default: 20)
- `min_depth`: Minimum read depth (default: 5)

**Output Structure:**
```
output_dir/
├── 01_clair3_discovery/
│   └── merge_output.vcf.gz          # Raw Clair3 calls
├── 02_filtered_variants/
│   └── clair3_filtered.vcf.gz       # Filtered by QUAL/DP
└── 03_final_genotypes/
    └── final_ploidy3.vcf.gz         # True polyploid genotypes
```

**Examples:**

Triploid with default filters:
```bash
./04a_call_variants_clair3_triploid_v2.sh \
    banana.sorted.bam \
    banana_ref.fasta \
    banana_variants \
    r1041_e82_400bps_sup_v500 \
    16 \
    3
```

Tetraploid with strict filters:
```bash
./04a_call_variants_clair3_triploid_v2.sh \
    potato.sorted.bam \
    potato_ref.fasta \
    potato_variants \
    r1041_e82_400bps_sup_v500 \
    16 \
    4 \
    30 \
    10
```

Diploid (uses Clair3 only, skips FreeBayes):
```bash
./04a_call_variants_clair3_triploid_v2.sh \
    sample.sorted.bam \
    reference.fasta \
    output \
    r1041_e82_400bps_sup_v500 \
    16 \
    2
```

**Use this for:** 
- Triploid bananas (3n)
- Tetraploid potatoes (4n)
- Hexaploid wheat (6n)
- Any polyploid organism

**Pipeline steps:**
1. **Clair3**: Discovers variant positions (fast, accurate)
2. **Filtering**: Removes low-quality variants
3. **FreeBayes**: Re-genotypes with correct ploidy model

**Genotype interpretations:**
- Diploid (2n): 0/0, 0/1, 1/1
- Triploid (3n): 0/0/0, 0/0/1, 0/1/1, 1/1/1
- Tetraploid (4n): 0/0/0/0, 0/0/0/1, 0/0/1/1, 0/1/1/1, 1/1/1/1
- etc.

---

### Script 4b: Pepper (Not Yet Tested)

### Script 4c: Longshot (Fast SNV Calling)

```bash
./04c_call_variants_longshot_updated.sh <bam> <reference> <output_vcf> <threads>
```

**Parameters:**
- `bam`: Sorted and indexed BAM
- `reference`: Reference genome FASTA
- `output_vcf`: Output VCF path
- `threads`: Number of threads

**Output:**
- VCF with SNVs only (no indels)

**Example:**
```bash
./04c_call_variants_longshot_updated.sh \
    aligned.sorted.bam \
    reference.fasta \
    longshot_snvs.vcf \
    16
```

**Use this for:** 
- Quick SNV calling
- Comparison with Clair3
- When indels are not needed

---

### Script 5: Variant Filtering (NOTE: when using Clair3 for polyploids, the filtering parameters are integrated into the script)

```bash
./05_filter_variants_updated.sh <input_vcf> <output_vcf> <min_qual> <min_depth>
```

**Parameters:**
- `input_vcf`: Input VCF file
- `output_vcf`: Output filtered VCF
- `min_qual`: Minimum QUAL score
- `min_depth`: Minimum read depth

**Recommended Thresholds:**

| Use Case | QUAL | DP | Description |
|----------|------|-----|-------------|
| Exploratory | 15 | 3 | Maximum sensitivity |
| Standard | 20 | 5 | Good balance (default) |
| High confidence | 30 | 10 | Fewer variants, very reliable |
| Publication | 40 | 15 | Strictest filtering |

**Example:**
```bash
# Standard filtering
./05_filter_variants_updated.sh variants.vcf.gz filtered.vcf.gz 20 5

# High confidence
./05_filter_variants_updated.sh variants.vcf.gz high_conf.vcf.gz 30 10
```

**Note:** For polyploid samples, filtering is integrated into the `04a_triploid` script.

---

### Script 6: Variant Statistics

```bash
./06_variant_stats_v2.sh <input_vcf> <reference> <output_prefix>
```

**Parameters:**
- `input_vcf`: VCF file to analyze
- `reference`: Reference genome FASTA
- `output_prefix`: Prefix for output files

**Output Files:**
- `<prefix>.summary.txt` - Overall statistics
- `<prefix>.by_type.txt` - Counts by variant type (SNP, Indel, MNP, Other)
- `<prefix>.by_chromosome.txt` - Distribution by chromosome
- `<prefix>.quality_dist.txt` - QUAL score distribution
- `<prefix>.depth_dist.txt` - Depth distribution

**Variant Types Explained:**
- **SNP**: Single nucleotide change (A→G)
- **Insertion**: Added bases (A→ATT)
- **Deletion**: Removed bases (ATT→A)
- **MNP**: Multiple consecutive changes (ATCG→GGCA)
- **Other**: Complex or symbolic variants

**Example:**
```bash
./06_variant_stats_v2.sh filtered.vcf.gz reference.fasta variant_stats
```

---

## Comparison Tools

NOTE: I have created additional comparison tools that may be useful.  For example, if you want to compare two or more SnpEff annotated VCFs, consider: https://github.com/bjtill/Annotated-VCF-Compare .  Another tool that does not require annotation is found here: https://github.com/bjtill/VCF-Compare-and-Consensus .  The following tool will compare three VCFs and create a venn diagram:  https://github.com/bjtill/ThreeVCF_Analysis_TVA-GUI

### Two-Way VCF Comparison

```bash
./compare_vcfs.sh <vcf1> <vcf2> <output_prefix>
```

**Example:**
```bash
./compare_vcfs.sh clair3.vcf.gz longshot.vcf.gz comparison
```

**Outputs:**
- Overlap statistics
- Unique variants per caller
- Concordance metrics (Jaccard index, F1 score)

### Three-Way VCF Comparison

```bash
./compare_three_vcfs.sh <vcf1> <vcf2> <vcf3> <output_prefix>
```

**Example:**
```bash
./compare_three_vcfs.sh \
    clair3.vcf.gz \
    longshot.vcf.gz \
    pepper.vcf.gz \
    three_way_comparison
```

**Outputs:**
- Venn diagram (text)
- All pairwise and three-way overlaps
- 7 separate VCF files (all combinations)

---

## Polyploid Genotyping Guide

### Supported Ploidy Levels

The pipeline supports **any ploidy level ≥ 2**:

| Ploidy | Organism Examples | Command |
|--------|------------------|---------|
| 2n | Most animals, many plants | `... 2` |
| 3n | Banana, seedless watermelon | `... 3` |
| 4n | Potato, cotton, peanut | `... 4` |
| 6n | Wheat, bread wheat | `... 6` |
| 8n | Strawberry (cultivated) | `... 8` |

### How It Works

1. **Clair3** discovers variants (assumes diploid, but positions are accurate)
2. **Filtering** removes low-quality calls
3. **FreeBayes** re-genotypes at known sites with correct ploidy:
   - Adjusts allele frequency expectations
   - Calculates proper genotype likelihoods
   - Reports correct multi-allelic genotypes

### Understanding Polyploid Genotypes

**Triploid (AAB example):**
- 0/0/0 = AAA (homozygous reference)
- 0/0/1 = AAB (one alternate allele)
- 0/1/1 = ABB (two alternate alleles)  
- 1/1/1 = BBB (homozygous alternate)

**Tetraploid (AABB example):**
- 0/0/0/0 = AAAA
- 0/0/0/1 = AAAB
- 0/0/1/1 = AABB
- 0/1/1/1 = ABBB
- 1/1/1/1 = BBBB

### Expected Variant Distributions in Polyploids

Polyploid samples often show higher MNP counts due to complex haplotype patterns:

| Variant Type | Diploid | Triploid | Tetraploid |
|--------------|---------|----------|------------|
| SNPs | 80-90% | 40-60% | 40-50% |
| Indels | 5-10% | 3-8% | 3-8% |
| MNPs | 3-8% | 30-40% | 35-45% |
| Other | 2-5% | 5-15% | 5-15% |

**This is normal!** Polyploid genomes have complex inheritance patterns that create more MNPs.

---

## Troubleshooting

### Clair3 Hangs or Fails

**Problem:** Conda environment resolution takes forever or Clair3 not found

**Solution:** Use Mamba instead of conda:
```bash
conda install -n base -c conda-forge mamba -y
mamba create -n clair3 -c bioconda clair3 -y
```

### NanoPlot Hangs

**Problem:** NanoPlot hangs when checking version

**Solution:** Reinstall with older version:
```bash
mamba remove nanoplot -y
mamba install -c bioconda "nanoplot<1.42" -y
```

### FreeBayes Takes Too Long

**Problem:** FreeBayes is very slow for large genomes

**Solutions:**
1. **Increase coverage threshold**: Use `--min-coverage 15` instead of 10
2. **Parallelize by chromosome**: Split VCF by chromosome, run in parallel
3. **Use variant sites only**: The integrated script already does this
4. **Skip FreeBayes for exploration**: Use the triploid estimation script instead

### Tab Characters Not Showing

**Problem:** Output files show `\t` instead of tabs

**Solution:** Use the v2 scripts which have this fixed. View files with:
```bash
column -t -s $'\t' file.txt  # Pretty print
cat -A file.txt              # Show all characters
```

### "No such INFO field: DP" Error

**Problem:** Filter script can't find DP field

**Solution:** Use updated filter script (`05_filter_variants_updated.sh`) which auto-detects INFO/DP vs FORMAT/DP

### GQ Type Warning from FreeBayes

**Warning:** `[W::bcf_hdr_check_sanity] GQ should be declared as Type=Integer`

**Answer:** This is harmless. FreeBayes declares GQ as Float instead of Integer. All tools handle this correctly. This warning is suppressed in the v2 scripts.

---

## Best Practices

### Sample Organization

```
project/
├── raw_data/
│   ├── sample1.fastq.gz
│   └── sample2.fastq.gz
├── reference/
│   └── genome.fasta
├── analysis/
│   ├── sample1/
│   │   ├── 01_qc/
│   │   ├── 02_alignment/
│   │   ├── 03_variants/
│   │   └── 04_stats/
│   └── sample2/
│       └── ...
└── scripts/
    └── nanopore_pipeline/
```

### Quality Control Checkpoints

1. **After QC (01)**: Check read quality, length distribution
2. **After alignment (02)**: Verify alignment rate >80%
3. **After filtering (05)**: Check Ti/Tv ratio (expect ~2.0 for most genomes)
4. **After stats (06)**: Review variant type distribution

### Multi-Sample Projects

For projects with many samples:
- Use consistent filtering thresholds
- Document Clair3 model used
- Save all QC metrics for comparison
- Consider joint genotyping for population studies

---

## Performance Expectations

**Typical runtimes for 50x coverage, 3 Gb genome, 32 cores:**

| Step | Time | Peak RAM | Notes |
|------|------|----------|-------|
| QC | 10-30 min | 4 GB | Fast |
| Alignment | 1-3 hours | 16 GB | Fast |
| BAM processing | 30-60 min | 8 GB | Fast |
| Clair3 | 4-12 hours | 32 GB | Slowest step |
| FreeBayes (triploid) | 2-8 hours | 8 GB | Only at known sites |
| Filtering | 5-10 min | 4 GB | Fast |
| Statistics | 10-20 min | 4 GB | Fast |

**Total time (diploid):** ~6-16 hours  
**Total time (triploid):** ~8-22 hours

---

## Citation

If you use this pipeline, please cite:

**Clair3:**
> Zheng et al. (2022). "Symphonizing pileup and full-alignment for deep learning-based long-read variant calling." *Nature Computational Science*, 2, 797-803.

**Minimap2:**
> Li, H. (2018). "Minimap2: pairwise alignment for nucleotide sequences." *Bioinformatics*, 34(18), 3094-3100.

**FreeBayes:**
> Garrison, E. & Marth, G. (2012). "Haplotype-based variant detection from short-read sequencing." *arXiv preprint* arXiv:1207.3907.

**NanoPlot:**
> De Coster, W. et al. (2018). "NanoPack: visualizing and processing long-read sequencing data." *Bioinformatics*, 34(15), 2666-2669.

---

## Support and Issues

For questions or issues:
1. Check the Troubleshooting section above
2. Verify all conda environments are correctly installed
3. Check log files in output directories
4. Ensure adequate disk space and RAM

---

## License

This pipeline is provided as-is for research use. Individual tools have their own licenses:
- Clair3: BSD-3-Clause
- minimap2: MIT
- FreeBayes: MIT
- samtools/bcftools: MIT/BSD
- NanoPlot: MIT

---

## Version History

- **v2.0** (2025-01): Added polyploid support, improved stats, fixed tab formatting


---

**Last Updated:** January 2, 2025
