# Isoform Switch Analysis in Hürthle Cell Carcinoma (HCC)

[![Status](https://img.shields.io/badge/Analysis-In--Progress-orange.svg)]()
> **Current Status:** Functional annotation steps (Pfam, DeepTMHMM) are currently running on the full dataset. Downstream isoform switch metrics are being integrated.

## Background & Motivation

Hürthle cell carcinoma (HCC) is a subtype of thyroid cancer, accounting for 3-5% of all thyroid malignancies. HCC is characterised by an abundance of malfunctioning mitochondria and poor response to radioiodine therapy. While prior studies have documented mitochondrial complex I DNA mutations and metabolomic vulnerabilities in HCC, the transcriptomic landscape 
remains largely unexplored. No studies to date have specifically characterised isoform switching events in HCC.

The analysis is performed on NCBI GEO dataset and explores the expression profiles of different isoforms in HCC tissues. The original study explored metabolomic profiles of HCC and identified that mitochondrial complex I loss along with lipid peroxide stress is a vulnerability in HCC. This analysis performed on a subset of samples explores the isoform profiles in HCC, identifies some major genes undergoing functional isoform switching and highlights alternative splicing mechanisms that may drive the pathogenesis of the disease. 

## Objective

* Quantify the expression of transcripts between normal & HCC tissues.
* Identify transcript isoforms in HCC tissues. 
* Predict functional consequences & alternative splicing event mechanisms in HCC tissues. 

## Dataset

The dataset for this analysis have been obtained from NCBI GEO with accession ID [GSE228870](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228870).
The dataset was subset to include 7 HCC samples & 7 matched normal thyroid samples. 

## Methodology 

1. **Data Import**

- Download paired-end raw reads (fastq files) from NCBI GEO. 

- Tool: fasterq-dump

2. **Initial QC**

- Check the quality of raw reads, including per base sequence quality, adapter sequences, GC content etc. 

- Tool: fastqc

3. **QC**

- All-in-one processing to remove low quality sequences, over represented sequences & adapter trimming. 

- Tool: fastp

4. **Quantification**

- Mapping & quantification of reads against a reference transcriptome (GRCh38). 

- Tool: Salmon

5. **Post-Alignment QC**

- Check the mapping quality & mapping rate of reads against the reference transcriptome. 

- Tool: MultiQC

6. **Isoform Switch Analysis**

- Isoform switches in tumor samples along with their functional consequences including Non-sense mediated decay (NMD) sensitivity, intron retention & coding potential of isoform was analyzed. 

- Tool: IsoformSwitchAnalyzeR

7. **Visualization**

- Statistically significant switches, alternative splicing events & consequence summary for different genes was visualized.

- Tool: IsoformSwitchAnalyzeR

## Analytical Decisions & Rationale

**Why 7 matched pairs instead of the full dataset?**

The full dataset contains samples lacking matched normal controls. For isoform switch analysis, tumour-normal pairing is essential to study transcriptomic variation. Computational constraints also informed this decision. Only samples with confirmed matched pairs were retained, yielding 7 HCC and 7 normal thyroid samples.

**Why Salmon over HISAT2 + featureCounts?**

Isoform-level quantification requires transcript-level resolution. HISAT2 + featureCounts is optimised for gene-level count matrices and would have required additional assembly steps (e.g. StringTie) to recover novel isoforms, which was outside the scope of this analysis. Salmon's quasi-mapping approach quantifies directly against the reference transcriptome at transcript resolution, is computationally efficient, and its output integrates directly with IsoformSwitchAnalyzeR. `--validateMappings` was enabled to improve mapping accuracy by removing invalid multi-mapping reads.

**Why GENCODE v44 as reference?**

GENCODE v44 provides the most comprehensive human transcript annotation, including complete genome assembly (GRCh38.p14) and a matched GTF file required for transcript-level quantification. Using a comprehensive annotation is critical for isoform analysis.

**Why isoform-level analysis over standard DEG?**

Standard DEG analysis (e.g. DESeq2 on gene-level counts) would not capture scenarios where total gene expression is stable but isoform usage shifts — a pattern that can have profound functional consequences. DIXDC1 is a direct example from this dataset: gene-level analysis would have missed it entirely.

## Results

This section will be updated once analysis is completed. 


## References
- [GEO DATASET](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228870)
- [Salmon Documentation](https://salmon.readthedocs.io/en/latest/salmon.html)
- [IsoformSwitchAnalyzeR Tool Paper](https://academic.oup.com/bioinformatics/article/35/21/4469/5466456)
- [Effectors Enabling Adaptation to Mitochondrial Complex I Loss in Hürthle Cell Carcinoma](https://pubmed.ncbi.nlm.nih.gov/37262067/)
- [The Molecular Landscape of Hürthle Cell Thyroid Cancer Is Associated with Altered Mitochondrial Function—A Comprehensive Review](https://www.mdpi.com/2073-4409/9/7/1570)
