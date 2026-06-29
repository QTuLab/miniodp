---
title: "Turquoise killifish"
type: "species"
---

## Overview

Turquoise killifish (*Nothobranchius furzeri*), also known as the African turquoise killifish, is a short-lived annual vertebrate from ephemeral pools in southeastern Africa. Rapid growth, early maturation, and embryonic diapause are central parts of its life cycle, and the species is widely used in studies of aging, developmental timing, regeneration, and life-history evolution.

The current data release combines baseline bulk transcriptomes with curated single-cell references. Public datasets provide a compact reference set for gene-centered exploration across sex, tissue, age, strain, and cell type in this species.

## Research Areas Covered

- **Aging and longevity**: age-series transcriptomes from multiple tissues, including brain, liver, skin, skeletal muscle, and fin clip.
- **Diapause and developmental timing**: a vertebrate system for studying suspended animation, embryo dormancy, and life-cycle transitions.
- **Tissue homeostasis and regeneration**: age-associated molecular changes in organs with strong regenerative or maintenance phenotypes.
- **Single-cell tissue atlas and sex dimorphism**: cell-type-resolved profiles from blood, kidney, liver, spleen, and brain support analysis of immune, metabolic, neural, and sex-biased programs.
- **Comparative genomics and life-history evolution**: a compact vertebrate model for linking genome evolution, sex chromosome biology, and rapid aging.

## Data Sources

##### Single cell data

| StudyID | Year | Journal | Title | PMID | Accessions | Sample | GeneID handling |
| ------- | ---- | ------- | ----- | ---- | ---------- | ------ | --------------- |
| 2026_SciData_Xu | 2026 | Sci Data | A single-cell RNA sequencing dataset of cardiac aging in African Turquoise Killifish | [42162007](https://pubmed.ncbi.nlm.nih.gov/42162007/) | [PRJNA1405634](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1405634); [GSE316983](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316983) | Heart | Converted (99.8%) |
| 2026_BioRxiv_Benayoun | 2026 | BioRxiv | A multi-omic atlas in the African turquoise killifish reveals increased glucocorticoid signaling as a hallmark of brain aging | [41993562](https://pubmed.ncbi.nlm.nih.gov/41993562/) | [PRJNA1105049](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1105049); [Figshare](https://doi.org/10.6084/m9.figshare.31847155) | Brain | Converted (81%) |
| 2026_KidneyInt_Paulmann | 2026 | Kidney Int | Sodium-glucose co-transporter 2 inhibition improves age-dependent kidney microvascular rarefaction | [41448458](https://pubmed.ncbi.nlm.nih.gov/41448458/) | [PRJNA1265492](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1265492); [GSE297623](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297623) | Kidney | Converted (62%) |
| 2024_AgingCell_Moons | 2024 | Aging Cell | Age-related dysregulation of the retinal transcriptome in African turquoise killifish | [38742929](https://pubmed.ncbi.nlm.nih.gov/38742929/) | [PRJNA1074745](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1074745); [GSE255364](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255364) | Retina | Converted (70%) |
| 2023_CellRep_Benayoun | 2023 | Cell Rep | Widespread sex dimorphism across single-cell transcriptomes of adult African turquoise killifish tissues | [37837621](https://pubmed.ncbi.nlm.nih.gov/37837621/) | [PRJNA952805](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA952805); [SCP2220](https://singlecell.broadinstitute.org/single_cell/study/SCP2220) | Blood, kidney, liver, spleen | Converted (83%) |

`GeneID handling` describes how single-cell gene identifiers are connected to the current portal database:

- `Original`: keeps the author’s original identifier system.
- `Converted`: keeps the author matrix but converts gene identifiers to the current miniodp database when a stable mapping is available. The percentage indicates the fraction of author features retained in the public matrix after stable mapping to the current miniodp database. Features without a confident mapping are excluded from public display, and renamed or merged annotations can introduce small differences relative to the author’s original feature set.
- `Rebuilt`: regenerates the matrix from raw reads against the current reference. This keeps identifiers consistent with the current miniodp database, but expression values, detected features, and cell annotations can differ from the author’s original workflow.

<div style="height: 0.75rem;"></div>

##### Bulk data

| StudyID | Year | Journal | Title | PMID | Accessions | Sample |
| ------- | ---- | ------- | ----- | ---- | ---------- | ------ |
| 2026_NatAging_WyssCoray | 2026 | Nat Aging | Multi-tissue transcriptomic aging atlas reveals predictive aging biomarkers in the killifish | [41776309](https://pubmed.ncbi.nlm.nih.gov/41776309/) | [PRJNA1274512](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1274512); [GSE308970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308970) | Multi-tissue atlas across bone, brain, fat, gut, heart, kidney, liver, muscle, ovary, skin, spinal cord, spleen, testis by sex and age; excludes FACS-sorted head kidney add-on samples |
| 2026_BioRxiv_Benayoun_ATAC | 2026 | BioRxiv | A multi-omic atlas in the African turquoise killifish reveals increased glucocorticoid signaling as a hallmark of brain aging | [41993562](https://pubmed.ncbi.nlm.nih.gov/41993562/) | [PRJNA1100604](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1100604); [SRP501900](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP501900) | ATAC: brain 06wk female, 06wk male, 16wk female, 16wk male |
| 2026_SciRep_Nishina | 2026 | Sci Rep | Identification of gradual aging and late-onset aging markers using male African turquoise killifish | [42020722](https://pubmed.ncbi.nlm.nih.gov/42020722/) | [PRJNA1399382](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1399382); [GSE315825](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315825) | Male liver young, adult, old |
| 2025_NatCommun_Englert | 2025 | Nat Commun | The master male sex determinant Gdf6Y of the turquoise killifish arose through allelic neofunctionalization | [39788971](https://pubmed.ncbi.nlm.nih.gov/39788971/) | [PRJEB5837](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB5837); [GSE263626](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263626) | Embryo 10dpf, trunk 0dph, trunk 3dph, gonad 3mph |
| 2024_Cell_Brunet_ATAC | 2024 | Cell | Evolution of diapause in the African turquoise killifish by remodeling the ancient gene regulatory landscape | [38810644](https://pubmed.ncbi.nlm.nih.gov/38810644/) | [PRJNA770947](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770947); [GSE185816](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185816) | ATAC: embryo pre-diapause, development, diapause day6, diapause 1month |
| 2024_AgingCell_Moons | 2024 | Aging Cell | Age-related dysregulation of the retinal transcriptome in African turquoise killifish | [38742929](https://pubmed.ncbi.nlm.nih.gov/38742929/) | [PRJNA1074746](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1074746); [GSE255363](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255363) | Retina 06wk, 12wk, 18wk |
| 2023_SciData_Benayoun | 2023 | Sci Data | Transcriptomes of aging brain, heart, muscle, and spleen from female and male African turquoise killifish | [37828039](https://pubmed.ncbi.nlm.nih.gov/37828039/) | [PRJNA952180](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA952180) | Brain, heart, muscle, spleen by sex at 06wk and 16wk |
| 2023_GenomeRes_Benayoun | 2023 | Genome Res | Dynamic regulation of gonadal transposon control across the lifespan of the naturally short-lived African turquoise killifish | [36577520](https://pubmed.ncbi.nlm.nih.gov/36577520/) | [PRJNA854614](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA854614) | Ovary and testis at young, middle age, old |
| 2020_Science_Brunet_ChIP | 2020 | Science | Vertebrate diapause preserves organisms long term through Polycomb complex members | [32079766](https://pubmed.ncbi.nlm.nih.gov/32079766/) | [PRJNA503701](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA503701) | ChIP: diapause and non-diapause embryo H3K27me3, H3K4me3, input |
| 2020_Science_Brunet | 2020 | Science | Vertebrate diapause preserves organisms long term through Polycomb complex members | [32079766](https://pubmed.ncbi.nlm.nih.gov/32079766/) | [PRJNA503701](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA503701) | Pre-diapause, diapause day3, diapause day6, diapause 1month, without diapause |
| 2020_Science_SanchezAlvarado_ChIP | 2020 | Science | Changes in regeneration-responsive enhancers shape regenerative capacities in vertebrates | [32883834](https://pubmed.ncbi.nlm.nih.gov/32883834/) | [PRJNA559885](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA559885); [SRP265421](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP265421) | ChIP: fin regeneration 0dpa and 1dpa H3K27ac, H3K4me3, input |
| 2020_Science_SanchezAlvarado | 2020 | Science | Changes in regeneration-responsive enhancers shape regenerative capacities in vertebrates | [32883834](https://pubmed.ncbi.nlm.nih.gov/32883834/) | [PRJNA559885](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA559885); [SRP265421](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP265421) | Fin regeneration 0dpa, 3hpa, 6hpa, 14hpa, 1dpa, 2dpa, 3dpa, 4dpa, 7dpa, 18dpa |
| 2020_MolSystBiol_Ori | 2020 | Mol Syst Biol | Reduced proteasome activity in the aging brain results in ribosome stoichiometry loss and aggregation | [32558274](https://pubmed.ncbi.nlm.nih.gov/32558274/) | [PRJNA631760](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA631760); [GSE150318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150318) | Fin clip 10wk, 20wk |
| 2019_Cells_Gaetano_ChIP | 2019 | Cells | Aging Triggers H3K27 Trimethylation Hoarding in the Chromatin of *Nothobranchius furzeri* Skeletal Muscle | [31569376](https://pubmed.ncbi.nlm.nih.gov/31569376/) | [PRJNA557479](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA557479); [GSE135129](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135129) | ChIP: skeletal muscle young and old H3K27me3, H3K9ac |
| 2019_Cells_Gaetano | 2019 | Cells | Aging Triggers H3K27 Trimethylation Hoarding in the Chromatin of *Nothobranchius furzeri* Skeletal Muscle | [31569376](https://pubmed.ncbi.nlm.nih.gov/31569376/) | [PRJNA557199](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA557199); [GSE135032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135032) | Skeletal muscle 05wk, 20wk, 35wk |
| 2017_JenAge_Brain | 2017 | - | Sequencing of *Nothobranchius furzeri* (strain: GRZ) brain in different age groups | - | [PRJNA400236](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA400236); [GSE103132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103132) | Brain 05wk, 07wk, 10wk, 14wk |
| 2017_JenAge_Skin | 2017 | - | Sequencing of *Nothobranchius furzeri* (strain: GRZ) skin in different age groups | - | [PRJNA400244](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA400244); [GSE103140](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103140) | Skin 05wk, 07wk, 10wk, 14wk |
| 2017_JenAge_Liver | 2017 | - | Sequencing of *Nothobranchius furzeri* (strain: GRZ) liver in different age groups | - | [PRJNA400247](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA400247); [GSE103137](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103137) | Liver 05wk, 07wk, 10wk, 14wk |
| 2015_Cell_Platzer_CrossSection | 2015 | Cell | Insights into Sex Chromosome Evolution and Aging from the Genome of a Short-Lived Fish | [26638077](https://pubmed.ncbi.nlm.nih.gov/26638077/) | [PRJNA284603](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA284603); [GSE69122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69122) | Brain, liver, skin at 12wk |
| 2015_Cell_Platzer_AgeSeries | 2015 | Cell | Insights into Sex Chromosome Evolution and Aging from the Genome of a Short-Lived Fish | [26638077](https://pubmed.ncbi.nlm.nih.gov/26638077/) | [PRJNA277768](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA277768); [GSE66712](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66712) | Liver and skin at 05wk, 12wk, 20wk, 27wk, 39wk |
| 2014_AgingCell_Cellerino | 2014 | Aging Cell | RNA-seq of the aging brain in the short-lived fish *N. furzeri* - conserved pathways and novel genes associated with neurogenesis | [25059688](https://pubmed.ncbi.nlm.nih.gov/25059688/) | [PRJNA229052](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA229052); [GSE52462](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52462) | Brain 05wk, 12wk, 20wk, 27wk, 39wk |
