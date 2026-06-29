---
title: "Zebrafish"
type: "species"
---

## Overview

Zebrafish (*Danio rerio*) is a premier model organism for developmental biology, genetics, and regenerative medicine research. This portal provides comprehensive genomic resources including bulk RNA-seq, ATAC-seq, ChIP-seq, single-cell RNA-seq and single-cell ATAC-seq datasets spanning embryonic development, organogenesis, and regeneration processes.

## Research Areas Covered

- **Embryonic Development**: Cell fate determination, pluripotency regulation, maternal-to-zygotic transition
- **Organogenesis**: Heart, brain, pancreas, retina, gonad, and other organ development
- **Hematopoiesis**: Stem cell expansion, lineage specification, blood cell development
- **Regeneration**: Heart, fin, brain, retinal, and muscle regeneration mechanisms
- **Disease Models**: Cancer, cardiovascular disease, toxicology, environmental stress responses
- **Epigenetics**: DNA methylation, histone modifications, chromatin accessibility, 3D genome organization
- **Single-cell Analysis**: Cell-type identification, developmental trajectories, lineage tracing, cell state transitions

## Data Sources

##### Single cell data

| StudyID | Year | Journal | Title | PMID | Sample | GeneID handling |
| ------- | ---- | ------- | ----- | ---- | ------ | --------------- |
| 2025_Elife_Chung | 2025 | Elife | Single-cell transcriptomes of zebrafish germline reveal progenitor types and feminization by Foxl2l | [40497446](https://pubmed.ncbi.nlm.nih.gov/40497446/) | Germ cells | Rebuilt |
| 2025_DevCell_Xiong | 2025 | Dev Cell | Cloche/Npas4l is a pro-regenerative platelet factor during zebrafish heart regeneration | [40602409](https://pubmed.ncbi.nlm.nih.gov/40602409/) | Heart | Converted |
| 2024_NatCellBiol_Lan | 2024 | Nat Cell Biol | Mapping the chromatin accessibility landscape of zebrafish embryogenesis at single-cell resolution by SPATAC-seq | [38977847](https://pubmed.ncbi.nlm.nih.gov/38977847/) | Embryo, scATAC | Converted |
| 2024_Development_Gagnon | 2024 | Development | Germ cell progression through zebrafish spermatogenesis declines with age | [39470160](https://pubmed.ncbi.nlm.nih.gov/39470160/) | Germ cells | Rebuilt |
| 2024_Cell_Royer | 2024 | Cell | A multimodal zebrafish developmental atlas reveals the state-transition dynamics of late-vertebrate pluripotent axial progenitors | [39454574](https://pubmed.ncbi.nlm.nih.gov/39454574/) | Embryo | Converted |
| 2024_BiolOpen_Gagnon | 2024 | Biol Open | Distinct features of the regenerating heart uncovered through comparative single-cell profiling | [38526188](https://pubmed.ncbi.nlm.nih.gov/38526188/) | Heart | Converted |
| 2023_Nature_Trapnell | 2023 | Nature | Embryo-scale reverse genetics at single-cell resolution | [37968389](https://pubmed.ncbi.nlm.nih.gov/37968389/) | Embryo | Converted |
| 2023_DevCell_Farrell | 2023 | Dev Cell | Single-cell analysis of shared signatures and transcriptional diversity during zebrafish development | [37995681](https://pubmed.ncbi.nlm.nih.gov/37995681/) | Embryo | Converted |
| 2022_NatCommun_Crump | 2022 | Nat Commun | Lifelong single-cell profiling of cranial neural crest diversification in zebrafish | [35013168](https://pubmed.ncbi.nlm.nih.gov/35013168/) | Cranial neural crest-derived cells (CNCCs), scRNA+scATAC | Converted |
| 2022_FrontiersinGenetics_Xie | 2022 | Front Genet | Single cell transcriptome sequencing of zebrafish testis revealed novel spermatogenesis marker genes and stronger leydig-germ cell paracrine interactions | [35360857](https://pubmed.ncbi.nlm.nih.gov/35360857/) | Testis | Rebuilt |
| 2022_Elife_Draper | 2022 | Elife | Single-cell transcriptome reveals insights into the development and function of the zebrafish ovary | [35588359](https://pubmed.ncbi.nlm.nih.gov/35588359/) | Germ cells | Rebuilt |
| 2020_Science_Blackshaw | 2020 | Science | Gene regulatory networks controlling vertebrate retinal regeneration | [33004674](https://pubmed.ncbi.nlm.nih.gov/33004674/) | Retina | Converted |
| 2020_EnvironSciTechnol_Wu | 2020 | Environ Sci Technol | Single-cell sequencing reveals heterogeneity effects of bisphenol A on zebrafish embryonic development | [32644799](https://pubmed.ncbi.nlm.nih.gov/32644799/) | Embryo | Converted |
| 2020_EMBORep_Sumeet | 2020 | EMBO Rep | Single-cell transcriptome analysis reveals thyrocyte diversity in the zebrafish thyroid gland | [33140917](https://pubmed.ncbi.nlm.nih.gov/33140917/) | Thyroid | Converted |
| 2020_DevCell_Riley | 2020 | Dev Cell | Functional heterogeneity within the developing zebrafish epicardium | [32084358](https://pubmed.ncbi.nlm.nih.gov/32084358/) | Epicardium | Converted |
| 2020_DevBiol_Miller | 2020 | Dev Biol | A single-cell transcriptome atlas for zebrafish development | [31782996](https://pubmed.ncbi.nlm.nih.gov/31782996/) | Embryo | Converted |
| 2019_CellRep_Liu | 2019 | Cell Rep | A 3D atlas of hematopoietic stem and progenitor cell expansion by multi-dimensional RNA-seq analysis | [31042481](https://pubmed.ncbi.nlm.nih.gov/31042481/) | Caudal hematopoietic tissue (CHT) | Converted |
| 2018_Science_Schier | 2018 | Science | Single-cell reconstruction of developmental trajectories during zebrafish embryogenesis | [29700225](https://pubmed.ncbi.nlm.nih.gov/29700225/) | Embryo | Converted |
| 2018_Science_Klein | 2018 | Science | Single-cell mapping of gene expression landscapes and lineage in the zebrafish embryo | [29700229](https://pubmed.ncbi.nlm.nih.gov/29700229/) | Embryo | Converted |
| 2018_NatBiotechnol_Schier | 2018 | Nat Biotechnol | Simultaneous single-cell profiling of lineages and cell types in the vertebrate brain | [29608178](https://pubmed.ncbi.nlm.nih.gov/29608178/) | Brain | Converted |

`GeneID handling` describes how single-cell gene identifiers are connected to the current portal database:

- `Original`: keeps the author’s original identifier system.
- `Converted`: keeps the author matrix but converts gene identifiers to the current database when a stable mapping is available. Some features may be omitted from public display when no confident mapping exists, and renamed or merged annotations can introduce small differences relative to the author’s original feature set.
- `Rebuilt`: regenerates the matrix from raw reads against the current reference. This keeps identifiers consistent with the current database, but expression values, detected features, and cell annotations can differ from the author’s original workflow.

<div style="height: 0.75rem;"></div>

##### Bulk data

| StudyID | Year | Journal | Title | PMID | Sample |
| ------- | ---- | ------- | ----- | ---- | ------ |
| 2024_DevCell_Simoes | 2024 | Dev Cell | Distinct epicardial gene regulatory programs drive development and regeneration of the zebrafish heart | [38237592](https://pubmed.ncbi.nlm.nih.gov/38237592/) | Epicardium control, larval, sham and cryoinjury states in tbx18 and tcf21 lineages |
| 2022_NatCommun_Onichtchouk | 2022 | Nat Commun | Pluripotency factors determine gene expression repertoire at zygotic genome activation | [35145080](https://pubmed.ncbi.nlm.nih.gov/35145080/) | Embryo 2.5 to 6 hpf, dome and oblong |
| 2022_NatCommun_Bessa | 2022 | Nat Commun | Multidimensional chromatin profiling of zebrafish pancreas to uncover and investigate disease-relevant enhancers | [35410466](https://pubmed.ncbi.nlm.nih.gov/35410466/) | Endocrine, exocrine, muscle and pancreas |
| 2022_Development | 2022 | Development | Wt1 transcription factor impairs cardiomyocyte specification and drives a phenotypic switch from myocardium to epicardium | [35312773](https://pubmed.ncbi.nlm.nih.gov/35312773/) | Heart tube, proepicardium, pericardium and wt1b overexpression |
| 2022_BMCGenomics | 2022 | BMC Genomics | High-throughput transcriptome sequencing reveals the key stages of cardiovascular development in zebrafish embryos | [35964013](https://pubmed.ncbi.nlm.nih.gov/35964013/) | Embryo 18, 24 and 42 hpf |
| 2021_Development | 2021 | Development | Extensive nuclear gyration and pervasive non-genic transcription during primordial germ cell development in zebrafish | [33298460](https://pubmed.ncbi.nlm.nih.gov/33298460/) | Embryo 0 hpf to 10 dpf and PGC RNA fractions |
| 2021_DevCell | 2021 | Dev Cell | Germ cell differentiation requires Tdrd7-dependent chromatin and transcriptome reprogramming marked by germ plasm relocalization | [33651978](https://pubmed.ncbi.nlm.nih.gov/33651978/) | PGC, soma and brain from 256-cell to prim5 and 10 somites |
| 2020_iScience | 2020 | iScience | Demarcation of topologically associating domains is uncoupled from enriched CTCF binding in developing zebrafish | [32334414](https://pubmed.ncbi.nlm.nih.gov/32334414/) | Embryo 24 hpf |
| 2020_eLife_Zhu | 2020 | eLife | An improved zebrafish transcriptome annotation for sensitive and comprehensive detection of cell type-specific genes | [32831172](https://pubmed.ncbi.nlm.nih.gov/32831172/) | Nr2f2/Y1 and pdgfrb positive and negative cell populations |
| 2020_eLife_Cornell | 2020 | eLife | Analysis of zebrafish periderm enhancers facilitates identification of a regulatory variant near human KRT8/18 | [32031521](https://pubmed.ncbi.nlm.nih.gov/32031521/) | krt4-GFP positive and negative embryo cells at 11 hpf |
| 2020_Science_SanchezAlvarado | 2020 | Science | Changes in regeneration-responsive enhancers shape regenerative capacities in vertebrates | [32883834](https://pubmed.ncbi.nlm.nih.gov/32883834/) | Fin regeneration 0 and 1 dpa, RNA and chromatin marks |
| 2020_Nature | 2020 | Nature | A map of cis-regulatory elements and 3D genome structures in zebrafish | [33239788](https://pubmed.ncbi.nlm.nih.gov/33239788/) | Adult tissues plus embryonic brain and trunk |
| 2020_GenomeBiol | 2020 | Genome Biol | Regenerating zebrafish fin epigenome is characterized by stable lineage-specific DNA methylation and dynamic chromatin accessibility | [32106888](https://pubmed.ncbi.nlm.nih.gov/32106888/) | Fin 0, 1, 2 and 4 dpa, sp7-positive and sp7-negative subsets |
| 2020_CirculationRes | 2020 | Circ Res | AP-1 contributes to chromatin accessibility to promote sarcomere disassembly and cardiomyocyte protrusion during zebrafish heart regeneration | [32312172](https://pubmed.ncbi.nlm.nih.gov/32312172/) | Heart 4 dpci A-Fos and gata4 states, uninjured overexpression controls |
| 2019_PLoSOne | 2019 | PLoS One | Transcriptomic profile of early zebrafish PGCs by single cell sequencing | [31412047](https://pubmed.ncbi.nlm.nih.gov/31412047/) | PGC 6, 11 and 24 hpf |
| 2019_NatCommun_Tena | 2019 | Nat Commun | Pioneer and repressive functions of p63 during zebrafish embryonic ectoderm specification | [31296872](https://pubmed.ncbi.nlm.nih.gov/31296872/) | Embryo 24 and 36 hpf, 8 somites, 80% epiboly |
| 2019_NatCommun_Bogdanovic | 2019 | Nat Commun | Retention of paternal DNA methylome in the developing zebrafish germline | [31296860](https://pubmed.ncbi.nlm.nih.gov/31296860/) | PGC and soma at 4 to 36 hpf, adult liver |
| 2019_InvestOphthalmolVisSci | 2019 | Invest Ophthalmol Vis Sci | Paradoxical changes underscore epigenetic reprogramming during adult zebrafish extraocular muscle regeneration | [31794598](https://pubmed.ncbi.nlm.nih.gov/31794598/) | Muscle control, 18 hpi, 9 dpi |
| 2019_GenomeRes | 2019 | Genome Res | Dynamics of cardiomyocyte transcriptome and chromatin landscape demarcates key events of heart development | [30760547](https://pubmed.ncbi.nlm.nih.gov/30760547/) | Heart GFP positive and negative cells, 24 to 72 hpf, mutants |
| 2019_Development | 2019 | Development | H3K27me3-mediated silencing of structural genes is required for zebrafish heart regeneration | [31427288](https://pubmed.ncbi.nlm.nih.gov/31427288/) | Heart and cardiomyocytes, injured and uninjured |
| 2019_DevBiol | 2019 | Dev Biol | foxc1 is required for embryonic head vascular smooth muscle differentiation in zebrafish | [31199900](https://pubmed.ncbi.nlm.nih.gov/31199900/) | acta2 and kdrl positive and negative vascular cells |
| 2019_BMCBiol | 2019 | BMC Biol | FoxH1 represses miR-430 during early embryonic development of zebrafish via non-canonical regulation | [31362746](https://pubmed.ncbi.nlm.nih.gov/31362746/) | Embryo 6 hpf |
| 2018_PLoSGenet | 2018 | PLoS Genet | Cooperation, cis-interactions, versatility and evolutionary plasticity of multiple cis-acting elements underlie krox20 hindbrain regulation | [30080860](https://pubmed.ncbi.nlm.nih.gov/30080860/) | Hindbrain and posterior tissues at 5 and 15 ss, 95% epiboly |
| 2018_Nature | 2018 | Nature | Amphioxus functional genomics and the origins of vertebrate gene regulation | [30464347](https://pubmed.ncbi.nlm.nih.gov/30464347/) | Embryo 2 to 72 hpf |
| 2018_GenomeRes | 2018 | Genome Res | Inherited DNA methylation primes the establishment of accessible chromatin during genome activation | [29844026](https://pubmed.ncbi.nlm.nih.gov/29844026/) | Embryo 64-cell to oblong, dome, 1k-cell |
| 2018_FrontEndocrinol_Lausanne | 2018 | Front Endocrinol (Lausanne) | Transcriptomic analysis for differentially expressed genes in ovarian follicle activation in the zebrafish | [30364302](https://pubmed.ncbi.nlm.nih.gov/30364302/) | Ovary PG and PV follicles |
| 2018_Development | 2018 | Development | Cohesin facilitates zygotic genome activation in zebrafish | [29158440](https://pubmed.ncbi.nlm.nih.gov/29158440/) | Embryo 2.5 to 10 hpf |
| 2018_Cell | 2018 | Cell | Placeholder nucleosomes underlie germline-to-embryo DNA methylation reprogramming | [29456083](https://pubmed.ncbi.nlm.nih.gov/29456083/) | Pre-ZGA and post-ZGA embryo, sperm |
| 2017_eLife_Stainier | 2017 | eLife | Reciprocal analyses in zebrafish and medaka reveal that harnessing the immune response promotes cardiac regeneration | [28632131](https://pubmed.ncbi.nlm.nih.gov/28632131/) | Heart regeneration 0 h to 5 dpci/dps, 6 hpci/hps |
| 2017_eLife_Giraldez | 2017 | eLife | A high-resolution mRNA expression time course of embryonic development in zebrafish | [29144233](https://pubmed.ncbi.nlm.nih.gov/29144233/) | Embryo 0 hpf to 5 dpf |
| 2017_Development | 2017 | Development | Heart morphogenesis gene regulatory networks revealed by temporal expression analysis | [28807900](https://pubmed.ncbi.nlm.nih.gov/28807900/) | Heart 30 to 72 hpf |
| 2017_DevCell | 2017 | Dev Cell | The vertebrate protein dead end maintains primordial germ cell fate by inhibiting somatic differentiation | [29257950](https://pubmed.ncbi.nlm.nih.gov/29257950/) | PGC 13 hpf |
| 2017_CellRep | 2017 | Cell Rep | Robust identification of developmentally active endothelial enhancers in zebrafish using FANS-assisted ATAC-seq | [28723572](https://pubmed.ncbi.nlm.nih.gov/28723572/) | fli1a-GFP positive and negative cells at 24 hpf |
| 2017_BMCGenomics | 2017 | BMC Genomics | Histological and transcriptomic effects of 17α-methyltestosterone on zebrafish gonad development | [28738802](https://pubmed.ncbi.nlm.nih.gov/28738802/) | Ovary and testis at 40 and 60 dpf |
| 2016_NatCommun | 2016 | Nat Commun | Repulsive cues combined with physical barriers and cell-cell adhesion determine progenitor cell positioning during organogenesis | [27088892](https://pubmed.ncbi.nlm.nih.gov/27088892/) | PGC and soma at 7 and 36 hpf |
| 2015_PNAS | 2015 | Proc Natl Acad Sci U S A | Deep conservation of wrist and digit enhancers in fish | [25535365](https://pubmed.ncbi.nlm.nih.gov/25535365/) | Embryo 24 hpf, fin 60 hpf |
| 2015_GenomicsData | 2015 | Genom Data | Genome-wide epigenetic cross-talk between DNA methylation and H3K27me3 in zebrafish embryos | [26697317](https://pubmed.ncbi.nlm.nih.gov/26697317/) | Embryo dome, 24 and 48 hpf |
| 2014_BMCBiol | 2014 | BMC Biol | Global identification of Smad2 and Eomesodermin targets in zebrafish identifies a conserved transcriptional network in mesendoderm and a novel role for Eomesodermin in repression of ectodermal gene expression | [25277163](https://pubmed.ncbi.nlm.nih.gov/25277163/) | Embryo blastula |
| 2013_StemCellReports | 2013 | Stem Cell Reports | A Cdx4-Sall4 regulatory module controls the transition from mesoderm formation to embryonic hematopoiesis | [24286030](https://pubmed.ncbi.nlm.nih.gov/24286030/) | Embryo bud stage |
| 2013_Science | 2013 | Science | Pou5f1 transcription factor controls zygotic gene activation in vertebrates | [23950494](https://pubmed.ncbi.nlm.nih.gov/23950494/) | Embryo 2.75 and 5 hpf |
| 2013_PLoSGenet | 2013 | PLoS Genet | Genome wide analysis reveals Zic3 interaction with distal regulatory elements of stage specific developmental genes in zebrafish | [24204288](https://pubmed.ncbi.nlm.nih.gov/24204288/) | Embryo 8 and 24 hpf |
| 2013_Cell | 2013 | Cell | Reprogramming the maternal zebrafish genome after fertilization to match the paternal methylation pattern | [23663776](https://pubmed.ncbi.nlm.nih.gov/23663776/) | Early embryo, egg, sperm, muscle |
| 2013_BMCGenomics | 2013 | BMC Genomics | Differential transcript isoform usage pre- and post-zygotic genome activation in zebrafish | [23676078](https://pubmed.ncbi.nlm.nih.gov/23676078/) | Embryo post-MBT |
| 2012_GenomeRes_Schier | 2012 | Genome Res | Systematic identification of long noncoding RNAs expressed during zebrafish embryogenesis | [22110045](https://pubmed.ncbi.nlm.nih.gov/22110045/) | Embryo 0.75 hpf to 5 dpf |
| 2012_GenomeRes_GomezSkarmeta | 2012 | Genome Res | Dynamics of enhancer chromatin signatures mark the transition from pluripotency to cell specification during embryogenesis | [22593555](https://pubmed.ncbi.nlm.nih.gov/22593555/) | Embryo dome, 80% epiboly, 24 and 48 hpf |
| 2012_GenomeRes_Bartel | 2012 | Genome Res | Extensive alternative polyadenylation during zebrafish development | [22722342](https://pubmed.ncbi.nlm.nih.gov/22722342/) | Embryo 1-cell to 72 hpf, adult brain, ovary, testis |
| 2012_DevCell | 2012 | Dev Cell | Nanog-like regulates endoderm formation through the Mxtx2-Nodal pathway | [22421047](https://pubmed.ncbi.nlm.nih.gov/22421047/) | Embryo 3.5 and 4.5 hpf |
| 2012_DevBiol | 2012 | Dev Biol | Zebrafish globin switching occurs in two developmental stages and is controlled by the LCR | [22537494](https://pubmed.ncbi.nlm.nih.gov/22537494/) | Adult erythroid cells and whole extract |
| 2011_DevBiol | 2011 | Dev Biol | Identification of cis regulatory features in the embryonic zebrafish genome through large-scale profiling of H3K4me1 and H3K4me3 binding sites | [21435340](https://pubmed.ncbi.nlm.nih.gov/21435340/) | Embryo 24 hpf |
| 2011_Cell | 2011 | Cell | Conserved function of lincRNAs in vertebrate embryonic development despite rapid sequence evolution | [22196729](https://pubmed.ncbi.nlm.nih.gov/22196729/) | Embryo 24 and 72 hpf, adult |
