# Automata-Based Probabilistic DNA Motif Detection
**Compiler Design & Theory of Computation Case Study Project**

An end-to-end bioinformatics pattern detection system implementing formal compiler design principles: **Alphabet ($\Sigma$)**, **Lexical DFA**, **Weighted Finite Automata (WFA)**, **Mismatch-Aware 2D Automata ($q(i, k)$)**, and **Lookahead Branch-and-Bound Pruning**.

---

## 👥 Project Architecture & 4-Member Work Allocation (>500 Lines Per Head)

The codebase is partitioned into **4 modular components** (each with at least 500 lines of verified, production-grade Python code) plus a modern **interactive Graphical User Interface (GUI)**:

| Member | Topic & Module File | Lines of Code | Core Compiler / Biological Responsibilities |
| :--- | :--- | :---: | :--- |
| **Member 1** | **Biological Motif Modeling & 20-Motif Bank**<br>[`module1_motif_models.py`](./module1_motif_models.py) | **867 lines** | • PFM $\to$ PPM $\to$ PSSM Log-Odds transformation.<br>• Shannon Entropy $H(X)$ & Information Content $R_{seq}$ (bits).<br>• Dynamic Programming (Staden Algorithm) for exact $P\text{-value} \leftrightarrow \text{Score}$ mapping.<br>• **20 Verified Transcription Factor Profiles** (GATA1, CTCF, SP1, TATA, NF-$\kappa$B, AP-1, p53, MYC, STAT3, FOXA1, ESR1, SOX2, OCT4, HNF4A, CREB1, E2F1, CEBPA, RUNX1, PAX6, MEF2A). |
| **Member 2** | **Formal Automata & Compiler Engines**<br>[`module2_automata_engine.py`](./module2_automata_engine.py) | **789 lines** | • Classical Exact Lexical DFA Compiler (Regular Expression / IUPAC baseline).<br>• Linear Weighted Finite Automaton (1D-WFA) with diagnostic step tracer.<br>• **Mismatch-Aware 2D-WFA**: $q(i, k)$ state space tracking mismatch budget $K$.<br>• Multi-Motif Trie Automaton & Table-driven Bytecode VM Compiler.<br>• `AutomatonVisualLayout` computing 2D geometric graph coordinates for GUI canvas. |
| **Member 3** | **Genomic DNA Scanner & 20-Genome Registry**<br>[`module3_scanner_pipeline.py`](./module3_scanner_pipeline.py) | **805 lines** | • **20 Preloaded Reference Genomes** (Human Chr11, Chr17, Chr8, Chr3, Chr19, Chr7, Chr6, Chr22, Chr1, Chr12; Mouse Hbb, Oct4; Zebrafish gata1a; Drosophila ftz; Yeast GAL; SARS-CoV-2; Lambda phage; E. coli lac operon; Adenovirus E1A; SV40).<br>• **Add / Append DNA Sequences**: User can add custom genomes or append new DNA fragments to any selected genome.<br>• Dual-Strand (+ and -) scanning & Reverse Chromosome Coordinate Mapping.<br>• Direct PSSM vs Lookahead Pruned WFA (100% equivalence proof).<br>• Standard Exporters: BED6, GFF3, CSV, and JSON. |
| **Member 4** | **Ground-Truth Evaluation & Scalability Profiler**<br>[`module4_eval_benchmark.py`](./module4_eval_benchmark.py) | **612 lines** | • Ground-Truth Synthetic Dataset Generator with biological mutation noise.<br>• Spatial overlap matching & Confusion Matrix (TP, FP, FN, TN).<br>• Precision, Recall, Specificity, F1-Score, MCC, Youden's $J$.<br>• Systematic Threshold Sweeps & ROC / PR Curve integration (AUC-ROC).<br>• Multi-Motif Comparative Benchmark Suite across 20 motifs.<br>• Scalability Profiler (Wall/CPU time, Peak Memory via `tracemalloc`, Throughput in Mbp/s). |
| **Desktop GUI** | **Interactive Modern Tkinter GUI**<br>[`gui_app.py`](./gui_app.py) | **756 lines** | • **Tab 1: Motif & Genome Explorer** (20 Motifs, 20 Genomes, Live Sequence Viewer, DNA Appender).<br>• **Tab 2: Automata Engine Visualizer** (Canvas graph drawing, Step-by-Step Tracer, Bytecode VM).<br>• **Tab 3: Genomic Scanner & Dual-Strand Results** (Interactive Treeview, Context Preview, Exporters).<br>• **Tab 4: Evaluation & Benchmarks** (Ground-truth tests, Confusion Matrix, Multi-Engine table). |
| **Unified** | **Main Pipeline Driver**<br>[`main_pipeline.py`](./main_pipeline.py) | **135 lines** | Orchestrates end-to-end multi-stage pipeline and provides CLI / GUI launcher. |
| **Total** | **4 Comprehensive Modules + GUI + Pipeline** | **3,964 lines** | **Full Production-Grade Case Study** |

---

## 🖥️ Launching the Interactive GUI Application

To launch the graphical desktop application:

```bash
python3 gui_app.py
# or
python3 main_pipeline.py --gui
```

### ✨ Features of the GUI:
1. **Motif Explorer**: Select from 20 transcription factors. View consensus sequence, PSSM matrix, and information content.
2. **Genome Registry**: Select from 20 preloaded reference genomes (human chromosomes, model organisms, viruses, bacteria).
3. **Add or Append DNA Sequence**:
   - **Append Mode**: Concatenate any extra DNA sequence (e.g. `AGATAAGGGGCGGGGGCACGGCACGGCACGG`) directly to the currently selected genome to re-scan it immediately!
   - **Custom Genome Mode**: Create and name your own custom genomic sequence.
4. **Automaton Graph Visualizer**: Interactive canvas rendering states ($q_0 \dots q_L$), weighted transitions, and diagnostic step traces for any custom test window.
5. **High-Throughput Dual-Strand Scanner**: Runs WFA traversal on $(+)$ and $(-)$ strands, displays hits with chromosome coordinates, scores, p-values, mismatches, and surrounding 20 bp nucleotide context.
6. **Bioinformatics Exporters**: Export hits directly into **BED6**, **GFF3**, **CSV**, or **JSON** files.

---

## 🚀 Command Line Execution

### 1. Run the Entire Integrated Pipeline
```bash
python3 main_pipeline.py --demo all
```

### 2. Run Individual Member Demos
```bash
# Member 1 Demo: 20 Motifs, Information Content & Staden DP P-values
python3 module1_motif_models.py

# Member 2 Demo: Formal DFAs, 1D/2D WFAs, Trie Automaton & Bytecode VM
python3 module2_automata_engine.py

# Member 3 Demo: 20 Genomes, Sequence Appending & Dual-Strand Scanner
python3 module3_scanner_pipeline.py

# Member 4 Demo: Ground-Truth benchmarks, Multi-Motif Suite & Profiler
python3 module4_eval_benchmark.py
```

---

## 🧬 Biological Transcription Factors in Motif Bank (20 Profiles)
1. **GATA1** (`MA0035.1`, 6 bp) — Erythroid & megakaryocyte differentiation
2. **CTCF** (`MA0139.1`, 19 bp) — Master 3D chromatin boundary & insulator
3. **SP1** (`MA0079.1`, 10 bp) — Housekeeping GC-box promoter activation
4. **TATA-Box** (`MA0108.1`, 15 bp) — RNA Pol II core promoter initiation
5. **NF-kB** (`MA0107.1`, 10 bp) — Master immune and inflammation regulator
6. **AP-1** (`MA0099.1`, 7 bp) — Stress response and cell proliferation
7. **p53** (`MA0106.1`, 20 bp) — Guardian of the genome tumor suppressor
8. **MYC** (`MA0147.1`, 8 bp) — E-box proto-oncogene transcription factor
9. **STAT3** (`MA0144.1`, 11 bp) — JAK/STAT cytokine signal transducer
10. **FOXA1** (`MA0148.1`, 12 bp) — Pioneer factor for hormone receptor recruitment
11. **ESR1** (`MA0112.1`, 15 bp) — Estrogen nuclear hormone receptor
12. **SOX2** (`MA0143.1`, 10 bp) — Core embryonic stem cell pluripotency
13. **OCT4** (`MA0785.1`, 12 bp) — Master ESC naive pluripotency factor
14. **HNF4A** (`MA0114.1`, 13 bp) — Hepatic metabolism and MODY1 etiology
15. **CREB1** (`MA0018.1`, 8 bp) — cAMP response element binding protein
16. **E2F1** (`MA0024.1`, 11 bp) — Cell cycle G1/S transition checkpoint
17. **CEBPA** (`MA0102.1`, 10 bp) — Granulocytic myeloid differentiation
18. **RUNX1** (`MA0002.1`, 11 bp) — Definitive hematopoiesis & AML factor
19. **PAX6** (`MA0069.1`, 14 bp) — Eye and brain neurodevelopmental regulator
20. **MEF2A** (`MA0052.1`, 10 bp) — Muscle cell differentiation and cardiac control

---

## 🌍 Reference Genomes & Genomic Loci in Registry (20 Loci)
1. **Human Chr11** — Beta-Globin Locus Control Region (HBB, 10,000 bp)
2. **Human Chr17** — TP53 Tumor Suppressor promoter and intron-1 locus (12,000 bp)
3. **Human Chr8** — MYC Oncogene promoter & super-enhancer domain (8,000 bp)
4. **Human Chr3** — SOX2 Core pluripotency regulatory region (9,500 bp)
5. **Human Chr19** — APOE Lipid metabolism and Alzheimer's risk locus (7,500 bp)
6. **Human Chr7** — EGFR Tyrosine kinase promoter region (11,000 bp)
7. **Human Chr6** — HLA-DRA MHC Class II immune regulatory locus (8,500 bp)
8. **Human Chr22** — BCR Chromatin insulator & boundary domain (15,000 bp)
9. **Human Chr1** — RUNX1 Hematopoietic regulatory region (9,000 bp)
10. **Human Chr12** — PAX6 Neuro-ocular master control locus (10,500 bp)
11. **Mouse Chr7** — *Mus musculus* Hbb-b1 adult hemoglobin locus (8,000 bp)
12. **Mouse Chr17** — *Mus musculus* Oct4 (Pou5f1) distal stem-cell enhancer (7,000 bp)
13. **Zebrafish Chr13** — *Danio rerio* gata1a erythroid promoter (6,500 bp)
14. **Fruit Fly Chr3R** — *Drosophila melanogaster* fushi tarazu (ftz) enhancer (5,000 bp)
15. **Baker's Yeast ChrII** — *S. cerevisiae* GAL1-GAL10 upstream activating sequence (4,500 bp)
16. **SARS-CoV-2** — Spike (S) Glycoprotein regulatory sequence (3,822 bp)
17. **Bacteriophage Lambda** — cI / Cro operator regulatory switch (3,200 bp)
18. **E. coli K-12** — lac operon promoter / operator regulatory region (2,800 bp)
19. **Adenovirus Type 5** — E1A early immediate oncoprotein transcription promoter (2,500 bp)
20. **Simian Virus 40** — SV40 Early bidirectional promoter with tandem 21-bp GC repeats (2,200 bp)
