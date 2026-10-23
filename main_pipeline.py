"""
Automata-Based Probabilistic DNA Motif Detection: Unified Main Pipeline Runner
==============================================================================
Compiler Design & Theory of Computation Case Study Project

Team Members & Architecture (4 Members - 500+ Lines Each):
- Member 1: Biological Motif Modeling, Information Theory & 20-Motif Bank (module1_motif_models.py)
- Member 2: Formal Automata & Compiler Engines: Linear WFA, 2D-WFA, DFA, Trie VM (module2_automata_engine.py)
- Member 3: Genomic DNA Scanning Engine, 20-Genome Registry & Coordinate Pipeline (module3_scanner_pipeline.py)
- Member 4: Ground-Truth Evaluation Framework & Scalability Profiler (module4_eval_benchmark.py)
- Interactive Desktop GUI Suite: (gui_app.py)
"""

from __future__ import annotations

import argparse
import os
import sys
import time

from module1_motif_models import MotifBank, BackgroundModel, BackgroundModelType
from module2_automata_engine import LinearWFA, MismatchAwareWFA, ClassicalExactDFA
from module3_scanner_pipeline import GenomeRegistry, DualStrandScannerPipeline
from module4_eval_benchmark import (
    GroundTruthGenerator,
    StatisticalEvaluator,
    ThresholdAnalyzer,
    ScalabilityProfiler,
    MultiMotifBenchmarkSuite
)


def run_full_integrated_pipeline() -> None:
    """Runs the complete end-to-end multi-stage pipeline connecting all 4 members' work."""
    print("=" * 80)
    print("   AUTOMATA-BASED PROBABILISTIC DNA MOTIF DETECTION - INTEGRATED PIPELINE   ")
    print("        Compiler Design & Theory of Computation Case Study Project         ")
    print("=" * 80)

    # -------------------------------------------------------------
    # STAGE 1 (Member 1): 20-Motif Selection & Information Modeling
    # -------------------------------------------------------------
    print("\n[STAGE 1 / MEMBER 1] Initializing Motif Profile & Information Content...")
    bank = MotifBank()
    all_profiles = bank.get_all_profiles()
    print(f"  [+] Loaded {len(all_profiles)} verified biological transcription factors into Motif Bank.")

    motif = bank.get("GATA1")
    print(f"  Target Transcription Factor: {motif.name} (JASPAR ID: {motif.motif_id})")
    print(f"  Organism:                    {motif.organism}")
    print(f"  Motif Length:                {motif.length} bp")
    print(f"  Strict Consensus:            {motif.consensus_sequence}")
    print(f"  IUPAC Consensus Pattern:     {motif.iupac_consensus}")
    print(f"  Total Information Content:   {motif.total_information_content:.3f} bits (Max {2*motif.length} bits)")
    print(f"  Score Range:                 [{motif.min_possible_score:.2f} .. {motif.max_possible_score:.2f}]")

    # -------------------------------------------------------------
    # STAGE 2 (Member 2): Automata Construction & Compiler Backend
    # -------------------------------------------------------------
    print("\n[STAGE 2 / MEMBER 2] Compiling Formal Automata...")
    dfa = ClassicalExactDFA(f"{motif.name}_DFA", motif.iupac_consensus)
    wfa = LinearWFA(motif)
    mismatch_wfa = MismatchAwareWFA(motif, max_mismatches=2)
    print(f"  [+] Compiled Classical Exact DFA ({len(dfa.states)} states)")
    print(f"  [+] Compiled Linear 1D-WFA ({len(wfa.states)} states, {sum(len(v) for v in wfa.transitions.values())} weighted transitions)")
    print(f"  [+] Compiled Mismatch-Aware 2D-WFA ({len(mismatch_wfa.states)} 2D states, K_max=2)")

    # -------------------------------------------------------------
    # STAGE 3 (Member 3): 20-Genome Registry & Dual-Strand Scanning
    # -------------------------------------------------------------
    print("\n[STAGE 3 / MEMBER 3] Loading Reference Genomes & Scanning Both Strands...")
    registry = GenomeRegistry()
    all_genomes = registry.get_all_records()
    print(f"  [+] Loaded {len(all_genomes)} reference genomes across human and model organisms.")

    target_genome = registry.get("GENOME_HUMAN_HBB_CHR11")
    print(f"  Selected Locus:  {target_genome.locus_name} ({target_genome.organism})")
    print(f"  Sequence Length: {target_genome.length:,} bp (GC: {target_genome.get_gc_content():.2%})")

    pipeline = DualStrandScannerPipeline(motif, max_mismatches=2)
    threshold = 4.0
    t0 = time.perf_counter()
    hits = pipeline.scan_genomic_record(target_genome, threshold=threshold, scan_both_strands=True)
    scan_time = (time.perf_counter() - t0) * 1000

    print(f"  [+] Scan Completed in {scan_time:.2f} ms")
    print(f"  [+] Total Candidate Hits Found (Score >= {threshold}): {len(hits)}")
    print(f"      - Forward (+) Strand Hits: {sum(1 for h in hits if h.strand == '+')}")
    print(f"      - Reverse (-) Strand Hits: {sum(1 for h in hits if h.strand == '-')}")
    print(f"      - Strict DFA Exact Hits:   {sum(1 for h in hits if h.is_exact_match)}")

    # -------------------------------------------------------------
    # STAGE 4 (Member 4): Ground-Truth Evaluation & Multi-Motif Suite
    # -------------------------------------------------------------
    print("\n[STAGE 4 / MEMBER 4] Ground-Truth Benchmarking & Multi-Motif Comparison...")
    syn_dataset = GroundTruthGenerator.generate_dataset(
        motif_profile=motif, sequence_length=15_000, num_sites=40, mutation_rate=0.10, seed=42
    )
    print(f"  [+] Generated Synthetic Benchmark: {syn_dataset.length:,} bp with {syn_dataset.total_true_sites} known true sites")

    eval_points = ThresholdAnalyzer.sweep_thresholds(
        dataset=syn_dataset, motif_profile=motif, min_threshold=1.0, max_threshold=8.0, num_steps=6
    )
    opt_pt = ThresholdAnalyzer.find_optimal_threshold(eval_points)
    print(f"  [+] Optimal Operating Threshold: {opt_pt.threshold:.2f} (Max F1-Score = {opt_pt.confusion_matrix.f1_score:.3f})")
    print(f"      Precision: {opt_pt.confusion_matrix.precision:.3f} | Recall: {opt_pt.confusion_matrix.recall:.3f} | MCC: {opt_pt.confusion_matrix.mcc:.3f}")

    # Comparative evaluation across 5 representative transcription factors
    print("\n  [+] Comparative Benchmark across 5 Diverse Motifs (6 to 20 bp):")
    comp_rows = MultiMotifBenchmarkSuite.evaluate_multiple_motifs(
        motif_names=["GATA1", "SP1", "AP-1", "MYC", "p53"],
        sequence_length=12_000,
        num_sites_per_motif=25
    )
    print(MultiMotifBenchmarkSuite.format_comparison_table(comp_rows))

    print("\n" + "=" * 80)
    print("   ALL 4 MODULES EXECUTED SUCCESSFULLY AND INTEGRATED COMPREHENSIVELY!   ")
    print("=" * 80)


def main() -> None:
    parser = argparse.ArgumentParser(description="Automata-Based DNA Motif Detection System")
    parser.add_argument("--demo", choices=["all", "m1", "m2", "m3", "m4"], default="all",
                        help="Select demonstration module to execute (m1=Member1, m2=Member2, m3=Member3, m4=Member4, all=Integrated)")
    parser.add_argument("--gui", action="store_true", help="Launch interactive graphical user interface (GUI)")
    args = parser.parse_args()

    if args.gui:
        import tkinter as tk
        from gui_app import MotifDetectionGUI
        print("[+] Launching interactive graphical user interface...")
        root = tk.Tk()
        app = MotifDetectionGUI(root)
        root.mainloop()
        return

    if args.demo == "m1":
        from module1_motif_models import run_member1_demonstration
        run_member1_demonstration()
    elif args.demo == "m2":
        from module2_automata_engine import run_member2_demonstration
        run_member2_demonstration()
    elif args.demo == "m3":
        from module3_scanner_pipeline import run_member3_demonstration
        run_member3_demonstration()
    elif args.demo == "m4":
        from module4_eval_benchmark import run_member4_demonstration
        run_member4_demonstration()
    else:
        run_full_integrated_pipeline()


if __name__ == "__main__":
    main()
