from __future__ import annotations

import math
import random
import time
import tracemalloc
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from module1_motif_models import DNA_BASES, BASE_TO_INDEX, MotifProfile, MotifBank
from module2_automata_engine import LinearWFA, MismatchAwareWFA, ClassicalExactDFA
from module3_scanner_pipeline import (
    GenomicSequence,
    MotifHit,
    DirectPSSMScanner,
    PrunedWFAScanner,
    DualStrandScannerPipeline,
    reverse_complement
)


# ============================================================================
# 1. GROUND-TRUTH SYNTHETIC DATASET GENERATOR
# ============================================================================

@dataclass
class GroundTruthSite:
    """Represents a known implanted true motif binding site."""
    site_id: str
    start: int
    end: int
    strand: str
    motif_id: str
    original_sequence: str
    mutated_sequence: str
    mutation_count: int
    mutation_details: List[Tuple[int, str, str]] = field(default_factory=list) # (pos, from, to)


@dataclass
class SyntheticBenchmarkDataset:
    """Container for synthetic test DNA sequences with known ground-truth annotations."""
    dataset_name: str
    sequence: str
    gc_content: float
    true_sites: List[GroundTruthSite] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def total_true_sites(self) -> int:
        return len(self.true_sites)


class GroundTruthGenerator:
    """
    Constructs high-fidelity synthetic genomic datasets with known motif implants
    to evaluate algorithmic precision and recall without circular logic.
    """
    @classmethod
    def generate_dataset(
        cls,
        motif_profile: MotifProfile,
        sequence_length: int = 20_000,
        num_sites: int = 50,
        gc_content: float = 0.45,
        mutation_rate: float = 0.10,
        seed: int = 42
    ) -> SyntheticBenchmarkDataset:
        """
        Generates a synthetic DNA sequence with implanted motifs at non-overlapping positions.
        """
        rng = random.Random(seed)
        at_prob = (1.0 - gc_content) / 2.0
        gc_prob = gc_content / 2.0
        weights = [at_prob, gc_prob, gc_prob, at_prob]

        # 1. Generate random background DNA
        background_bases = rng.choices(DNA_BASES, weights=weights, k=sequence_length)

        # 2. Select non-overlapping implantation coordinates
        L = motif_profile.length
        min_spacing = L + 10
        possible_positions = list(range(10, sequence_length - L - 10, min_spacing))
        rng.shuffle(possible_positions)

        num_to_sample = min(num_sites, len(possible_positions))
        selected_positions = sorted(possible_positions[:num_to_sample])

        true_sites: List[GroundTruthSite] = []

        # 3. Sample and implant motifs
        for idx, pos in enumerate(selected_positions):
            # Sample a sequence from the motif PPM
            sampled_chars = []
            for j in range(L):
                col_probs = [motif_profile.ppm[b][j] for b in range(4)]
                base = rng.choices(DNA_BASES, weights=col_probs)[0]
                sampled_chars.append(base)
            orig_seq = "".join(sampled_chars)

            # Apply mutation noise
            mutated_chars = list(sampled_chars)
            mutations = []
            for j in range(L):
                if rng.random() < mutation_rate:
                    # Mutate to another base
                    current_b = mutated_chars[j]
                    other_bases = [b for b in DNA_BASES if b != current_b]
                    new_b = rng.choice(other_bases)
                    mutated_chars[j] = new_b
                    mutations.append((j, current_b, new_b))

            final_motif = "".join(mutated_chars)
            strand = "+" if rng.random() >= 0.5 else "-"

            if strand == "+":
                background_bases[pos:pos + L] = list(final_motif)
            else:
                rc_motif = reverse_complement(final_motif)
                background_bases[pos:pos + L] = list(rc_motif)

            site = GroundTruthSite(
                site_id=f"site_{idx+1:03d}",
                start=pos,
                end=pos + L,
                strand=strand,
                motif_id=motif_profile.motif_id,
                original_sequence=orig_seq,
                mutated_sequence=final_motif,
                mutation_count=len(mutations),
                mutation_details=mutations
            )
            true_sites.append(site)

        full_seq = "".join(background_bases)
        return SyntheticBenchmarkDataset(
            dataset_name=f"Synthetic_{motif_profile.name}_L{sequence_length}_N{len(true_sites)}",
            sequence=full_seq,
            gc_content=gc_content,
            true_sites=true_sites
        )


# ============================================================================
# 2. STATISTICAL EVALUATION ENGINE & CONFUSION MATRIX
# ============================================================================

@dataclass
class ConfusionMatrix:
    """Complete statistical evaluation metrics for binary classification."""
    tp: int = 0  # True Positives: Predicted hits matching ground truth sites
    fp: int = 0  # False Positives: Predicted hits not matching any ground truth
    fn: int = 0  # False Negatives: Ground truth sites missed by predictions
    tn: int = 0  # True Negatives: Non-motif windows correctly rejected

    @property
    def total_predictions(self) -> int:
        return self.tp + self.fp

    @property
    def total_ground_truth(self) -> int:
        return self.tp + self.fn

    @property
    def precision(self) -> float:
        """Precision (Positive Predictive Value) = TP / (TP + FP)"""
        return self.tp / (self.tp + self.fp) if (self.tp + self.fp) > 0 else 0.0

    @property
    def recall(self) -> float:
        """Recall (Sensitivity / True Positive Rate) = TP / (TP + FN)"""
        return self.tp / (self.tp + self.fn) if (self.tp + self.fn) > 0 else 0.0

    @property
    def specificity(self) -> float:
        """Specificity (True Negative Rate) = TN / (TN + FP)"""
        return self.tn / (self.tn + self.fp) if (self.tn + self.fp) > 0 else 0.0

    @property
    def fpr(self) -> float:
        """False Positive Rate (1 - Specificity) = FP / (FP + TN)"""
        return self.fp / (self.fp + self.tn) if (self.fp + self.tn) > 0 else 0.0

    @property
    def f1_score(self) -> float:
        """F1-Score (Harmonic mean of precision and recall)"""
        p = self.precision
        r = self.recall
        return (2.0 * p * r) / (p + r) if (p + r) > 0 else 0.0

    @property
    def mcc(self) -> float:
        """Matthews Correlation Coefficient (MCC)"""
        numerator = (self.tp * self.tn) - (self.fp * self.fn)
        denominator = math.sqrt(
            (self.tp + self.fp) * (self.tp + self.fn) * (self.tn + self.fp) * (self.tn + self.fn)
        )
        return numerator / denominator if denominator > 0 else 0.0

    @property
    def youden_j(self) -> float:
        """Youden's J Statistic = Sensitivity + Specificity - 1"""
        return self.recall + self.specificity - 1.0


class StatisticalEvaluator:
    """
    Evaluates detected genomic hits against ground truth with spatial overlap matching.
    """
    @classmethod
    def evaluate_hits(
        cls,
        predicted_hits: List[MotifHit],
        ground_truth: List[GroundTruthSite],
        total_sequence_length: int,
        motif_length: int,
        tolerance_bp: int = 2
    ) -> ConfusionMatrix:
        """
        Computes TP, FP, FN, TN by matching predicted windows against ground truth sites.
        A predicted hit is a TP if |predicted_start - true_start| <= tolerance_bp and strand matches.
        """
        matched_true_indices = set()
        tp = 0
        fp = 0

        for pred in predicted_hits:
            match_found = False
            for gt_idx, gt in enumerate(ground_truth):
                if gt_idx in matched_true_indices:
                    continue
                # Spatial proximity and strand check
                if abs(pred.start - gt.start) <= tolerance_bp and pred.strand == gt.strand:
                    match_found = True
                    matched_true_indices.add(gt_idx)
                    tp += 1
                    break
            if not match_found:
                fp += 1

        fn = len(ground_truth) - len(matched_true_indices)

        # Total possible windows in sequence on both strands
        total_windows = 2 * (total_sequence_length - motif_length + 1)
        # Non-motif windows
        non_motif_windows = max(0, total_windows - len(ground_truth))
        tn = max(0, non_motif_windows - fp)

        return ConfusionMatrix(tp=tp, fp=fp, fn=fn, tn=tn)


# ============================================================================
# 3. THRESHOLD ANALYSIS & ROC / PR CURVE ENGINE
# ============================================================================

@dataclass
class ThresholdPoint:
    """Performance metrics evaluated at a specific score threshold."""
    threshold: float
    confusion_matrix: ConfusionMatrix
    p_value: float


class ThresholdAnalyzer:
    """
    Performs comprehensive sensitivity/specificity sweeps across candidate score thresholds
    to generate ROC and Precision-Recall curves and locate optimal operational cutoffs.
    """
    @classmethod
    def sweep_thresholds(
        cls,
        dataset: SyntheticBenchmarkDataset,
        motif_profile: MotifProfile,
        min_threshold: Optional[float] = None,
        max_threshold: Optional[float] = None,
        num_steps: int = 30
    ) -> List[ThresholdPoint]:
        """Sweeps score thresholds and records full performance curve."""
        pipeline = DualStrandScannerPipeline(motif_profile)
        rec = GenomicSequence(
            sequence_id=dataset.dataset_name,
            description="Synthetic Test",
            sequence=dataset.sequence
        )

        min_t = min_threshold if min_threshold is not None else max(0.0, motif_profile.min_possible_score)
        max_t = max_threshold if max_threshold is not None else motif_profile.max_possible_score - 0.5
        step_size = (max_t - min_t) / max(1, num_steps - 1)

        points: List[ThresholdPoint] = []
        for step in range(num_steps):
            t = min_t + step * step_size
            hits = pipeline.scan_genomic_record(rec, threshold=t, scan_both_strands=True)
            cm = StatisticalEvaluator.evaluate_hits(
                predicted_hits=hits,
                ground_truth=dataset.true_sites,
                total_sequence_length=dataset.length,
                motif_length=motif_profile.length
            )
            pval = motif_profile.score_to_pvalue(t)
            points.append(ThresholdPoint(threshold=t, confusion_matrix=cm, p_value=pval))

        return points

    @classmethod
    def compute_auc_roc(cls, points: List[ThresholdPoint]) -> float:
        """Calculates Area Under ROC Curve (AUC-ROC) using trapezoidal numerical integration."""
        # Sort points by FPR ascending
        sorted_pts = sorted(points, key=lambda p: (p.confusion_matrix.fpr, p.confusion_matrix.recall))
        auc = 0.0
        for i in range(len(sorted_pts) - 1):
            x1 = sorted_pts[i].confusion_matrix.fpr
            x2 = sorted_pts[i+1].confusion_matrix.fpr
            y1 = sorted_pts[i].confusion_matrix.recall
            y2 = sorted_pts[i+1].confusion_matrix.recall
            auc += (x2 - x1) * (y1 + y2) / 2.0
        return max(0.0, min(1.0, auc))

    @classmethod
    def find_optimal_threshold(cls, points: List[ThresholdPoint]) -> ThresholdPoint:
        """Finds operating threshold that maximizes F1-Score."""
        return max(points, key=lambda p: p.confusion_matrix.f1_score)


# ============================================================================
# 4. PERFORMANCE & SCALABILITY PROFILER SUITE (UPGRADE 8)
# ============================================================================

@dataclass
class BenchmarkResult:
    """Performance metrics for an individual engine benchmark run."""
    engine_name: str
    sequence_length: int
    num_hits_found: int
    wall_time_ms: float
    cpu_time_ms: float
    peak_memory_kb: float
    throughput_mbps: float  # Megabases per second


class ScalabilityProfiler:
    """
    Measures algorithmic runtime scaling, CPU utilization, peak memory allocation,
    and throughput across varying sequence scales (10KB to 1MB).
    """
    @classmethod
    def profile_scanner(
        cls,
        scanner_fn: Callable[[str, float], List[Any]],
        sequence: str,
        threshold: float,
        engine_name: str,
        warmup_runs: int = 0,
        benchmark_runs: int = 1
    ) -> BenchmarkResult:
        """Executes high-precision timing and memory profiling of a scanner function."""
        # Warmup
        for _ in range(warmup_runs):
            _ = scanner_fn(sequence, threshold)

        tracemalloc.start()
        start_wall = time.perf_counter()
        start_cpu = time.process_time()

        hits = []
        for _ in range(benchmark_runs):
            hits = scanner_fn(sequence, threshold)

        end_wall = time.perf_counter()
        end_cpu = time.process_time()
        _, peak_mem_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        avg_wall_s = (end_wall - start_wall) / benchmark_runs
        avg_cpu_s = (end_cpu - start_cpu) / benchmark_runs
        wall_ms = avg_wall_s * 1000.0
        cpu_ms = avg_cpu_s * 1000.0
        peak_kb = peak_mem_bytes / 1024.0

        # Throughput in Megabases / sec
        mb_scanned = len(sequence) / 1_000_000.0
        throughput = mb_scanned / avg_wall_s if avg_wall_s > 0 else 0.0

        return BenchmarkResult(
            engine_name=engine_name,
            sequence_length=len(sequence),
            num_hits_found=len(hits),
            wall_time_ms=wall_ms,
            cpu_time_ms=cpu_ms,
            peak_memory_kb=peak_kb,
            throughput_mbps=throughput
        )

    @classmethod
    def run_multi_engine_scalability_suite(
        cls,
        motif_profile: MotifProfile,
        lengths: List[int] = [10_000, 50_000, 100_000, 500_000]
    ) -> List[BenchmarkResult]:
        """Runs comparative scalability benchmarks across all five compiler engines."""
        results: List[BenchmarkResult] = []
        pipeline = DualStrandScannerPipeline(motif_profile)
        threshold = 3.5

        for length in lengths:
            # Generate deterministic benchmark sequence
            rng = random.Random(length)
            seq = "".join(rng.choices(DNA_BASES, k=length))

            # 1. Classical DFA Baseline
            exact_dfa = ClassicalExactDFA(f"{motif_profile.name}_Exact", motif_profile.iupac_consensus)
            def run_dfa(s: str, t: float) -> List[int]:
                hits = []
                L = exact_dfa.length
                for i in range(len(s) - L + 1):
                    if exact_dfa.match_window(s[i:i+L]):
                        hits.append(i)
                return hits

            res_dfa = cls.profile_scanner(run_dfa, seq, threshold, "Classical Exact DFA")
            results.append(res_dfa)

            # 2. Direct PSSM Scanner
            def run_direct(s: str, t: float) -> List[Any]:
                return pipeline.direct_scanner.scan_sequence(s, threshold=t)
            res_direct = cls.profile_scanner(run_direct, seq, threshold, "Direct PSSM Scanner")
            results.append(res_direct)

            # 3. Lookahead Pruned WFA Scanner
            def run_pruned(s: str, t: float) -> List[Any]:
                return pipeline.pruned_scanner.scan_sequence(s, threshold=t)
            res_pruned = cls.profile_scanner(run_pruned, seq, threshold, "Pruned Lookahead WFA")
            results.append(res_pruned)

            # 4. 2D Mismatch-Aware WFA (K=2)
            def run_2d(s: str, t: float) -> List[Any]:
                hits = []
                L = motif_profile.length
                for i in range(len(s) - L + 1):
                    score, k = pipeline.mismatch_wfa.score_window(s[i:i+L])
                    if score >= t and k >= 0:
                        hits.append((i, score, k))
                return hits
            res_2d = cls.profile_scanner(run_2d, seq, threshold, "2D Mismatch-Aware WFA")
            results.append(res_2d)

        return results


# ============================================================================
# 5. MULTI-MOTIF COMPARATIVE BENCHMARK SUITE & REPORT ENGINE
# ============================================================================

@dataclass
class MultiMotifComparisonRow:
    """Stores comparative evaluation metrics for a transcription factor."""
    motif_name: str
    motif_id: str
    motif_length: int
    optimal_threshold: float
    true_positives: int
    false_positives: int
    precision: float
    recall: float
    f1_score: float
    mcc: float
    auc_roc: float


class MultiMotifBenchmarkSuite:
    """
    Executes automated comparative benchmarking across multiple transcription factors
    from the 20-motif bank, measuring performance across variable motif lengths and information content.
    """
    @classmethod
    def evaluate_multiple_motifs(
        cls,
        motif_names: List[str] = ["GATA1", "SP1", "AP-1", "MYC", "p53"],
        sequence_length: int = 15_000,
        num_sites_per_motif: int = 30,
        mutation_rate: float = 0.10
    ) -> List[MultiMotifComparisonRow]:
        """Runs threshold sweeps and evaluation across multiple target motifs."""
        bank = MotifBank()
        rows: List[MultiMotifComparisonRow] = []

        for name in motif_names:
            profile = bank.get(name)
            dataset = GroundTruthGenerator.generate_dataset(
                motif_profile=profile,
                sequence_length=sequence_length,
                num_sites=num_sites_per_motif,
                mutation_rate=mutation_rate,
                seed=42 + profile.length
            )

            # Sweep 6 thresholds
            pts = ThresholdAnalyzer.sweep_thresholds(
                dataset=dataset,
                motif_profile=profile,
                num_steps=6
            )
            opt_pt = ThresholdAnalyzer.find_optimal_threshold(pts)
            auc = ThresholdAnalyzer.compute_auc_roc(pts)

            cm = opt_pt.confusion_matrix
            row = MultiMotifComparisonRow(
                motif_name=profile.name,
                motif_id=profile.motif_id,
                motif_length=profile.length,
                optimal_threshold=opt_pt.threshold,
                true_positives=cm.tp,
                false_positives=cm.fp,
                precision=cm.precision,
                recall=cm.recall,
                f1_score=cm.f1_score,
                mcc=cm.mcc,
                auc_roc=auc
            )
            rows.append(row)

        return rows

    @staticmethod
    def format_comparison_table(rows: List[MultiMotifComparisonRow]) -> str:
        """Formats multi-motif comparison metrics as an aligned ASCII table."""
        headers = ["Motif", "JASPAR ID", "Len", "Opt Thresh", "TP", "FP", "Precision", "Recall", "F1-Score", "MCC", "AUC-ROC"]
        lines = [
            "-" * 95,
            f"{headers[0]:<10} {headers[1]:<11} {headers[2]:<5} {headers[3]:<12} {headers[4]:<5} {headers[5]:<5} {headers[6]:<10} {headers[7]:<8} {headers[8]:<9} {headers[9]:<7} {headers[10]:<8}",
            "-" * 95
        ]
        for r in rows:
            lines.append(
                f"{r.motif_name:<10} {r.motif_id:<11} {r.motif_length:<5} {r.optimal_threshold:<12.2f} "
                f"{r.true_positives:<5} {r.false_positives:<5} {r.precision:<10.3f} {r.recall:<8.3f} "
                f"{r.f1_score:<9.3f} {r.mcc:<7.3f} {r.auc_roc:<8.3f}"
            )
        lines.append("-" * 95)
        return "\n".join(lines)


# ============================================================================
# STANDALONE VERIFICATION & DEMONSTRATION SUITE
# ============================================================================

def run_member4_demonstration() -> None:
    """Executes a complete demonstration of all Member 4 functionality."""
    print("=" * 70)
    print("MEMBER 4 MODULE DEMO: Evaluation Framework & Scalability Profiler")
    print("=" * 70)

    # 1. Load GATA1 Motif
    bank = MotifBank()
    gata1 = bank.get("GATA1")

    # 2. Generate Ground-Truth Synthetic Dataset
    print("\n[+] Generating Synthetic Dataset with Implanted Ground-Truth Motifs:")
    dataset = GroundTruthGenerator.generate_dataset(
        motif_profile=gata1,
        sequence_length=15_000,
        num_sites=40,
        gc_content=0.45,
        mutation_rate=0.10,
        seed=101
    )
    print(f"    Dataset:           {dataset.dataset_name}")
    print(f"    Total Length:      {dataset.length:,} bp")
    print(f"    True Sites Impl.:  {dataset.total_true_sites}")

    # 3. Systematic Threshold Sweep & ROC / PR Analysis
    print("\n[+] Executing Threshold Sweep & Statistical Evaluation:")
    points = ThresholdAnalyzer.sweep_thresholds(
        dataset=dataset,
        motif_profile=gata1,
        min_threshold=1.0,
        max_threshold=8.0,
        num_steps=8
    )

    print(f"{'Threshold':<11} {'TP':<5} {'FP':<6} {'FN':<5} {'Precision':<11} {'Recall':<9} {'F1-Score':<10} {'MCC':<8}")
    print("-" * 75)
    for p in points:
        cm = p.confusion_matrix
        print(f"{p.threshold:<11.2f} {cm.tp:<5} {cm.fp:<6} {cm.fn:<5} "
              f"{cm.precision:<11.3f} {cm.recall:<9.3f} {cm.f1_score:<10.3f} {cm.mcc:<8.3f}")

    optimal_pt = ThresholdAnalyzer.find_optimal_threshold(points)
    print(f"\n[+] Optimal Operating Threshold: {optimal_pt.threshold:.2f} (Max F1 = {optimal_pt.confusion_matrix.f1_score:.3f})")

    # 4. Scalability & Performance Benchmarking Suite
    print("\n[+] Running Performance & Scalability Suite across Sequence Scales:")
    bench_results = ScalabilityProfiler.run_multi_engine_scalability_suite(
        motif_profile=gata1,
        lengths=[5_000, 10_000]
    )

    print(f"\n{'Engine':<26} {'Length':<10} {'Hits':<8} {'Time (ms)':<12} {'Peak Mem':<12} {'Throughput'}")
    print("-" * 80)
    for r in bench_results:
        print(f"{r.engine_name:<26} {r.sequence_length:<10,d} {r.num_hits_found:<8} "
              f"{r.wall_time_ms:<12.3f} {r.peak_memory_kb:<10.2f} KB  {r.throughput_mbps:.2f} Mbp/s")

    print("=" * 70)


if __name__ == "__main__":
    run_member4_demonstration()
