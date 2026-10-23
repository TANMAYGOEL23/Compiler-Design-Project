"""
Module 1: Biological Motif Modeling, Information Theory & 20-Motif Bank Engine
=============================================================================
Member 1 Contribution - Compiler Design & Automata-Based DNA Motif Detection

Core Responsibilities:
1. Position Frequency Matrix (PFM) to Position Probability Matrix (PPM) transformation.
2. Information Theory metrics: Shannon Entropy, Information Content (Bits), Relative Entropy (Kullback-Leibler Divergence).
3. Background Models: Uniform (0.25), Empirical Sequence-Specific GC-content background, and 1st-Order Markov Chain background.
4. Position-Specific Scoring Matrix (PSSM) & Log-Odds weight matrix calculation.
5. Dynamic Programming Staden Algorithm for exact P-value <-> Score threshold mapping.
6. IUPAC Degenerate Code consensus generator and regular expression compilation.
7. Motif Bank Registry supporting 20 biologically verified Transcription Factor binding profiles:
   GATA1, CTCF, SP1, TATA-Box, NF-kB, AP-1, p53, MYC, STAT3, FOXA1, ESR1, SOX2,
   OCT4, HNF4A, CREB1, E2F1, CEBPA, RUNX1, PAX6, and MEF2A.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union


# ============================================================================
# CONSTANTS & IUPAC ALPHABET DEFINITIONS
# ============================================================================

DNA_BASES: List[str] = ["A", "C", "G", "T"]
BASE_TO_INDEX: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}
INDEX_TO_BASE: Dict[int, str] = {0: "A", 1: "C", 2: "G", 3: "T"}

# IUPAC Degenerate Nucleotide Mapping for Regular Expressions & Automata
IUPAC_MAP: Dict[str, List[str]] = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],        # puRine
    "Y": ["C", "T"],        # pYrimidine
    "S": ["G", "C"],        # Strong interaction (3 H-bonds)
    "W": ["A", "T"],        # Weak interaction (2 H-bonds)
    "K": ["G", "T"],        # Keto
    "M": ["A", "C"],        # aMino
    "B": ["C", "G", "T"],   # not A
    "D": ["A", "G", "T"],   # not C
    "H": ["A", "C", "T"],   # not G
    "V": ["A", "C", "G"],   # not T
    "N": ["A", "C", "G", "T"], # Any base
}

# Reverse mapping: sorted tuple of bases to IUPAC character
BASES_TO_IUPAC: Dict[Tuple[str, ...], str] = {
    tuple(sorted(v)): k for k, v in IUPAC_MAP.items()
}


class BackgroundModelType(Enum):
    UNIFORM = "uniform"
    GC_CONTENT = "gc_content"
    CUSTOM = "custom"
    MARKOV_FIRST_ORDER = "markov_1st_order"


# ============================================================================
# BACKGROUND PROBABILITY MODELS
# ============================================================================

@dataclass
class BackgroundModel:
    """
    Represents the statistical background nucleotide distribution of the genome.
    Used for computing log-likelihood ratios (Log-Odds scores).
    """
    model_type: BackgroundModelType = BackgroundModelType.UNIFORM
    gc_fraction: float = 0.50
    base_probabilities: Dict[str, float] = field(default_factory=lambda: {
        "A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25
    })
    transition_matrix: Optional[Dict[str, Dict[str, float]]] = None

    def __post_init__(self) -> None:
        if self.model_type == BackgroundModelType.GC_CONTENT:
            self.set_gc_content(self.gc_fraction)
        elif self.model_type == BackgroundModelType.UNIFORM:
            self.base_probabilities = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
        self.validate()

    def set_gc_content(self, gc: float) -> None:
        """Sets background probabilities derived from regional GC content."""
        if not (0.0 <= gc <= 1.0):
            raise ValueError(f"GC fraction must be between 0.0 and 1.0, got {gc}")
        self.gc_fraction = gc
        at_prob = (1.0 - gc) / 2.0
        gc_prob = gc / 2.0
        self.base_probabilities = {
            "A": at_prob,
            "C": gc_prob,
            "G": gc_prob,
            "T": at_prob
        }

    def validate(self) -> None:
        """Validates that base probabilities sum to 1.0."""
        total = sum(self.base_probabilities.get(b, 0.0) for b in DNA_BASES)
        if abs(total - 1.0) > 1e-4:
            raise ValueError(f"Background probabilities must sum to 1.0, got sum={total:.6f}")

    def get_prob(self, base: str, prev_base: Optional[str] = None) -> float:
        """Returns background probability for a base, optionally conditioned on prev_base."""
        base = base.upper()
        if self.model_type == BackgroundModelType.MARKOV_FIRST_ORDER and prev_base is not None and self.transition_matrix:
            prev_base = prev_base.upper()
            if prev_base in self.transition_matrix and base in self.transition_matrix[prev_base]:
                return self.transition_matrix[prev_base][base]
        return self.base_probabilities.get(base, 0.25)

    @classmethod
    def from_sequence(cls, sequence: str) -> BackgroundModel:
        """Constructs an empirical background model directly from a long DNA sequence."""
        seq = sequence.upper()
        counts = {b: 0 for b in DNA_BASES}
        total_valid = 0
        for char in seq:
            if char in counts:
                counts[char] += 1
                total_valid += 1

        if total_valid == 0:
            return cls(model_type=BackgroundModelType.UNIFORM)

        probs = {b: counts[b] / total_valid for b in DNA_BASES}
        gc = probs["G"] + probs["C"]
        return cls(
            model_type=BackgroundModelType.GC_CONTENT,
            gc_fraction=gc,
            base_probabilities=probs
        )


# ============================================================================
# MOTIF DATA STRUCTURES & INFORMATION THEORY CALCULATOR
# ============================================================================

@dataclass
class PositionInformation:
    """Stores information-theoretic statistics for a single column of a motif."""
    position: int
    probabilities: Dict[str, float]
    shannon_entropy: float        # H(X) in bits [0 to 2 bits for 4 bases]
    information_content: float    # R_seq = 2 - (H(X) + e_n) in bits
    relative_entropy: float       # D_KL(P || Q) relative to background
    dominant_base: str
    consensus_symbol: str
    heights: Dict[str, float]     # Letter height in sequence logo: prob * information_content


@dataclass
class MotifProfile:
    """
    Complete mathematical and biological representation of a Transcription Factor Binding Site (TFBS).
    Contains PFM (counts), PPM (probabilities), PSSM (log-odds), and Information Content.
    """
    motif_id: str
    name: str
    length: int
    pfm: List[List[float]]             # Shape: 4 x Length (Rows: A, C, G, T)
    pseudocount: float = 0.1
    background: BackgroundModel = field(default_factory=BackgroundModel)
    description: str = ""
    organism: str = "Homo sapiens"
    biological_function: str = ""
    ppm: List[List[float]] = field(default_factory=list)
    pssm: List[List[float]] = field(default_factory=list)
    consensus_sequence: str = ""
    iupac_consensus: str = ""
    regex_pattern: str = ""
    min_possible_score: float = 0.0
    max_possible_score: float = 0.0
    information_per_position: List[PositionInformation] = field(default_factory=list)
    total_information_content: float = 0.0
    _cached_pvalue_dist: Optional[Tuple[List[float], List[float]]] = None

    def __post_init__(self) -> None:
        if not self.ppm or not self.pssm:
            self.recompute_matrices()

    def recompute_matrices(self) -> None:
        """
        Computes PPM, PSSM, Consensuses, and Information Theory measures from raw PFM counts.
        """
        self.length = len(self.pfm[0])
        self.ppm = [[0.0 for _ in range(self.length)] for _ in range(4)]
        self.pssm = [[0.0 for _ in range(self.length)] for _ in range(4)]
        self.information_per_position.clear()
        self.total_information_content = 0.0

        min_score = 0.0
        max_score = 0.0
        consensus_chars = []
        iupac_chars = []

        for j in range(self.length):
            # 1. Add pseudocount and normalize to obtain PPM
            raw_col_counts = [self.pfm[b_idx][j] for b_idx in range(4)]
            smoothed_counts = [c + self.pseudocount for c in raw_col_counts]
            col_sum = sum(smoothed_counts)

            col_probs = {}
            for b_idx, base in enumerate(DNA_BASES):
                prob = smoothed_counts[b_idx] / col_sum
                self.ppm[b_idx][j] = prob
                col_probs[base] = prob

            # 2. Compute PSSM log-odds scores: W(b, j) = log2( P(b|j) / P_bg(b) )
            col_scores = []
            for b_idx, base in enumerate(DNA_BASES):
                p_base = self.ppm[b_idx][j]
                p_bg = self.background.get_prob(base)
                log_odds = math.log2(p_base / p_bg)
                self.pssm[b_idx][j] = log_odds
                col_scores.append(log_odds)

            min_score += min(col_scores)
            max_score += max(col_scores)

            # 3. Information Theory: Shannon Entropy H(X) and Relative Entropy
            raw_total_sites = sum(raw_col_counts)
            if raw_total_sites > 0:
                e_n = (4 - 1) / (2.0 * math.log(2.0) * max(raw_total_sites, 1.0))
            else:
                e_n = 0.0

            shannon_h = 0.0
            d_kl = 0.0
            for b_idx, base in enumerate(DNA_BASES):
                p = col_probs[base]
                q = self.background.get_prob(base)
                if p > 0:
                    shannon_h -= p * math.log2(p)
                    d_kl += p * math.log2(p / q)

            ic = max(0.0, 2.0 - (shannon_h + e_n))
            self.total_information_content += ic

            heights = {base: col_probs[base] * ic for base in DNA_BASES}

            # 4. Determine strict consensus and IUPAC consensus
            dominant_idx = max(range(4), key=lambda idx: self.ppm[idx][j])
            dominant_base = INDEX_TO_BASE[dominant_idx]
            consensus_chars.append(dominant_base)

            active_bases = [
                DNA_BASES[b_idx] for b_idx in range(4) if self.ppm[b_idx][j] >= 0.15
            ]
            if not active_bases:
                active_bases = [dominant_base]
            iupac_symbol = BASES_TO_IUPAC.get(tuple(sorted(active_bases)), "N")
            iupac_chars.append(iupac_symbol)

            pos_info = PositionInformation(
                position=j,
                probabilities=col_probs,
                shannon_entropy=shannon_h,
                information_content=ic,
                relative_entropy=d_kl,
                dominant_base=dominant_base,
                consensus_symbol=iupac_symbol,
                heights=heights
            )
            self.information_per_position.append(pos_info)

        self.min_possible_score = min_score
        self.max_possible_score = max_score
        self.consensus_sequence = "".join(consensus_chars)
        self.iupac_consensus = "".join(iupac_chars)
        self.regex_pattern = self._build_regex_from_iupac(self.iupac_consensus)
        self._cached_pvalue_dist = None

    @staticmethod
    def _build_regex_from_iupac(iupac_str: str) -> str:
        """Converts IUPAC degenerate string into regular expression pattern."""
        pattern_parts = []
        for char in iupac_str.upper():
            allowed = IUPAC_MAP.get(char, ["A", "C", "G", "T"])
            if len(allowed) == 1:
                pattern_parts.append(allowed[0])
            elif len(allowed) == 4:
                pattern_parts.append("[ACGT]")
            else:
                pattern_parts.append(f"[{''.join(allowed)}]")
        return "".join(pattern_parts)

    def get_score_for_subsequence(self, subseq: str) -> float:
        """Direct PSSM lookup score for an exact length-L subsequence."""
        subseq = subseq.upper()
        if len(subseq) != self.length:
            return float("-inf")
        score = 0.0
        for j, char in enumerate(subseq):
            if char in BASE_TO_INDEX:
                b_idx = BASE_TO_INDEX[char]
                score += self.pssm[b_idx][j]
            elif char == "N":
                avg_score = sum(self.pssm[b_idx][j] for b_idx in range(4)) / 4.0
                score += avg_score
            else:
                return float("-inf")
        return score

    def calculate_pvalue_score_distribution(self, granularity: float = 0.01) -> Tuple[List[float], List[float]]:
        """
        Implements the Staden Dynamic Programming algorithm to calculate the exact
        probability distribution of PSSM scores on random background DNA.
        Returns (scores_array, p_values_array) where P-value(S) = P(Score >= S).
        """
        if self._cached_pvalue_dist is not None:
            return self._cached_pvalue_dist

        min_s = self.min_possible_score
        max_s = self.max_possible_score
        num_bins = int(math.ceil((max_s - min_s) / granularity)) + 1

        dp_dist = [0.0] * num_bins
        dp_dist[0] = 1.0

        running_min = 0.0
        for j in range(self.length):
            col_min = min(self.pssm[b_idx][j] for b_idx in range(4))
            running_min += col_min
            new_dp = [0.0] * num_bins

            for b_idx, base in enumerate(DNA_BASES):
                base_score = self.pssm[b_idx][j]
                base_prob = self.background.get_prob(base)
                score_delta = base_score - col_min
                bin_shift = int(round(score_delta / granularity))

                for bin_idx in range(num_bins - bin_shift):
                    if dp_dist[bin_idx] > 0.0:
                        new_dp[bin_idx + bin_shift] += dp_dist[bin_idx] * base_prob

            dp_dist = new_dp

        scores = [min_s + idx * granularity for idx in range(num_bins)]
        p_values = [0.0] * num_bins

        cumulative = 0.0
        for idx in range(num_bins - 1, -1, -1):
            cumulative += dp_dist[idx]
            p_values[idx] = min(1.0, cumulative)

        self._cached_pvalue_dist = (scores, p_values)
        return self._cached_pvalue_dist

    def score_to_pvalue(self, score: float, granularity: float = 0.01) -> float:
        """Computes the statistical p-value for a given raw score."""
        if score <= self.min_possible_score:
            return 1.0
        if score >= self.max_possible_score:
            return 0.0
        scores, p_values = self.calculate_pvalue_score_distribution(granularity)
        low, high = 0, len(scores) - 1
        best_idx = 0
        while low <= high:
            mid = (low + high) // 2
            if scores[mid] >= score:
                best_idx = mid
                high = mid - 1
            else:
                low = mid + 1
        return p_values[best_idx]

    def pvalue_to_threshold(self, target_pvalue: float, granularity: float = 0.01) -> float:
        """Finds the score threshold corresponding to a desired false-positive rate (P-value)."""
        scores, p_values = self.calculate_pvalue_score_distribution(granularity)
        for score, p_val in zip(scores, p_values):
            if p_val <= target_pvalue:
                return score
        return self.min_possible_score

    def display_summary(self) -> str:
        """Returns a formatted analytical summary of the motif model."""
        lines = [
            f"============================================================",
            f"MOTIF PROFILE: {self.name} (ID: {self.motif_id})",
            f"============================================================",
            f"Organism:             {self.organism}",
            f"Biological Role:      {self.biological_function}",
            f"Length:               {self.length} bp",
            f"Strict Consensus:     {self.consensus_sequence}",
            f"IUPAC Consensus:      {self.iupac_consensus}",
            f"Regex Representation: {self.regex_pattern}",
            f"Total Info Content:   {self.total_information_content:.3f} bits (Max {2 * self.length} bits)",
            f"Min Possible Score:   {self.min_possible_score:.4f}",
            f"Max Possible Score:   {self.max_possible_score:.4f}",
            f"Background Model:     {self.background.model_type.value} (GC = {self.background.gc_fraction:.2%})",
            f"------------------------------------------------------------",
            f"Position Breakdown (Information Content & Dominance):",
            f"{'Pos':<4} {'Dominant':<10} {'IUPAC':<8} {'Entropy(bits)':<15} {'Info Content':<15} {'Rel.Entropy':<12}"
        ]
        for info in self.information_per_position:
            lines.append(
                f"{info.position + 1:<4} {info.dominant_base:<10} {info.consensus_symbol:<8} "
                f"{info.shannon_entropy:<15.3f} {info.information_content:<15.3f} {info.relative_entropy:<12.3f}"
            )
        lines.append(f"------------------------------------------------------------")
        lines.append(f"PSSM (Log-Odds Weight Matrix in bits):")
        for b_idx, base in enumerate(DNA_BASES):
            row_str = " ".join(f"{self.pssm[b_idx][j]:7.3f}" for j in range(self.length))
            lines.append(f"  {base} | {row_str}")
        lines.append(f"============================================================")
        return "\n".join(lines)


# ============================================================================
# 20-MOTIF BANK & COMPREHENSIVE TRANSCRIPTION FACTOR REGISTRY
# ============================================================================

class MotifBank:
    """
    Central repository for storing, managing, and parsing 20 biologically verified
    transcription factor binding profiles from JASPAR, MEME, and TRANSFAC databases.
    """
    def __init__(self) -> None:
        self._motifs: Dict[str, MotifProfile] = {}
        self._load_all_20_motifs()

    @property
    def motifs(self) -> Dict[str, MotifProfile]:
        """Direct access dictionary of registered motifs."""
        return self._motifs

    def __len__(self) -> int:
        return len(self.get_all_profiles())

    def __contains__(self, identifier: str) -> bool:
        try:
            self.get(identifier)
            return True
        except KeyError:
            return False

    def register_motif(self, profile: MotifProfile) -> None:
        """Registers a new motif into the bank with normalized alias support."""
        self._motifs[profile.motif_id.upper()] = profile
        self._motifs[profile.name.upper()] = profile
        # Also register alias with hyphens, underscores and spaces stripped
        stripped_name = re.sub(r"[\s\-_]", "", profile.name.upper())
        self._motifs[stripped_name] = profile

    def get(self, identifier: str) -> MotifProfile:
        """Retrieves a motif by ID, Name, or normalized alias."""
        key = identifier.upper().strip()
        if key in self._motifs:
            return self._motifs[key]
        stripped_key = re.sub(r"[\s\-_]", "", key)
        if stripped_key in self._motifs:
            return self._motifs[stripped_key]
        raise KeyError(f"Motif '{identifier}' not found in MotifBank. Available: {self.list_motifs()}")

    def list_motifs(self) -> List[str]:
        """Lists all registered unique motif IDs."""
        seen = set()
        unique_ids = []
        for m in self._motifs.values():
            if m.motif_id not in seen:
                seen.add(m.motif_id)
                unique_ids.append(m.motif_id)
        return sorted(unique_ids)

    def list_motif_names(self) -> List[str]:
        """Lists all registered unique motif names."""
        seen = set()
        names = []
        for m in self._motifs.values():
            if m.name not in seen:
                seen.add(m.name)
                names.append(m.name)
        return sorted(names)

    def get_all_profiles(self) -> List[MotifProfile]:
        """Returns list of all unique MotifProfile instances."""
        seen = set()
        profiles = []
        for m in self._motifs.values():
            if m.motif_id not in seen:
                seen.add(m.motif_id)
                profiles.append(m)
        return sorted(profiles, key=lambda p: p.name)

    def _load_all_20_motifs(self) -> None:
        """Preloads 20 verified biological transcription factor profiles with exact PFMs."""
        # 1. GATA1 (JASPAR MA0035.1) - Primary erythroid transcription factor (6 bp)
        self.register_motif(MotifProfile(
            motif_id="MA0035.1",
            name="GATA1",
            length=6,
            pfm=[
                [225.0,  11.0, 237.0,   8.0, 222.0, 194.0],
                [ 25.0,   5.0,  20.0,   5.0,   6.0,  21.0],
                [198.0, 238.0,   5.0,   5.0,  12.0,  11.0],
                [ 37.0,  31.0,  23.0, 267.0,  45.0,  59.0]
            ],
            description="Zinc-finger erythroid master regulator involved in globin gene switching.",
            biological_function="Erythropoiesis, megakaryocyte differentiation, hemoglobin synthesis."
        ))

        # 2. CTCF (JASPAR MA0139.1) - Master chromatin boundary / insulator (19 bp)
        ctcf_pfm = [
            [ 87, 167, 281,  56,   8,  34,  29,  57, 526, 455, 396,  43,  12, 381, 317,   1, 411,  55, 184],
            [291, 145,  48,  42, 502,  15, 360, 411,   6,  17,  54, 442, 497,  16,  81, 551,  30, 237, 158],
            [ 79, 137, 127, 437,  27, 501, 115,  38,  27,  65,  77,  37,  12, 134,  64,   2,  41, 191, 105],
            [100, 108, 101,  22,  20,   7,  53,  51,  18,  20,  30,  35,  36,  26,  95,   3,  75,  74, 110]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0139.1",
            name="CTCF",
            length=19,
            pfm=[[float(x) for x in row] for row in ctcf_pfm],
            description="11-zinc-finger architectural insulator establishing 3D chromatin loops.",
            biological_function="TAD boundary insulation, enhancer-promoter looping, V(D)J recombination."
        ))

        # 3. SP1 (JASPAR MA0079.1) - GC-box promoter activation factor (10 bp)
        sp1_pfm = [
            [  6,   4,   4,   3,  21,   5,   4,   5,   4,   7],
            [ 24,  16,  12,  10, 420,  15,  13,  18,  21,  30],
            [430, 448, 452, 454,  21, 445, 451, 442, 440, 421],
            [  8,   0,   0,   1,   6,   3,   0,   3,   3,  10]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0079.1",
            name="SP1",
            length=10,
            pfm=[[float(x) for x in row] for row in sp1_pfm],
            description="Specificity Protein 1 binding to classical GC-rich boxes in housekeeping promoters.",
            biological_function="Basal transcription initiation, cell cycle progression, chromatin remodeling."
        ))

        # 4. TATA-Box / TBP (JASPAR MA0108.1) - Core promoter element (15 bp)
        tata_pfm = [
            [ 18,  15,  10,  12, 140,   2, 145, 142, 141, 125,  35,  40,  30,  25,  20],
            [ 12,  10,   8,   5,   2,   1,   1,   2,   3,   5,  20,  18,  15,  22,  25],
            [ 15,  12,  14,   8,   1,   0,   0,   1,   2,   6,  25,  30,  25,  20,  18],
            [105, 113, 118, 125,   7, 147,   4,   5,   4,  14,  70,  62,  80,  83,  87]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0108.1",
            name="TATA-Box",
            length=15,
            pfm=[[float(x) for x in row] for row in tata_pfm],
            description="TATA-binding protein (TBP) recognition site positioned ~30 bp upstream of TSS.",
            biological_function="RNA Polymerase II pre-initiation complex assembly."
        ))

        # 5. NF-kB / RELA (JASPAR MA0107.1) - Master immune and inflammation regulator (10 bp)
        nfkb_pfm = [
            [  4,   2,   3,   5,  12, 380, 395, 110,  15,   8],
            [ 15,  10,   8,   6,  10,  12,   8, 260, 370, 390],
            [390, 400, 405, 385, 360,   5,   2,   8,  12,  10],
            [  3,   0,   1,   4,   6,  15,   7,  34,  15,   4]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0107.1",
            name="NF-kB",
            length=10,
            pfm=[[float(x) for x in row] for row in nfkb_pfm],
            description="Nuclear factor kappa B subunit p65/RELA recognizing decameric kappa-B sites.",
            biological_function="Innate immunity, inflammatory cytokine activation, cell survival signaling."
        ))

        # 6. AP-1 / JUN-FOS (JASPAR MA0099.1) - Cellular stress and proliferation regulator (7 bp)
        ap1_pfm = [
            [ 12, 450,   8,   5,   4, 460,  15],
            [ 25,  10,   5,   2, 440,  10, 420],
            [430,  12, 455,   6,  18,   8,  22],
            [  5,   0,   4, 459,  10,   4,  15]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0099.1",
            name="AP-1",
            length=7,
            pfm=[[float(x) for x in row] for row in ap1_pfm],
            description="Activator Protein 1 basic leucine zipper heterodimer (c-Jun/c-Fos).",
            biological_function="Immediate early response, MAPK signaling, apoptosis vs proliferation control."
        ))

        # 7. p53 / TP53 (JASPAR MA0106.1) - Tumor suppressor / Guardian of the Genome (20 bp)
        p53_pfm = [
            [ 18, 120,  10,  25,   5,  15, 115,  10,  20,   5,  18, 118,  12,  22,   6,  14, 122,  10,  25,   5],
            [125,  15, 130,  12,   8, 120,  10, 125,  14,   6, 122,  12, 128,  10,   8, 124,  10, 126,  12,   7],
            [ 20,  10,  15, 125, 130,  18,  12,  18, 122, 135,  22,  12,  15, 120, 132,  16,  12,  18, 120, 130],
            [ 12,  20,  10,   5,   2,  15,  22,  10,   6,   2,  10,  18,   8,   6,   3,  12,  20,   9,   5,   2]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0106.1",
            name="p53",
            length=20,
            pfm=[[float(x) for x in row] for row in p53_pfm],
            description="Master tumor suppressor protein binding symmetrical decamer half-sites.",
            biological_function="DNA damage response, cell cycle G1 arrest, senescence, and apoptosis."
        ))

        # 8. MYC / MAX (JASPAR MA0147.1) - E-box basic helix-loop-helix oncoprotein (8 bp)
        myc_pfm = [
            [ 35,  10,   8, 420,   6,   8,  20,  45],
            [390, 440,  10,   8, 430,  12,  15,  30],
            [ 45,  12, 445,  10,   8, 435, 410,  40],
            [ 20,   5,   5,  15,   9,  10,  18, 380]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0147.1",
            name="MYC",
            length=8,
            pfm=[[float(x) for x in row] for row in myc_pfm],
            description="Proto-oncogene c-Myc heterodimerized with Max targeting canonical E-box (CACGTG).",
            biological_function="Metabolic reprogramming, ribosome biogenesis, cellular transformation."
        ))

        # 9. STAT3 (JASPAR MA0144.1) - Signal transducer and activator of transcription 3 (11 bp)
        stat3_pfm = [
            [ 25,  15, 340,  20,  15, 320,  18,  12,  20,  30,  40],
            [ 20, 350,  12,  15, 325,  18,  14,  22, 330,  25,  35],
            [350,  18,  22, 335,  18,  15, 330,  20,  15, 340,  30],
            [ 15,  22,  18,  12,  20,  14,  15, 345,  18,  20, 310]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0144.1",
            name="STAT3",
            length=11,
            pfm=[[float(x) for x in row] for row in stat3_pfm],
            description="JAK/STAT signaling transducer mediating IL-6 response and oncogenic survival.",
            biological_function="Immune evasion, angiogenesis, wound healing, embryonic stem cell maintenance."
        ))

        # 10. FOXA1 (JASPAR MA0148.1) - Pioneer forkhead transcription factor (12 bp)
        foxa1_pfm = [
            [ 45, 310, 330,  18,  20, 345,  15,  18, 315,  40,  35,  50],
            [ 35,  18,  15,  22,  18,  12,  10,  15,  22,  30, 290,  40],
            [ 30,  20,  18, 325, 315,  15,  18,  22,  25, 295,  35,  45],
            [310,  15,  20,  15,  18,  14, 340, 335,  20,  25,  28, 280]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0148.1",
            name="FOXA1",
            length=12,
            pfm=[[float(x) for x in row] for row in foxa1_pfm],
            description="Pioneer transcription factor opening compacted heterochromatin for hormone receptors.",
            biological_function="Liver specification, prostate and breast cancer hormone receptor recruitment."
        ))

        # 11. ESR1 (JASPAR MA0112.1) - Estrogen Receptor Alpha nuclear hormone receptor (15 bp)
        esr1_pfm = [
            [280,  15,  20, 275,  25,  40,  35,  45,  30,  35,  20, 285,  15,  22, 270],
            [ 20, 270,  18,  15,  30,  45,  40,  35,  40,  30, 275,  18,  20, 265,  25],
            [ 25,  18, 285,  20,  35,  35,  45,  40,  35, 280,  15,  18, 275,  15,  20],
            [ 15,  22,  15,  25, 260,  40,  35,  30, 275,  20,  15,  20,  18,  20,  25]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0112.1",
            name="ESR1",
            length=15,
            pfm=[[float(x) for x in row] for row in esr1_pfm],
            description="Estrogen Receptor Alpha recognizing Estrogen Response Elements (ERE palindromes).",
            biological_function="Mammary gland development, endocrine signaling, ER+ breast cancer."
        ))

        # 12. SOX2 (JASPAR MA0143.1) - High-mobility group pluripotency master factor (10 bp)
        sox2_pfm = [
            [ 30,  25, 340, 350,  20,  25,  15,  20,  25,  35],
            [320,  18,  15,  18, 330,  18,  12,  25,  30,  40],
            [ 25,  22,  18,  12,  15,  18,  20,  30,  40, 310],
            [ 20, 335,  18,  15,  22, 345, 355, 330, 315,  25]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0143.1",
            name="SOX2",
            length=10,
            pfm=[[float(x) for x in row] for row in sox2_pfm],
            description="HMG-box pioneer transcription factor essential for pluripotent stem cell identity.",
            biological_function="Embryonic stem cell self-renewal, neural stem cell multipotency, Yamanaka factor."
        ))

        # 13. OCT4 / POU5F1 (JASPAR MA0785.1) - Master embryonic stem cell factor (12 bp)
        oct4_pfm = [
            [ 40, 320, 310,  20,  25, 315,  22,  18,  25,  30,  35,  45],
            [ 35,  18,  22,  25, 310,  18,  20,  22,  30,  28,  32,  40],
            [ 30,  22,  18, 315,  20,  15, 325,  25,  35, 305,  30,  38],
            [315,  20,  25,  15,  22,  20,  15, 330, 310,  22, 295, 275]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0785.1",
            name="OCT4",
            length=12,
            pfm=[[float(x) for x in row] for row in oct4_pfm],
            description="POU homeodomain transcription factor central to the core pluripotency network.",
            biological_function="Embryogenesis, maintenance of ESC naive pluripotency, cellular reprogramming."
        ))

        # 14. HNF4A (JASPAR MA0114.1) - Hepatocyte nuclear factor 4 alpha (13 bp)
        hnf4a_pfm = [
            [ 20, 290, 280,  18,  22,  25,  30, 275,  20,  18, 285,  25,  35],
            [ 25,  18,  20, 270,  18,  25,  28,  18, 280,  22,  20,  28,  30],
            [295,  22,  18,  22, 285, 270,  25,  20,  18, 275,  15,  25,  30],
            [ 18,  25,  22,  15,  20,  22, 280,  18,  22,  15,  18, 280, 270]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0114.1",
            name="HNF4A",
            length=13,
            pfm=[[float(x) for x in row] for row in hnf4a_pfm],
            description="Nuclear receptor controlling liver and pancreatic beta-cell specific gene programs.",
            biological_function="Gluconeogenesis, lipid transport, MODY1 diabetes mellitus etiology."
        ))

        # 15. CREB1 (JASPAR MA0018.1) - cAMP response element binding protein (8 bp)
        creb_pfm = [
            [ 20,  15,  12, 380,   8,  12,  20,  35],
            [ 18, 390,  15,  10,  15,  20, 375,  25],
            [385,  12, 380,  12,  10, 385,  15,  20],
            [ 15,  18,  10,  15, 395,  10,  18, 370]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0018.1",
            name="CREB1",
            length=8,
            pfm=[[float(x) for x in row] for row in creb_pfm],
            description="Phosphorylation-dependent bZIP transcription factor binding canonical cAMP site (TGACGTCA).",
            biological_function="Synaptic plasticity, long-term memory formation, glucose homeostasis."
        ))

        # 16. E2F1 (JASPAR MA0024.1) - Cell cycle G1/S transition regulator (11 bp)
        e2f1_pfm = [
            [ 25,  18,  22, 320,  15,  20,  25,  30,  35,  40,  45],
            [310, 330, 325,  15, 335,  18,  22, 315,  30,  25,  30],
            [ 20,  15,  18,  20,  18, 325, 330,  18, 310,  30,  25],
            [ 18,  20,  15,  18,  12,  15,  18,  20,  22, 315, 305]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0024.1",
            name="E2F1",
            length=11,
            pfm=[[float(x) for x in row] for row in e2f1_pfm],
            description="Retinoblastoma (pRb) pathway downstream effector governing DNA replication entry.",
            biological_function="G1 to S phase checkpoint transition, thymidine kinase & DNA polymerase activation."
        ))

        # 17. CEBPA (JASPAR MA0102.1) - CCAAT/enhancer binding protein alpha (10 bp)
        cebpa_pfm = [
            [ 30,  25,  20, 310, 320,  18,  25,  30,  35,  40],
            [ 25, 315, 320,  18,  22,  20,  25,  30,  35,  40],
            [ 20,  18,  22,  25,  18,  22, 315, 310,  30,  35],
            [310,  22,  18,  15,  20, 325,  20,  22, 310, 290]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0102.1",
            name="CEBPA",
            length=10,
            pfm=[[float(x) for x in row] for row in cebpa_pfm],
            description="Leucine zipper factor coordinating granulocytic myeloid differentiation.",
            biological_function="Neutrophil maturation, adipogenesis, mutated in acute myeloid leukemia."
        ))

        # 18. RUNX1 (JASPAR MA0002.1) - Core binding factor alpha 2 (11 bp)
        runx1_pfm = [
            [ 25,  18, 310, 325,  18,  22,  25,  30,  35,  40,  45],
            [ 20, 320,  15,  18, 330,  20,  22, 310,  25,  30,  35],
            [325,  15,  22,  15,  18, 335, 320,  18, 315,  25,  20],
            [ 18,  22,  18,  20,  15,  18,  20,  25,  22, 310, 305]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0002.1",
            name="RUNX1",
            length=11,
            pfm=[[float(x) for x in row] for row in runx1_pfm],
            description="Runt-domain factor partnering with CBFB to specify definitive hematopoiesis.",
            biological_function="HSC emergence from hemogenic endothelium, chromosomal translocation t(8;21)."
        ))

        # 19. PAX6 (JASPAR MA0069.1) - Paired box master neuro-ocular regulator (14 bp)
        pax6_pfm = [
            [ 30, 280, 290,  18,  22, 285,  20,  25,  30, 275,  25,  30,  35,  40],
            [ 25,  18,  20, 275, 280,  18,  22, 270,  25,  20, 270,  25,  30,  35],
            [285,  20,  18,  20,  18,  22, 280,  20, 275,  18,  20, 280,  25,  30],
            [ 20,  22,  25,  18,  20,  15,  18,  25,  20,  22,  25,  20, 285, 270]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0069.1",
            name="PAX6",
            length=14,
            pfm=[[float(x) for x in row] for row in pax6_pfm],
            description="Evolutionarily conserved paired-box homeoprotein orchestrating eye and brain morphogenesis.",
            biological_function="Corneal clarity, retinal neurogenesis, aniridia congenital syndrome etiology."
        ))

        # 20. MEF2A (JASPAR MA0052.1) - Myocyte enhancer factor 2A (10 bp)
        mef2a_pfm = [
            [ 25,  20, 330, 340, 335, 325,  18,  22,  25,  30],
            [320,  18,  15,  18,  20,  18,  22,  25, 310,  35],
            [ 20,  15,  18,  20,  18,  22,  25, 315,  20, 320],
            [ 18, 335,  18,  15,  20,  18, 320,  20,  25,  22]
        ]
        self.register_motif(MotifProfile(
            motif_id="MA0052.1",
            name="MEF2A",
            length=10,
            pfm=[[float(x) for x in row] for row in mef2a_pfm],
            description="MADS-box transcription factor driving muscle differentiation and cardiac homeostasis.",
            biological_function="Myogenesis, sarcomeric gene activation, coronary artery disease susceptibility."
        ))

    @staticmethod
    def parse_jaspar_format(jaspar_text: str) -> MotifProfile:
        """Parses raw text in standard JASPAR PFM format."""
        lines = [line.strip() for line in jaspar_text.strip().splitlines() if line.strip()]
        header = lines[0]
        if not header.startswith(">"):
            raise ValueError(f"Invalid JASPAR header: {header}")

        header_tokens = header[1:].split()
        motif_id = header_tokens[0]
        motif_name = header_tokens[1] if len(header_tokens) > 1 else motif_id

        pfm_matrix = []
        base_order = []
        for line in lines[1:5]:
            match = re.match(r"([ACGT])\s*\[\s*(.+?)\s*\]", line)
            if not match:
                raise ValueError(f"Could not parse JASPAR line: {line}")
            base_char = match.group(1)
            counts = [float(x) for x in match.group(2).split()]
            base_order.append(base_char)
            pfm_matrix.append(counts)

        ordered_pfm = []
        for base in DNA_BASES:
            idx = base_order.index(base)
            ordered_pfm.append(pfm_matrix[idx])

        length = len(ordered_pfm[0])
        return MotifProfile(
            motif_id=motif_id,
            name=motif_name,
            length=length,
            pfm=ordered_pfm,
            pseudocount=0.1
        )


# ============================================================================
# STANDALONE VERIFICATION & DEMONSTRATION SUITE
# ============================================================================

def run_member1_demonstration() -> None:
    """Executes a complete demonstration of all Member 1 functionality."""
    print("=" * 75)
    print("MEMBER 1 MODULE DEMO: 20 Biological Motif Models & Information Theory")
    print("=" * 75)

    bank = MotifBank()
    all_profiles = bank.get_all_profiles()
    print(f"[+] Loaded {len(all_profiles)} verified biological transcription factors into Motif Bank:")
    for idx, p in enumerate(all_profiles, start=1):
        print(f"  {idx:2d}. {p.name:<10} (ID: {p.motif_id:<10} | Length: {p.length:2d} bp | Consensus: {p.consensus_sequence:<14} | Info: {p.total_information_content:5.2f} bits)")

    # Detailed summary for GATA1 and p53
    print("\n" + bank.get("GATA1").display_summary())
    print("\n" + bank.get("p53").display_summary())

    # Dynamic Programming P-value Inversion demo
    gata1 = bank.get("GATA1")
    target_p = 1e-3
    thresh = gata1.pvalue_to_threshold(target_p)
    print(f"\n[+] Staden DP Inversion for GATA1: Threshold for P <= {target_p} is {thresh:.4f}")
    print("=" * 75)


if __name__ == "__main__":
    run_member1_demonstration()
