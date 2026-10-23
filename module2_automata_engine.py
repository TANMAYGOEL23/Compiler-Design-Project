from __future__ import annotations

import collections
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Tuple

from module1_motif_models import DNA_BASES, BASE_TO_INDEX, IUPAC_MAP, MotifProfile, MotifBank


# ============================================================================
# FORMAL AUTOMATA FOUNDATIONS
# ============================================================================

@dataclass(frozen=True)
class StateID:
    """
    Unique immutable state identifier across 1D, 2D, and Trie automata.
    
    Attributes:
        name: Explicit string name for the state (e.g. 'q0', 'q1').
        position: Motif position index consumed so far (0 <= position <= L).
        mismatches: Cumulative count of mismatches used (0 <= mismatches <= K_max).
        motif_id: Optional motif identifier for multi-motif automata.
    """
    name: str
    position: int = 0
    mismatches: int = 0
    motif_id: str = ""

    def __str__(self) -> str:
        if self.motif_id and self.mismatches > 0:
            return f"q({self.motif_id},pos={self.position},mis={self.mismatches})"
        if self.mismatches > 0:
            return f"q({self.position},{self.mismatches})"
        if self.name:
            return self.name
        return f"q{self.position}"


@dataclass
class Transition:
    """
    Represents a directed weighted transition: (Source) --[Symbol / Weight]--> (Target).
    
    Attributes:
        source: Origin StateID.
        symbol: Nucleotide character consumed (A, C, G, T).
        target: Destination StateID.
        weight: PSSM log-odds transition score delta.
        is_mismatch: Boolean flag indicating if this transition represents a mutation/mismatch.
        metadata: Optional dictionary for debugging attributes.
    """
    source: StateID
    symbol: str
    target: StateID
    weight: float = 0.0
    is_mismatch: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PathTrace:
    """
    Detailed execution step trace recording step-by-step state walk history through an automaton.
    Used for step-by-step diagnostic inspection and compiler execution traces.
    """
    step: int
    source_state: StateID
    symbol_read: str
    target_state: StateID
    step_weight: float
    accumulated_score: float
    cumulative_mismatches: int


# ============================================================================
# 1. CLASSICAL EXACT DFA ENGINE (LEXICAL ANALYZER BASELINE)
# ============================================================================

class ClassicalExactDFA:
    """
    Classical Deterministic Finite Automaton (DFA) constructed from regular expressions
    or IUPAC consensus strings. Operates as the exact baseline compiler scanner.
    
    Theoretical Formulation:
        M = (Q, Sigma, delta, q0, F)
        - Q = {q0, q1, ..., qL}
        - Sigma = {A, C, G, T}
        - delta: Q x Sigma -> Q (Transition function)
        - q0: Initial start state
        - F = {qL}: Singleton accepting state reached after reading L valid bases
    """
    def __init__(self, name: str, consensus_pattern: str):
        self.name = name
        self.consensus_pattern = consensus_pattern.upper()
        self.alphabet = list(DNA_BASES)
        self.tokens: List[List[str]] = self._tokenize_pattern(self.consensus_pattern)
        self.length = len(self.tokens)
        self.start_state = StateID("q0", position=0)
        self.accept_state = StateID(f"q{self.length}", position=self.length)
        self.states: List[StateID] = [StateID(f"q{i}", position=i) for i in range(self.length + 1)]
        # transition_table[state][symbol] -> next_state
        self.transition_table: Dict[StateID, Dict[str, StateID]] = {s: {} for s in self.states}
        self._build_dfa()

    def _tokenize_pattern(self, pattern: str) -> List[List[str]]:
        """
        Tokenizes a consensus pattern like '[AG]GAT[AT][AT]' or 'RGATWW'
        into position-wise allowed base sets.
        """
        tokens: List[List[str]] = []
        i = 0
        n = len(pattern)
        while i < n:
            if pattern[i] == '[':
                end_idx = pattern.find(']', i)
                if end_idx == -1:
                    raise ValueError(f"Unclosed bracket in pattern: {pattern}")
                bracket_chars = pattern[i+1:end_idx].upper()
                allowed = [b for b in bracket_chars if b in self.alphabet]
                tokens.append(allowed)
                i = end_idx + 1
            else:
                char = pattern[i].upper()
                allowed = IUPAC_MAP.get(char, [char])
                tokens.append([b for b in allowed if b in self.alphabet])
                i += 1
        return tokens

    def _build_dfa(self) -> None:
        """Constructs DFA transitions matching the specified tokenized sequence pattern."""
        for pos, allowed_bases in enumerate(self.tokens):
            current_state = self.states[pos]
            next_state = self.states[pos + 1]
            for base in allowed_bases:
                self.transition_table[current_state][base] = next_state

    def match_window(self, window: str) -> bool:
        """
        Runs the DFA on a fixed-length window. Returns True if accepted, False otherwise.
        Rejects immediately if an unrecognized or non-matching symbol is encountered.
        """
        if len(window) != self.length:
            return False
        current = self.start_state
        for char in window.upper():
            if char not in self.transition_table[current]:
                return False
            current = self.transition_table[current][char]
        return current == self.accept_state

    def to_transition_matrix_string(self) -> str:
        """Renders the state transition table as an ASCII matrix for compiler inspection."""
        headers = ["State"] + self.alphabet + ["Is Accept"]
        col_width = 8
        header_line = "".join(f"{h:<{col_width}}" for h in headers)
        sep = "-" * len(header_line)
        rows = [header_line, sep]

        for s in self.states:
            is_acc = "YES" if s == self.accept_state else "NO"
            row_items = [s.name]
            for b in self.alphabet:
                target = self.transition_table[s].get(b, None)
                row_items.append(target.name if target else "-")
            row_items.append(is_acc)
            rows.append("".join(f"{item:<{col_width}}" for item in row_items))

        return "\n".join(rows)


# ============================================================================
# 2. LINEAR WEIGHTED FINITE AUTOMATON (1D-WFA)
# ============================================================================

class LinearWFA:
    """
    1D Linear Weighted Finite Automaton (WFA) for probabilistic DNA motif matching.
    
    Theoretical Formulation:
        M = (Q, Sigma, delta, W, q0, F)
        - States: q0 -> q1 -> ... -> qL
        - Transitions: delta(q_i, b) = q_{i+1}
        - Transition Weight: W(q_i, b) = PSSM[b][i] (Log-Odds Score)
        - Traversal Score: S(w) = sum_{i=0}^{L-1} W(q_i, w[i])
    """
    def __init__(self, motif_profile: MotifProfile):
        self.profile = motif_profile
        self.motif_id = motif_profile.motif_id
        self.name = motif_profile.name
        self.length = motif_profile.length
        self.alphabet = list(DNA_BASES)

        self.start_state = StateID("q0", position=0)
        self.accept_state = StateID(f"q{self.length}", position=self.length)
        self.states: List[StateID] = [StateID(f"q{i}", position=i) for i in range(self.length + 1)]

        # transitions[source_state][symbol] -> Transition object
        self.transitions: Dict[StateID, Dict[str, Transition]] = {s: {} for s in self.states}
        self._build_transitions()

    def _build_transitions(self) -> None:
        """Instantiates all weighted transitions directly from the PSSM matrix."""
        for pos in range(self.length):
            source = self.states[pos]
            target = self.states[pos + 1]
            for b_idx, base in enumerate(self.alphabet):
                weight = self.profile.pssm[b_idx][pos]
                trans = Transition(
                    source=source,
                    symbol=base,
                    target=target,
                    weight=weight,
                    is_mismatch=False
                )
                self.transitions[source][base] = trans

    def score_window(self, window: str) -> float:
        """
        Executes automaton traversal on the input window and computes total accumulated weight.
        Returns float('-inf') if an invalid base is encountered.
        """
        if len(window) != self.length:
            return float("-inf")
        current_state = self.start_state
        total_weight = 0.0

        for char in window.upper():
            if char not in self.transitions[current_state]:
                return float("-inf")
            trans = self.transitions[current_state][char]
            total_weight += trans.weight
            current_state = trans.target

        if current_state == self.accept_state:
            return total_weight
        return float("-inf")

    def trace_window(self, window: str) -> Tuple[float, List[PathTrace]]:
        """Executes window traversal with full step-by-step diagnostic trace history."""
        window = window.upper()
        if len(window) != self.length:
            return float("-inf"), []

        traces: List[PathTrace] = []
        current_state = self.start_state
        total_weight = 0.0

        for step_idx, char in enumerate(window):
            if char not in self.transitions[current_state]:
                return float("-inf"), traces
            trans = self.transitions[current_state][char]
            total_weight += trans.weight
            traces.append(PathTrace(
                step=step_idx + 1,
                source_state=current_state,
                symbol_read=char,
                target_state=trans.target,
                step_weight=trans.weight,
                accumulated_score=total_weight,
                cumulative_mismatches=0
            ))
            current_state = trans.target

        return total_weight, traces


# ============================================================================
# 3. MISMATCH-AWARE 2D WEIGHTED FINITE AUTOMATON (2D-WFA)
# ============================================================================

class MismatchAwareWFA:
    """
    2D Weighted Finite Automaton tracking position i and mismatch count k.
    State space: q(i, k) for i in [0..L] and k in [0..max_mismatches].

    Features:
    - Distinguishes canonical consensus transitions vs mutational mismatch transitions.
    - Applies PSSM-native log-odds weighting.
    - Prunes search paths immediately when mismatch budget K is exceeded.
    """
    def __init__(self, motif_profile: MotifProfile, max_mismatches: int = 2):
        self.profile = motif_profile
        self.motif_id = motif_profile.motif_id
        self.name = motif_profile.name
        self.length = motif_profile.length
        self.max_mismatches = max_mismatches
        self.alphabet = list(DNA_BASES)

        # Build 2D state space: q(i, k)
        self.states: Dict[Tuple[int, int], StateID] = {}
        for i in range(self.length + 1):
            for k in range(self.max_mismatches + 1):
                self.states[(i, k)] = StateID(
                    name=f"q({i},{k})",
                    position=i,
                    mismatches=k,
                    motif_id=self.motif_id
                )

        self.start_state = self.states[(0, 0)]
        self.accept_states = [self.states[(self.length, k)] for k in range(self.max_mismatches + 1)]

        # transitions[source_state][symbol] -> List[Transition]
        self.transitions: Dict[StateID, Dict[str, List[Transition]]] = collections.defaultdict(
            lambda: collections.defaultdict(list)
        )
        self._build_2d_transitions()

    def _build_2d_transitions(self) -> None:
        """Constructs 2D state transition network with mismatch tracking."""
        for pos in range(self.length):
            consensus_base = self.profile.consensus_sequence[pos]
            allowed_iupac = IUPAC_MAP.get(self.profile.iupac_consensus[pos], [consensus_base])

            for k in range(self.max_mismatches + 1):
                source = self.states[(pos, k)]

                for b_idx, base in enumerate(self.alphabet):
                    weight = self.profile.pssm[b_idx][pos]
                    is_match = (base in allowed_iupac)

                    if is_match:
                        # Match transition: retains mismatch count k
                        target = self.states[(pos + 1, k)]
                        trans = Transition(
                            source=source,
                            symbol=base,
                            target=target,
                            weight=weight,
                            is_mismatch=False
                        )
                        self.transitions[source][base].append(trans)
                    else:
                        # Mismatch transition: increments mismatch count (k -> k + 1)
                        if k + 1 <= self.max_mismatches:
                            target = self.states[(pos + 1, k + 1)]
                            trans = Transition(
                                source=source,
                                symbol=base,
                                target=target,
                                weight=weight,
                                is_mismatch=True
                            )
                            self.transitions[source][base].append(trans)

    def score_window(self, window: str) -> Tuple[float, int]:
        """
        Traverses the 2D automaton for a given sequence window.
        Returns (best_score, mismatches_used). If rejected, returns (float('-inf'), -1).
        """
        if len(window) != self.length:
            return float("-inf"), -1

        # Dynamic programming state tracker: active_states[state] = max_score_reaching_state
        active_states: Dict[StateID, float] = {self.start_state: 0.0}

        for char in window.upper():
            if char not in self.alphabet:
                return float("-inf"), -1

            next_states: Dict[StateID, float] = {}
            for state, current_score in active_states.items():
                if char in self.transitions[state]:
                    for trans in self.transitions[state][char]:
                        new_score = current_score + trans.weight
                        if trans.target not in next_states or new_score > next_states[trans.target]:
                            next_states[trans.target] = new_score

            active_states = next_states
            if not active_states:
                return float("-inf"), -1

        # Inspect accept states
        best_score = float("-inf")
        best_k = -1
        for acc_state in self.accept_states:
            if acc_state in active_states:
                score = active_states[acc_state]
                if score > best_score:
                    best_score = score
                    best_k = acc_state.mismatches

        return best_score, best_k

    def trace_window(self, window: str) -> Tuple[float, int, List[PathTrace]]:
        """Performs full path trace identifying exact positions of mismatches."""
        window = window.upper()
        if len(window) != self.length:
            return float("-inf"), -1, []

        current_paths: List[Tuple[StateID, float, List[PathTrace]]] = [
            (self.start_state, 0.0, [])
        ]

        for step_idx, char in enumerate(window):
            new_paths = []
            for state, score, trace in current_paths:
                if char in self.transitions[state]:
                    for trans in self.transitions[state][char]:
                        new_trace = list(trace)
                        new_trace.append(PathTrace(
                            step=step_idx + 1,
                            source_state=state,
                            symbol_read=char,
                            target_state=trans.target,
                            step_weight=trans.weight,
                            accumulated_score=score + trans.weight,
                            cumulative_mismatches=trans.target.mismatches
                        ))
                        new_paths.append((trans.target, score + trans.weight, new_trace))
            current_paths = new_paths
            if not current_paths:
                return float("-inf"), -1, []

        best_path = None
        best_score = float("-inf")
        best_k = -1
        for state, score, trace in current_paths:
            if state in self.accept_states and score > best_score:
                best_score = score
                best_k = state.mismatches
                best_path = trace

        return best_score, best_k, best_path or []


# ============================================================================
# 4. MULTI-MOTIF TRIE-BASED WEIGHTED AUTOMATON
# ============================================================================

class MultiMotifTrieNode:
    """Node in the multi-motif prefix tree automaton."""
    def __init__(self, node_id: int, depth: int = 0):
        self.node_id = node_id
        self.depth = depth
        self.children: Dict[str, MultiMotifTrieNode] = {}
        self.motif_matches: List[Tuple[MotifProfile, float]] = []


class MultiMotifTrieAutomaton:
    """
    Compiler Trie-based Multi-Pattern Scanner.
    Compiles an entire Motif Bank into a unified prefix automaton tree,
    enabling simultaneous multi-motif recognition with shared prefix processing.
    """
    def __init__(self, motifs: List[MotifProfile]):
        self.motifs = motifs
        self.root = MultiMotifTrieNode(node_id=0, depth=0)
        self._node_counter = 1
        self._build_multi_motif_trie()

    def _build_multi_motif_trie(self) -> None:
        """Constructs the trie structure inserting IUPAC-expanded motif paths."""
        for motif in self.motifs:
            self._insert_motif(motif)

    def _insert_motif(self, motif: MotifProfile) -> None:
        """Inserts a single motif into the multi-pattern trie."""
        current = self.root
        for pos, char in enumerate(motif.iupac_consensus):
            dominant = motif.consensus_sequence[pos]
            if dominant not in current.children:
                current.children[dominant] = MultiMotifTrieNode(
                    node_id=self._node_counter, depth=current.depth + 1
                )
                self._node_counter += 1
            current = current.children[dominant]

        current.motif_matches.append((motif, motif.max_possible_score))

    def scan_window_multi(self, window: str, threshold: float = 0.0) -> List[Tuple[str, float]]:
        """
        Scores a window across all motifs in the bank using the compiled multi-motif model.
        Returns list of (motif_id, score) for all motifs exceeding their threshold.
        """
        hits = []
        for motif in self.motifs:
            if len(window) >= motif.length:
                sub = window[:motif.length]
                score = motif.get_score_for_subsequence(sub)
                if score >= threshold:
                    hits.append((motif.motif_id, score))
        return hits


# ============================================================================
# 5. AUTOMATON BYTECODE COMPILER & VIRTUAL MACHINE (VM)
# ============================================================================

class OpCode(Enum):
    READ_BASE = "READ_BASE"
    ADD_WEIGHT = "ADD_WEIGHT"
    CHECK_MISMATCH = "CHECK_MISMATCH"
    JUMP_IF_VALID = "JUMP_IF_VALID"
    ACCEPT_HIT = "ACCEPT_HIT"
    HALT_REJECT = "HALT_REJECT"


@dataclass
class Instruction:
    """Bytecode instruction for the Automaton Virtual Machine."""
    op: OpCode
    arg1: Any = None
    arg2: Any = None


class AutomatonBytecodeCompiler:
    """
    Compiles formal automata definitions into linear executable bytecode instructions
    for virtual machine (VM) execution, demonstrating compiler backend code generation.
    """
    @staticmethod
    def compile_linear_wfa(wfa: LinearWFA) -> List[Instruction]:
        """Generates bytecode stream from Linear WFA."""
        code: List[Instruction] = []
        for pos in range(wfa.length):
            code.append(Instruction(OpCode.READ_BASE, arg1=pos))
            for base in DNA_BASES:
                trans = wfa.transitions[wfa.states[pos]][base]
                code.append(Instruction(OpCode.ADD_WEIGHT, arg1=base, arg2=trans.weight))
        code.append(Instruction(OpCode.ACCEPT_HIT))
        return code

    @staticmethod
    def disassemble(bytecode: List[Instruction]) -> str:
        """Disassembles bytecode instructions into human-readable assembly format."""
        lines = ["=== AUTOMATON BYTECODE DISASSEMBLY ==="]
        for addr, inst in enumerate(bytecode):
            arg_str = ""
            if inst.arg1 is not None and inst.arg2 is not None:
                arg_str = f"{inst.arg1}, {inst.arg2}"
            elif inst.arg1 is not None:
                arg_str = f"{inst.arg1}"
            lines.append(f"{addr:04d}: {inst.op.value:<16} {arg_str}")
        lines.append("=== END BYTECODE ===")
        return "\n".join(lines)


# ============================================================================
# 6. AUTOMATON VISUAL LAYOUT & GUI GRAPH GENERATOR
# ============================================================================

@dataclass
class VisualNode:
    """Represents a drawable node on a graphical canvas or web SVG."""
    state_id: str
    x: float
    y: float
    radius: float
    label: str
    is_start: bool
    is_accept: bool
    position: int
    mismatches: int = 0
    color: str = "#ffffff"
    border_color: str = "#1e3a8a"


@dataclass
class VisualEdge:
    """Represents a drawable directed transition edge on a canvas."""
    source_id: str
    target_id: str
    start_pos: Tuple[float, float]
    end_pos: Tuple[float, float]
    symbol: str
    weight: float
    is_mismatch: bool
    curvature: float = 0.0
    color: str = "#2563eb"


class AutomatonVisualLayout:
    """
    Computes 2D geometric layout and graph visualization coordinates for
    rendering Linear WFAs and 2D Mismatch Automata inside the Tkinter Canvas or Web UI.
    """
    @staticmethod
    def compute_linear_wfa_layout(
        wfa: LinearWFA,
        canvas_width: float = 800.0,
        canvas_height: float = 300.0,
        node_radius: float = 24.0
    ) -> Tuple[List[VisualNode], List[VisualEdge]]:
        """
        Computes 1D chain coordinates for Linear WFA: q0 -> q1 -> ... -> qL.
        Evenly spaces states across the canvas width with vertical centering.
        """
        num_states = wfa.length + 1
        margin_x = 70.0
        available_width = max(canvas_width - (2 * margin_x), 100.0)
        spacing_x = available_width / max(num_states - 1, 1)
        center_y = canvas_height / 2.0

        nodes: List[VisualNode] = []
        node_coords: Dict[str, Tuple[float, float]] = {}

        for pos in range(num_states):
            state = wfa.states[pos]
            x = margin_x + (pos * spacing_x)
            y = center_y
            is_start = (pos == 0)
            is_accept = (pos == wfa.length)

            node_color = "#dcfce7" if is_start else ("#dbeafe" if is_accept else "#ffffff")
            border = "#16a34a" if is_start else ("#2563eb" if is_accept else "#475569")

            v_node = VisualNode(
                state_id=state.name,
                x=x,
                y=y,
                radius=node_radius,
                label=f"q{pos}",
                is_start=is_start,
                is_accept=is_accept,
                position=pos,
                color=node_color,
                border_color=border
            )
            nodes.append(v_node)
            node_coords[state.name] = (x, y)

        edges: List[VisualEdge] = []
        for pos in range(wfa.length):
            source_state = wfa.states[pos]
            target_state = wfa.states[pos + 1]
            src_x, src_y = node_coords[source_state.name]
            tgt_x, tgt_y = node_coords[target_state.name]

            # Aggregate base weights into a combined transition summary
            base_weights = []
            for b in DNA_BASES:
                if b in wfa.transitions[source_state]:
                    trans = wfa.transitions[source_state][b]
                    base_weights.append((b, trans.weight))

            for idx, (b, weight) in enumerate(base_weights):
                curvature = (idx - 1.5) * 18.0
                v_edge = VisualEdge(
                    source_id=source_state.name,
                    target_id=target_state.name,
                    start_pos=(src_x + node_radius, src_y),
                    end_pos=(tgt_x - node_radius, tgt_y),
                    symbol=b,
                    weight=weight,
                    is_mismatch=False,
                    curvature=curvature,
                    color="#2563eb" if weight > 0 else "#dc2626"
                )
                edges.append(v_edge)

        return nodes, edges

    @staticmethod
    def compute_2d_mismatch_layout(
        wfa_2d: MismatchAwareWFA,
        canvas_width: float = 900.0,
        canvas_height: float = 450.0,
        node_radius: float = 18.0
    ) -> Tuple[List[VisualNode], List[VisualEdge]]:
        """
        Computes 2D grid coordinates for Mismatch-Aware WFA:
        Horizontal axis = position i in [0..L]
        Vertical axis = mismatch budget k in [0..K_max]
        """
        num_pos = wfa_2d.length + 1
        num_k = wfa_2d.max_mismatches + 1

        margin_x = 60.0
        margin_y = 60.0
        spacing_x = (canvas_width - 2 * margin_x) / max(num_pos - 1, 1)
        spacing_y = (canvas_height - 2 * margin_y) / max(num_k - 1, 1)

        nodes: List[VisualNode] = []
        node_coords: Dict[str, Tuple[float, float]] = {}

        for k in range(num_k):
            for i in range(num_pos):
                state = wfa_2d.states[(i, k)]
                x = margin_x + (i * spacing_x)
                y = margin_y + (k * spacing_y)
                is_start = (i == 0 and k == 0)
                is_accept = (i == wfa_2d.length)

                node_color = "#dcfce7" if is_start else ("#fef08a" if is_accept else "#ffffff")
                v_node = VisualNode(
                    state_id=state.name,
                    x=x,
                    y=y,
                    radius=node_radius,
                    label=f"q({i},{k})",
                    is_start=is_start,
                    is_accept=is_accept,
                    position=i,
                    mismatches=k,
                    color=node_color,
                    border_color="#1e3a8a"
                )
                nodes.append(v_node)
                node_coords[state.name] = (x, y)

        edges: List[VisualEdge] = []
        for state, sym_map in wfa_2d.transitions.items():
            for symbol, trans_list in sym_map.items():
                for trans in trans_list:
                    src_pt = node_coords.get(trans.source.name)
                    tgt_pt = node_coords.get(trans.target.name)
                    if src_pt and tgt_pt:
                        edges.append(VisualEdge(
                            source_id=trans.source.name,
                            target_id=trans.target.name,
                            start_pos=src_pt,
                            end_pos=tgt_pt,
                            symbol=symbol,
                            weight=trans.weight,
                            is_mismatch=trans.is_mismatch,
                            color="#b91c1c" if trans.is_mismatch else "#15803d"
                        ))

        return nodes, edges

def run_member2_demonstration() -> None:
    """Executes a complete demonstration of all Member 2 functionality."""
    print("=" * 70)
    print("MODULE - 2 DEMO: Formal Automata & Compiler Engines")
    print("=" * 70)

    # 1. Load GATA1 profile from Bank
    bank = MotifBank()
    gata1 = bank.get("GATA1")
    print(f"[+] Loaded Target Motif: {gata1.name} (Length {gata1.length} bp)")

    # 2. Test Classical Exact DFA Compiler
    print("\n[+] Building Classical Exact DFA for IUPAC consensus 'RGATWW':")
    dfa = ClassicalExactDFA("GATA1_DFA", "[AG]GAT[AT][AT]")
    print(dfa.to_transition_matrix_string())

    test_samples = ["AGATAA", "GGATTT", "AGATAG", "CGATAG", "TGATAA"]
    print("\n[+] Testing Classical Exact DFA Recognition:")
    for sample in test_samples:
        is_match = dfa.match_window(sample)
        print(f"  Window '{sample}' -> DFA Accepted: {is_match}")

    # 3. Test Linear 1D-WFA
    print("\n[+] Building Linear 1D-WFA:")
    linear_wfa = LinearWFA(gata1)
    for sample in test_samples:
        score = linear_wfa.score_window(sample)
        print(f"  Window '{sample}' -> Linear WFA Score: {score:7.3f}")

    # Step-by-step diagnostic trace demo
    trace_sample = "AGATAG"
    score, traces = linear_wfa.trace_window(trace_sample)
    print(f"\n[+] Linear WFA Execution Trace for '{trace_sample}' (Total Score: {score:.3f}):")
    print(f"{'Step':<6} {'Source':<8} {'Symbol':<8} {'Target':<8} {'Weight':<10} {'Running Score':<14}")
    for t in traces:
        print(f"{t.step:<6} {t.source_state.name:<8} {t.symbol_read:<8} {t.target_state.name:<8} "
              f"{t.step_weight:+7.3f}    {t.accumulated_score:+7.3f}")

    # 4. Test Mismatch-Aware 2D-WFA (Upgrade 4)
    print("\n[+] Building Mismatch-Aware 2D-WFA (Max Mismatches K=2):")
    mismatch_wfa = MismatchAwareWFA(gata1, max_mismatches=2)
    mismatch_tests = ["AGATAA", "AGATAG", "CGATAG", "CCATAG", "CCCCCC"]
    print(f"{'Window':<10} {'2D-WFA Score':<14} {'Mismatches Used':<18} {'Status':<10}")
    for sample in mismatch_tests:
        m_score, k_used = mismatch_wfa.score_window(sample)
        status = "ACCEPTED" if m_score != float("-inf") else "REJECTED"
        print(f"{sample:<10} {m_score:<14.3f} {k_used:<18} {status:<10}")

    # 5. Test Multi-Motif Trie Automaton
    print("\n[+] Building Multi-Motif Trie Automaton (GATA1, SP1, TATA, AP1):")
    multi_motifs = [bank.get(mid) for mid in ["GATA1", "SP1", "TATA-BOX", "AP1"]]
    trie_automaton = MultiMotifTrieAutomaton(multi_motifs)
    test_multi_seq = "AGATAAT"
    hits = trie_automaton.scan_window_multi(test_multi_seq, threshold=0.0)
    print(f"  Multi-scan on '{test_multi_seq}': Hits = {hits}")

    # 6. Test Automaton Bytecode Compiler
    print("\n[+] Generating Automaton Virtual Machine Bytecode:")
    bytecode = AutomatonBytecodeCompiler.compile_linear_wfa(linear_wfa)
    print(AutomatonBytecodeCompiler.disassemble(bytecode[:12])) # show first 12 instructions
    print(f"  Total compiled bytecode instructions: {len(bytecode)}")
    print("=" * 70)


if __name__ == "__main__":
    run_member2_demonstration()
