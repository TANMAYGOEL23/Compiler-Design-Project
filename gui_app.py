"""
Automata-Based DNA Motif Detection: Interactive Graphical User Interface (GUI)
==============================================================================
Modern Tkinter/TTK Desktop Application for Compiler Design & DNA Motif Detection

Features:
1. 20-Motif Bank Explorer: Inspect PFMs, PSSMs, Information Content, and Staden DP P-values.
2. 20-Genome Reference Registry: Inspect real human & model organism loci.
3. Custom DNA Sequence Editor & Sequence Appender ("Add / Append DNA to Selected Genome").
4. Interactive Automaton Graph Visualizer & Diagnostic Single-Window Step Tracer.
5. High-Throughput Dual-Strand (+ and -) Scanner with Treeview Results and Surrounding Context Viewer.
6. Multi-format Exporters: BED6, GFF3, CSV, and JSON.
7. Ground-Truth Confusion Matrix & Multi-Engine Scalability Profiler.
"""

from __future__ import annotations

import math
import os
import sys
import time
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

from module1_motif_models import DNA_BASES, BASE_TO_INDEX, MotifProfile, MotifBank
from module2_automata_engine import (
    LinearWFA,
    MismatchAwareWFA,
    ClassicalExactDFA,
    AutomatonVisualLayout,
    AutomatonBytecodeCompiler
)
from module3_scanner_pipeline import (
    GenomicSequence,
    GenomeRegistry,
    DualStrandScannerPipeline,
    clean_dna_sequence
)
from module4_eval_benchmark import (
    GroundTruthGenerator,
    StatisticalEvaluator,
    ThresholdAnalyzer,
    ScalabilityProfiler
)


class MotifDetectionGUI:
    """Main Application Window for Automata-Based DNA Motif Detection."""
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("Automata-Based Probabilistic DNA Motif Detection — Compiler Design Suite")
        self.root.geometry("1200x820")
        self.root.minsize(1050, 720)

        # Initialize Core Registries
        self.motif_bank = MotifBank()
        self.genome_registry = GenomeRegistry()
        self.current_motif: MotifProfile = self.motif_bank.get("GATA1")
        self.current_genome: GenomicSequence = self.genome_registry.get("GENOME_HUMAN_HBB_CHR11")
        self.current_hits = []

        # Configure TTK Modern Style
        self._configure_styles()

        # Build UI Layout
        self._build_header()
        self._build_notebook()
        self._build_statusbar()

        # Initial Load
        self._on_motif_selected()
        self._on_genome_selected()

    def _configure_styles(self) -> None:
        """Applies modern clean styling to TTK widgets."""
        self.style = ttk.Style()
        try:
            self.style.theme_use("clam")
        except Exception:
            pass

        # Color Palette
        bg_main = "#f8fafc"
        navy_dark = "#0f172a"
        blue_primary = "#2563eb"
        blue_light = "#dbeafe"

        self.root.configure(bg=bg_main)

        self.style.configure(".", background=bg_main, foreground="#1e293b", font=("Helvetica", 10))
        self.style.configure("TNotebook", background=bg_main, tabmargins=[4, 4, 4, 0])
        self.style.configure("TNotebook.Tab", background="#e2e8f0", foreground="#334155",
                             padding=[14, 7], font=("Helvetica", 10, "bold"))
        self.style.map("TNotebook.Tab",
                       background=[("selected", blue_primary)],
                       foreground=[("selected", "#ffffff")])

        self.style.configure("Header.TFrame", background=navy_dark)
        self.style.configure("Card.TFrame", background="#ffffff", relief="solid", borderwidth=1)
        self.style.configure("Primary.TButton", background=blue_primary, foreground="#ffffff",
                             font=("Helvetica", 10, "bold"), padding=[10, 5])
        self.style.map("Primary.TButton",
                       background=[("active", "#1d4ed8"), ("disabled", "#94a3b8")])

        self.style.configure("Accent.TButton", background="#10b981", foreground="#ffffff",
                             font=("Helvetica", 10, "bold"), padding=[10, 5])
        self.style.map("Accent.TButton", background=[("active", "#059669")])

        self.style.configure("Treeview.Heading", font=("Helvetica", 10, "bold"), background="#e2e8f0")
        self.style.configure("Treeview", rowheight=24, font=("Consolas", 9))

    def _build_header(self) -> None:
        """Renders the top branding header banner."""
        header_frame = tk.Frame(self.root, bg="#0f172a", height=65)
        header_frame.pack(side="top", fill="x")

        title_lbl = tk.Label(header_frame, text="🧬 Automata-Based Probabilistic DNA Motif Detection",
                             bg="#0f172a", fg="#ffffff", font=("Helvetica", 15, "bold"))
        title_lbl.pack(side="left", padx=16, pady=8)

        subtitle_lbl = tk.Label(header_frame,
                                text="Compiler Design & Theory of Computation Case Study  |  20 Motifs & 20 Genomes",
                                bg="#0f172a", fg="#94a3b8", font=("Helvetica", 9))
        subtitle_lbl.pack(side="left", padx=8, pady=12)

        # Quick Status Badge
        self.badge_lbl = tk.Label(header_frame, text="Active: GATA1 (6 bp) on Human Chr11",
                                  bg="#1e293b", fg="#38bdf8", font=("Helvetica", 9, "bold"),
                                  padx=10, pady=4, relief="groove")
        self.badge_lbl.pack(side="right", padx=16, pady=10)

    def _build_statusbar(self) -> None:
        """Bottom status indicator bar."""
        self.status_bar = tk.Label(self.root, text="Ready. Select a transcription factor motif and genomic locus to begin.",
                                   bd=1, relief="sunken", anchor="w", bg="#f1f5f9", fg="#475569",
                                   font=("Helvetica", 9), padx=10, pady=4)
        self.status_bar.pack(side="bottom", fill="x")

    def _build_notebook(self) -> None:
        """Creates the main multi-tab container."""
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=12, pady=8)

        # Tab 1: Motif & Genome Selector + Sequence Appender
        self.tab_explorer = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_explorer, text="  1. Motif & Genome Explorer  ")
        self._build_tab_explorer()

        # Tab 2: Automata Engine Visualizer
        self.tab_automata = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_automata, text="  2. Automata Engine Visualizer  ")
        self._build_tab_automata()

        # Tab 3: Dual-Strand Scanner & Results
        self.tab_scanner = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_scanner, text="  3. Genomic Scanner & Results  ")
        self._build_tab_scanner()

        # Tab 4: Evaluation & Performance
        self.tab_eval = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_eval, text="  4. Evaluation & Benchmarks  ")
        self._build_tab_eval()

    # =========================================================================
    # TAB 1: MOTIF & GENOME EXPLORER WITH DNA APPENDING
    # =========================================================================

    def _build_tab_explorer(self) -> None:
        """Constructs Tab 1 for exploring 20 Motifs, 20 Genomes, and Appending Sequences."""
        paned = ttk.PanedWindow(self.tab_explorer, orient="horizontal")
        paned.pack(fill="both", expand=True, padx=8, pady=8)

        # --- LEFT PANEL: 20-MOTIF SELECTION ---
        left_frame = ttk.LabelFrame(paned, text=" Transcription Factor Motif (20 Options) ", padding=10)
        paned.add(left_frame, weight=1)

        ttk.Label(left_frame, text="Select Target Motif:", font=("Helvetica", 10, "bold")).pack(anchor="w")

        # Combobox of 20 Motifs
        all_motifs = self.motif_bank.get_all_profiles()
        motif_labels = [f"{p.name} ({p.motif_id}) - {p.length} bp" for p in all_motifs]
        self.motif_combo = ttk.Combobox(left_frame, values=motif_labels, state="readonly", width=32)
        self.motif_combo.current(7) # default to GATA1
        self.motif_combo.pack(fill="x", pady=5)
        self.motif_combo.bind("<<ComboboxSelected>>", self._on_motif_combo_change)

        # Motif Details Display
        self.motif_info_txt = tk.Text(left_frame, height=18, font=("Consolas", 9), wrap="word",
                                      bg="#ffffff", relief="solid", bd=1)
        self.motif_info_txt.pack(fill="both", expand=True, pady=6)

        # --- RIGHT PANEL: 20-GENOME SELECTION & CUSTOM DNA APPENDER ---
        right_frame = ttk.LabelFrame(paned, text=" Genomic Target & Custom Sequence Editor ", padding=10)
        paned.add(right_frame, weight=2)

        ttk.Label(right_frame, text="Select Genomic Target (20 Preloaded Options):",
                  font=("Helvetica", 10, "bold")).pack(anchor="w")

        all_genomes = self.genome_registry.get_all_records()
        genome_labels = [f"{g.locus_name} [{g.organism}] ({g.length:,} bp)" for g in all_genomes]
        self.genome_combo = ttk.Combobox(right_frame, values=genome_labels, state="readonly", width=45)
        self.genome_combo.current(0)
        self.genome_combo.pack(fill="x", pady=5)
        self.genome_combo.bind("<<ComboboxSelected>>", self._on_genome_combo_change)

        # Genome Info Banner
        self.genome_meta_lbl = ttk.Label(right_frame, text="", font=("Helvetica", 9, "italic"))
        self.genome_meta_lbl.pack(anchor="w", pady=2)

        # Live Sequence Display
        ttk.Label(right_frame, text="Sequence Viewer (5' -> 3'):", font=("Helvetica", 9, "bold")).pack(anchor="w")
        self.seq_display_txt = tk.Text(right_frame, height=8, font=("Consolas", 9), wrap="char",
                                       bg="#f8fafc", relief="solid", bd=1)
        self.seq_display_txt.pack(fill="x", pady=4)

        # --- DNA SEQUENCE ADDER / APPENDER PANEL ---
        append_box = ttk.LabelFrame(right_frame, text=" ➕ Add or Append DNA Sequence ", padding=8)
        append_box.pack(fill="both", expand=True, pady=8)

        ctrl_row = ttk.Frame(append_box)
        ctrl_row.pack(fill="x", pady=2)

        self.append_mode_var = tk.StringVar(value="append")
        ttk.Radiobutton(ctrl_row, text="Append sequence to current genome",
                        variable=self.append_mode_var, value="append").pack(side="left", padx=5)
        ttk.Radiobutton(ctrl_row, text="Create new custom genome",
                        variable=self.append_mode_var, value="new").pack(side="left", padx=10)

        ttk.Label(append_box, text="Enter / Paste DNA nucleotides (A, C, G, T):", font=("Helvetica", 8)).pack(anchor="w")
        self.append_input_txt = tk.Text(append_box, height=4, font=("Consolas", 9), wrap="char",
                                        bg="#ffffff", relief="solid", bd=1)
        self.append_input_txt.pack(fill="x", pady=4)

        btn_row = ttk.Frame(append_box)
        btn_row.pack(fill="x", pady=4)

        ttk.Button(btn_row, text="Paste Example Regulatory DNA",
                   command=self._load_sample_dna_snippet).pack(side="left", padx=4)
        ttk.Button(btn_row, text="Apply / Append Sequence", style="Primary.TButton",
                   command=self._apply_sequence_update).pack(side="right", padx=4)

    def _on_motif_combo_change(self, event=None) -> None:
        idx = self.motif_combo.current()
        all_motifs = self.motif_bank.get_all_profiles()
        if 0 <= idx < len(all_motifs):
            self.current_motif = all_motifs[idx]
            self._on_motif_selected()

    def _on_genome_combo_change(self, event=None) -> None:
        idx = self.genome_combo.current()
        all_genomes = self.genome_registry.get_all_records()
        if 0 <= idx < len(all_genomes):
            self.current_genome = all_genomes[idx]
            self._on_genome_selected()

    def _on_motif_selected(self) -> None:
        """Refreshes UI panels with the currently selected motif."""
        p = self.current_motif
        self.badge_lbl.config(text=f"Active: {p.name} ({p.length} bp) on {self.current_genome.locus_name}")

        self.motif_info_txt.delete("1.0", "end")
        self.motif_info_txt.insert("end", f"TRANSCRIPTION FACTOR: {p.name}\n")
        self.motif_info_txt.insert("end", f"JASPAR Profile:   {p.motif_id}\n")
        self.motif_info_txt.insert("end", f"Organism:         {p.organism}\n")
        self.motif_info_txt.insert("end", f"Length:           {p.length} bp\n")
        self.motif_info_txt.insert("end", f"Strict Consensus: {p.consensus_sequence}\n")
        self.motif_info_txt.insert("end", f"IUPAC Consensus:  {p.iupac_consensus}\n")
        self.motif_info_txt.insert("end", f"Regex Pattern:    {p.regex_pattern}\n")
        self.motif_info_txt.insert("end", f"Info Content:     {p.total_information_content:.3f} bits (Max {2*p.length})\n")
        self.motif_info_txt.insert("end", f"Score Range:      [{p.min_possible_score:.2f} .. {p.max_possible_score:.2f}]\n\n")
        self.motif_info_txt.insert("end", f"BIOLOGICAL ROLE:\n{p.biological_function}\n\n")
        self.motif_info_txt.insert("end", "POSITION-SPECIFIC SCORING MATRIX (PSSM):\n")
        for b_idx, base in enumerate(DNA_BASES):
            row_str = " ".join(f"{p.pssm[b_idx][j]:6.2f}" for j in range(p.length))
            self.motif_info_txt.insert("end", f"  {base} | {row_str}\n")

        # Update scanner threshold recommendation
        if hasattr(self, "thresh_entry"):
            opt_thresh = max(3.0, round(p.max_possible_score * 0.60, 2))
            self.thresh_entry.delete(0, "end")
            self.thresh_entry.insert(0, str(opt_thresh))

        self._draw_automaton_canvas()

    def _on_genome_selected(self) -> None:
        """Refreshes UI panels with the currently selected genome."""
        g = self.current_genome
        self.badge_lbl.config(text=f"Active: {self.current_motif.name} ({self.current_motif.length} bp) on {g.locus_name}")

        self.genome_meta_lbl.config(
            text=f"Organism: {g.organism}  |  Length: {g.length:,} bp  |  GC Content: {g.get_gc_content():.2%}  |  Coords: {g.start_coordinate:,} - {g.end_coordinate:,}"
        )

        self.seq_display_txt.delete("1.0", "end")
        preview_text = g.sequence[:500] + ("..." if len(g.sequence) > 500 else "")
        self.seq_display_txt.insert("end", preview_text)

    def _load_sample_dna_snippet(self) -> None:
        """Loads a realistic regulatory DNA segment into the text editor."""
        sample_fragment = "AGATAAGGGGCGGGGGCACGGCACGGCACGGTTTTATAAAATTTTTCACGTG"
        self.append_input_txt.delete("1.0", "end")
        self.append_input_txt.insert("end", sample_fragment)

    def _apply_sequence_update(self) -> None:
        """Applies sequence addition or appending as requested by user."""
        raw_input = self.append_input_txt.get("1.0", "end").strip()
        cleaned = clean_dna_sequence(raw_input)

        if not cleaned:
            messagebox.showwarning("Empty Input", "Please enter valid DNA bases (A, C, G, T).")
            return

        mode = self.append_mode_var.get()
        if mode == "append":
            prev_len = self.current_genome.length
            new_len = self.current_genome.append_sequence(cleaned)
            self._on_genome_selected()
            self.status_bar.config(text=f"Successfully appended {len(cleaned):,} bp to {self.current_genome.locus_name} (Length: {prev_len:,} -> {new_len:,} bp).")
            messagebox.showinfo("Sequence Appended",
                                f"Appended {len(cleaned):,} bp to {self.current_genome.locus_name}.\nNew Total Length: {new_len:,} bp.")
        else:
            custom_id = f"CUSTOM_GENOME_{int(time.time())}"
            new_record = self.genome_registry.add_custom_genome(
                sequence_id=custom_id,
                sequence_data=cleaned,
                description="User Added Custom Genomic Sequence"
            )
            self.current_genome = new_record
            # Refresh genome dropdown
            all_genomes = self.genome_registry.get_all_records()
            labels = [f"{g.locus_name} [{g.organism}] ({g.length:,} bp)" for g in all_genomes]
            self.genome_combo.config(values=labels)
            self.genome_combo.current(len(labels) - 1)
            self._on_genome_selected()
            self.status_bar.config(text=f"Created custom genome '{custom_id}' ({len(cleaned):,} bp).")
            messagebox.showinfo("Genome Created", f"Successfully registered new genome: {custom_id} ({len(cleaned):,} bp).")

        self.append_input_txt.delete("1.0", "end")

    # =========================================================================
    # TAB 2: AUTOMATA ENGINE VISUALIZER & DIAGNOSTIC STEP TRACER
    # =========================================================================

    def _build_tab_automata(self) -> None:
        """Constructs Tab 2 for visual state transition inspection."""
        top_ctrl = ttk.Frame(self.tab_automata, padding=8)
        top_ctrl.pack(fill="x")

        ttk.Label(top_ctrl, text="Test Single Sequence Window:", font=("Helvetica", 10, "bold")).pack(side="left", padx=6)
        self.test_window_entry = ttk.Entry(top_ctrl, font=("Consolas", 10), width=24)
        self.test_window_entry.insert(0, "AGATAG")
        self.test_window_entry.pack(side="left", padx=6)

        ttk.Button(top_ctrl, text="Trace Automaton Walk", style="Primary.TButton",
                   command=self._trace_test_window).pack(side="left", padx=8)

        self.trace_score_lbl = ttk.Label(top_ctrl, text="Score: -- | Status: --", font=("Helvetica", 10, "bold"))
        self.trace_score_lbl.pack(side="left", padx=15)

        # Canvas for Graph Drawing
        canvas_box = ttk.LabelFrame(self.tab_automata, text=" Automaton State Transition Network ", padding=6)
        canvas_box.pack(fill="x", padx=8, pady=4)

        self.canvas = tk.Canvas(canvas_box, height=220, bg="#ffffff", highlightthickness=1, highlightbackground="#cbd5e1")
        self.canvas.pack(fill="both", expand=True)

        # Step Trace Table & Disassembly Split
        split_frame = ttk.PanedWindow(self.tab_automata, orient="horizontal")
        split_frame.pack(fill="both", expand=True, padx=8, pady=6)

        # Left: Step Trace Table
        trace_frame = ttk.LabelFrame(split_frame, text=" Step-by-Step Diagnostic Trace ", padding=6)
        split_frame.add(trace_frame, weight=3)

        cols = ("Step", "Source", "Symbol", "Target", "Weight", "Running Score")
        self.trace_tree = ttk.Treeview(trace_frame, columns=cols, show="headings", height=8)
        for c in cols:
            self.trace_tree.heading(c, text=c)
            self.trace_tree.column(c, width=100, anchor="center")
        self.trace_tree.pack(fill="both", expand=True)

        # Right: Bytecode VM Disassembly
        vm_frame = ttk.LabelFrame(split_frame, text=" Automaton VM Bytecode Disassembly ", padding=6)
        split_frame.add(vm_frame, weight=2)

        self.vm_text = tk.Text(vm_frame, font=("Consolas", 9), bg="#f8fafc", relief="solid", bd=1)
        self.vm_text.pack(fill="both", expand=True)

    def _draw_automaton_canvas(self, active_path_states: Optional[List[str]] = None) -> None:
        """Renders the automaton states and transitions on the Tkinter canvas."""
        if not hasattr(self, "canvas"):
            return
        self.canvas.delete("all")

        wfa = LinearWFA(self.current_motif)
        w = max(self.canvas.winfo_width(), 850)
        h = max(self.canvas.winfo_height(), 200)

        nodes, edges = AutomatonVisualLayout.compute_linear_wfa_layout(wfa, canvas_width=w, canvas_height=h)

        # Draw transitions (edges)
        for e in edges:
            x1, y1 = e.start_pos
            x2, y2 = e.end_pos
            cy = y1 + e.curvature

            # Draw curved line using quadratic Bezier approximation
            self.canvas.create_line(x1, y1, (x1 + x2)/2, cy, x2, y2,
                                    smooth=True, arrow=tk.LAST, arrowshape=(8, 10, 3),
                                    fill="#94a3b8" if abs(e.curvature) > 10 else e.color, width=1.5)
            # Label on top
            if e.symbol in ["A", "T"]:
                self.canvas.create_text((x1 + x2)/2, cy - 8, text=f"{e.symbol}/{e.weight:+.1f}",
                                        font=("Helvetica", 7), fill="#334155")

        # Draw states (nodes)
        for n in nodes:
            r = n.radius
            fill_c = n.color
            if active_path_states and n.state_id in active_path_states:
                fill_c = "#fde047" # Highlight in yellow

            self.canvas.create_oval(n.x - r, n.y - r, n.x + r, n.y + r,
                                    fill=fill_c, outline=n.border_color, width=2)
            if n.is_accept:
                # Double circle for accept state
                self.canvas.create_oval(n.x - r + 3, n.y - r + 3, n.x + r - 3, n.y + r - 3,
                                        outline=n.border_color, width=1.5)

            self.canvas.create_text(n.x, n.y, text=n.label, font=("Helvetica", 9, "bold"), fill="#0f172a")

        # Update bytecode disassembly
        if hasattr(self, "vm_text"):
            bytecode = AutomatonBytecodeCompiler.compile_linear_wfa(wfa)
            disasm = AutomatonBytecodeCompiler.disassemble(bytecode[:18])
            self.vm_text.delete("1.0", "end")
            self.vm_text.insert("end", disasm)

    def _trace_test_window(self) -> None:
        """Performs step-by-step diagnostic traversal of an input window."""
        window = self.test_window_entry.get().strip().upper()
        if len(window) != self.current_motif.length:
            messagebox.showwarning("Length Mismatch",
                                   f"Input window length ({len(window)}) must match motif length ({self.current_motif.length} bp).")
            return

        wfa = LinearWFA(self.current_motif)
        score, traces = wfa.trace_window(window)

        # Clear table
        for item in self.trace_tree.get_children():
            self.trace_tree.delete(item)

        active_states = [wfa.start_state.name]
        for t in traces:
            self.trace_tree.insert("", "end", values=(
                t.step, t.source_state.name, t.symbol_read, t.target_state.name,
                f"{t.step_weight:+.3f}", f"{t.accumulated_score:+.3f}"
            ))
            active_states.append(t.target_state.name)

        if score != float("-inf"):
            self.trace_score_lbl.config(text=f"Score: {score:.3f} | Status: ACCEPTED", foreground="#16a34a")
        else:
            self.trace_score_lbl.config(text="Score: -inf | Status: REJECTED", foreground="#dc2626")

        self._draw_automaton_canvas(active_path_states=active_states)

    # =========================================================================
    # TAB 3: GENOMIC SCANNER & DUAL-STRAND RESULTS TABLE
    # =========================================================================

    def _build_tab_scanner(self) -> None:
        """Constructs Tab 3 for running dual-strand genomic scans and inspecting results."""
        ctrl_card = ttk.LabelFrame(self.tab_scanner, text=" Scanning Engine Parameters ", padding=10)
        ctrl_card.pack(fill="x", padx=8, pady=6)

        row1 = ttk.Frame(ctrl_card)
        row1.pack(fill="x", pady=2)

        ttk.Label(row1, text="Score Threshold:", font=("Helvetica", 10, "bold")).pack(side="left", padx=5)
        self.thresh_entry = ttk.Entry(row1, width=8, font=("Helvetica", 10))
        self.thresh_entry.insert(0, "4.0")
        self.thresh_entry.pack(side="left", padx=5)

        ttk.Button(row1, text="Use Optimal Staden DP Cutoff (P <= 0.001)",
                   command=self._set_optimal_threshold).pack(side="left", padx=10)

        ttk.Label(row1, text="Mismatch Budget K:", font=("Helvetica", 10, "bold")).pack(side="left", padx=15)
        self.k_combo = ttk.Combobox(row1, values=["0 (Exact)", "1 Mismatch", "2 Mismatches", "3 Mismatches"],
                                    state="readonly", width=14)
        self.k_combo.current(2)
        self.k_combo.pack(side="left", padx=5)

        row2 = ttk.Frame(ctrl_card)
        row2.pack(fill="x", pady=6)

        self.dual_strand_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(row2, text="Scan Both Strands (+ and -)", variable=self.dual_strand_var).pack(side="left", padx=5)

        self.pruning_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(row2, text="Lookahead Branch-and-Bound Pruning", variable=self.pruning_var).pack(side="left", padx=15)

        self.use_2d_wfa_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(row2, text="Use 2D Mismatch-Aware WFA Engine", variable=self.use_2d_wfa_var).pack(side="left", padx=15)

        ttk.Button(row2, text="  ⚡ RUN AUTOMATA SCAN  ", style="Primary.TButton",
                   command=self._run_genomic_scan).pack(side="right", padx=6)

        # Metrics Bar
        self.metrics_frame = tk.Frame(self.tab_scanner, bg="#f1f5f9", height=32, bd=1, relief="groove")
        self.metrics_frame.pack(fill="x", padx=8, pady=2)
        self.hits_summary_lbl = tk.Label(self.metrics_frame, text="Total Hits: 0  |  Forward (+): 0  |  Reverse (-): 0  |  Exact DFA: 0  |  Time: 0.00 ms",
                                         bg="#f1f5f9", fg="#0f172a", font=("Helvetica", 9, "bold"))
        self.hits_summary_lbl.pack(side="left", padx=10, pady=4)

        # Results Table
        table_frame = ttk.Frame(self.tab_scanner)
        table_frame.pack(fill="both", expand=True, padx=8, pady=4)

        cols = ("#", "Start", "End", "Strand", "Sequence", "Score", "P-Value", "Mismatches", "Match Type")
        self.results_tree = ttk.Treeview(table_frame, columns=cols, show="headings", selectmode="browse")

        widths = [45, 95, 95, 60, 110, 80, 95, 85, 120]
        for c, w in zip(cols, widths):
            self.results_tree.heading(c, text=c)
            self.results_tree.column(c, width=w, anchor="center")

        tree_scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.results_tree.yview)
        self.results_tree.configure(yscrollcommand=tree_scroll.set)
        self.results_tree.pack(side="left", fill="both", expand=True)
        tree_scroll.pack(side="right", fill="y")

        self.results_tree.bind("<<TreeviewSelect>>", self._on_hit_row_select)

        # Genomic Context Preview & Export Bar
        bottom_bar = ttk.Frame(self.tab_scanner, padding=4)
        bottom_bar.pack(fill="x", padx=8, pady=4)

        self.context_lbl = ttk.Label(bottom_bar, text="Click any hit to preview genomic context (+/- 20 bp): --",
                                     font=("Consolas", 9, "bold"), foreground="#1e3a8a")
        self.context_lbl.pack(side="left", padx=4)

        ttk.Button(bottom_bar, text="Export BED6", command=self._export_bed).pack(side="right", padx=3)
        ttk.Button(bottom_bar, text="Export GFF3", command=self._export_gff3).pack(side="right", padx=3)
        ttk.Button(bottom_bar, text="Export CSV", command=self._export_csv).pack(side="right", padx=3)
        ttk.Button(bottom_bar, text="Export JSON", command=self._export_json).pack(side="right", padx=3)

    def _set_optimal_threshold(self) -> None:
        """Sets the operating threshold to the Staden DP calculated P <= 0.001 value."""
        t = self.current_motif.pvalue_to_threshold(1e-3)
        self.thresh_entry.delete(0, "end")
        self.thresh_entry.insert(0, f"{t:.2f}")

    def _run_genomic_scan(self) -> None:
        """Executes the dual-strand genomic scanning pipeline."""
        try:
            threshold = float(self.thresh_entry.get().strip())
        except ValueError:
            messagebox.showerror("Invalid Threshold", "Please enter a numeric score threshold.")
            return

        k_val = self.k_combo.current()
        scan_both = self.dual_strand_var.get()
        use_prune = self.pruning_var.get()
        use_2d = self.use_2d_wfa_var.get()

        pipeline = DualStrandScannerPipeline(self.current_motif, max_mismatches=k_val)

        t0 = time.perf_counter()
        hits = pipeline.scan_genomic_record(
            self.current_genome,
            threshold=threshold,
            scan_both_strands=scan_both,
            use_pruning=use_prune,
            use_mismatch_wfa=use_2d
        )
        scan_ms = (time.perf_counter() - t0) * 1000
        self.current_hits = hits

        # Populate Results Treeview
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)

        fwd_count = sum(1 for h in hits if h.strand == "+")
        rev_count = sum(1 for h in hits if h.strand == "-")
        exact_count = sum(1 for h in hits if h.is_exact_match)

        for idx, h in enumerate(hits, start=1):
            m_type = "Exact Match" if h.is_exact_match else "Probabilistic WFA"
            self.results_tree.insert("", "end", values=(
                idx, f"{h.start:,}", f"{h.end:,}", h.strand, h.matched_sequence,
                f"{h.score:.3f}", f"{h.p_value:.2e}", h.mismatches, m_type
            ))

        self.hits_summary_lbl.config(
            text=f"Total Hits: {len(hits)}  |  Forward (+): {fwd_count}  |  Reverse (-): {rev_count}  |  Exact DFA: {exact_count}  |  Scan Time: {scan_ms:.2f} ms"
        )
        self.status_bar.config(text=f"Scan complete: Found {len(hits)} candidate sites for {self.current_motif.name} in {self.current_genome.locus_name}.")

    def _on_hit_row_select(self, event=None) -> None:
        """Extracts and displays the surrounding genomic nucleotide context."""
        selected = self.results_tree.selection()
        if not selected:
            return
        item_vals = self.results_tree.item(selected[0], "values")
        hit_idx = int(item_vals[0]) - 1
        if 0 <= hit_idx < len(self.current_hits):
            h = self.current_hits[hit_idx]
            # Calculate local sequence offset
            local_start = h.start - self.current_genome.start_coordinate
            local_end = h.end - self.current_genome.start_coordinate
            full_seq = self.current_genome.sequence

            flank = 15
            left_flank = full_seq[max(0, local_start - flank):local_start]
            hit_seq = full_seq[local_start:local_end]
            right_flank = full_seq[local_end:min(len(full_seq), local_end + flank)]

            context_str = f"Context: ...{left_flank} [{hit_seq}] {right_flank}... ({h.strand} strand)"
            self.context_lbl.config(text=context_str)

    def _export_bed(self) -> None:
        if not self.current_hits:
            messagebox.showwarning("No Hits", "Please run a scan before exporting.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".bed", filetypes=[("BED files", "*.bed")])
        if path:
            pipeline = DualStrandScannerPipeline(self.current_motif)
            pipeline.export_hits_to_bed(self.current_hits, path)
            messagebox.showinfo("Export Successful", f"Saved {len(self.current_hits)} hits to BED:\n{path}")

    def _export_gff3(self) -> None:
        if not self.current_hits:
            messagebox.showwarning("No Hits", "Please run a scan before exporting.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".gff3", filetypes=[("GFF3 files", "*.gff3")])
        if path:
            pipeline = DualStrandScannerPipeline(self.current_motif)
            pipeline.export_hits_to_gff3(self.current_hits, path)
            messagebox.showinfo("Export Successful", f"Saved {len(self.current_hits)} hits to GFF3:\n{path}")

    def _export_csv(self) -> None:
        if not self.current_hits:
            messagebox.showwarning("No Hits", "Please run a scan before exporting.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if path:
            pipeline = DualStrandScannerPipeline(self.current_motif)
            pipeline.export_hits_to_csv(self.current_hits, path)
            messagebox.showinfo("Export Successful", f"Saved {len(self.current_hits)} hits to CSV:\n{path}")

    def _export_json(self) -> None:
        if not self.current_hits:
            messagebox.showwarning("No Hits", "Please run a scan before exporting.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON files", "*.json")])
        if path:
            pipeline = DualStrandScannerPipeline(self.current_motif)
            pipeline.export_hits_to_json(self.current_hits, path)
            messagebox.showinfo("Export Successful", f"Saved {len(self.current_hits)} hits to JSON:\n{path}")

    # =========================================================================
    # TAB 4: EVALUATION, CONFUSION MATRIX & MULTI-ENGINE PROFILER
    # =========================================================================

    def _build_tab_eval(self) -> None:
        """Constructs Tab 4 for ground-truth benchmarks and multi-engine profiling."""
        ctrl_frame = ttk.LabelFrame(self.tab_eval, text=" Ground-Truth Synthetic Evaluation Engine ", padding=10)
        ctrl_frame.pack(fill="x", padx=8, pady=6)

        row = ttk.Frame(ctrl_frame)
        row.pack(fill="x")

        ttk.Label(row, text="Synthetic DNA Length:", font=("Helvetica", 9, "bold")).pack(side="left", padx=5)
        self.syn_len_entry = ttk.Entry(row, width=10)
        self.syn_len_entry.insert(0, "15000")
        self.syn_len_entry.pack(side="left", padx=5)

        ttk.Label(row, text="Implanted True Sites:", font=("Helvetica", 9, "bold")).pack(side="left", padx=10)
        self.syn_sites_entry = ttk.Entry(row, width=6)
        self.syn_sites_entry.insert(0, "40")
        self.syn_sites_entry.pack(side="left", padx=5)

        ttk.Label(row, text="Mutation Rate:", font=("Helvetica", 9, "bold")).pack(side="left", padx=10)
        self.syn_mut_entry = ttk.Entry(row, width=6)
        self.syn_mut_entry.insert(0, "0.10")
        self.syn_mut_entry.pack(side="left", padx=5)

        ttk.Button(row, text="Run Statistical Evaluation", style="Primary.TButton",
                   command=self._run_statistical_eval).pack(side="right", padx=6)

        # Results Text Box
        self.eval_text = tk.Text(self.tab_eval, font=("Consolas", 9), bg="#ffffff", relief="solid", bd=1)
        self.eval_text.pack(fill="both", expand=True, padx=8, pady=6)

    def _run_statistical_eval(self) -> None:
        """Executes ground-truth benchmark and scalability comparison."""
        try:
            seq_len = int(self.syn_len_entry.get().strip())
            num_sites = int(self.syn_sites_entry.get().strip())
            mut_rate = float(self.syn_mut_entry.get().strip())
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter valid numerical benchmark parameters.")
            return

        self.eval_text.delete("1.0", "end")
        self.eval_text.insert("end", f"=== STATISTICAL EVALUATION REPORT: {self.current_motif.name} ===\n")
        self.eval_text.insert("end", f"Generating synthetic DNA ({seq_len:,} bp) with {num_sites} implanted mutated sites...\n\n")

        dataset = GroundTruthGenerator.generate_dataset(
            motif_profile=self.current_motif,
            sequence_length=seq_len,
            num_sites=num_sites,
            mutation_rate=mut_rate,
            seed=42
        )

        pts = ThresholdAnalyzer.sweep_thresholds(dataset, self.current_motif, num_steps=7)
        opt_pt = ThresholdAnalyzer.find_optimal_threshold(pts)
        auc = ThresholdAnalyzer.compute_auc_roc(pts)

        self.eval_text.insert("end", "THRESHOLD SWEEP & CONFUSION MATRIX:\n")
        self.eval_text.insert("end", f"{'Threshold':<11} {'TP':<5} {'FP':<6} {'FN':<5} {'Precision':<11} {'Recall':<9} {'F1-Score':<10} {'MCC':<8}\n")
        self.eval_text.insert("end", "-" * 75 + "\n")
        for p in pts:
            cm = p.confusion_matrix
            self.eval_text.insert("end",
                f"{p.threshold:<11.2f} {cm.tp:<5} {cm.fp:<6} {cm.fn:<5} "
                f"{cm.precision:<11.3f} {cm.recall:<9.3f} {cm.f1_score:<10.3f} {cm.mcc:<8.3f}\n"
            )

        self.eval_text.insert("end", "-" * 75 + "\n")
        self.eval_text.insert("end", f"[+] Optimal Operating Threshold: {opt_pt.threshold:.2f} (Max F1 = {opt_pt.confusion_matrix.f1_score:.3f})\n")
        self.eval_text.insert("end", f"[+] Area Under ROC Curve (AUC-ROC): {auc:.3f}\n\n")

        # Multi-Engine Scalability
        self.eval_text.insert("end", "SCALABILITY & ENGINE COMPARISON (5,000 bp):\n")
        bench = ScalabilityProfiler.run_multi_engine_scalability_suite(self.current_motif, lengths=[5_000])
        self.eval_text.insert("end", f"{'Engine':<26} {'Length':<10} {'Hits':<8} {'Time (ms)':<12} {'Peak Mem':<12} {'Throughput'}\n")
        self.eval_text.insert("end", "-" * 80 + "\n")
        for r in bench:
            self.eval_text.insert("end",
                f"{r.engine_name:<26} {r.sequence_length:<10,d} {r.num_hits_found:<8} "
                f"{r.wall_time_ms:<12.3f} {r.peak_memory_kb:<10.2f} KB  {r.throughput_mbps:.2f} Mbp/s\n"
            )


def main() -> None:
    root = tk.Tk()
    app = MotifDetectionGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
