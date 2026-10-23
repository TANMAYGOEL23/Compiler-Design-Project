"""
Automata-Based DNA Motif Detection: Interactive Streamlit Web Application
========================================================================
A modern, aesthetic bioinformatics web application implementing formal compiler design:
Weighted Finite Automata (WFA), Mismatch-Aware 2D Automata, Staden DP P-values,
20 JASPAR Transcription Factors, 20 Reference Genomes, and Custom DNA Appending.
"""

from __future__ import annotations

import json
import os
import sys
import time
import pandas as pd
import numpy as np
import streamlit as st

# Import the 4 underlying compiler modules
from module1_motif_models import DNA_BASES, BASE_TO_INDEX, MotifProfile, MotifBank
from module2_automata_engine import (
    LinearWFA,
    MismatchAwareWFA,
    ClassicalExactDFA,
    AutomatonBytecodeCompiler,
    AutomatonVisualLayout
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
    ScalabilityProfiler,
    MultiMotifBenchmarkSuite
)


# ============================================================================
# PAGE CONFIGURATION & CUSTOM MODERN CSS
# ============================================================================

st.set_page_config(
    page_title="Automata DNA Motif Scanner | Compiler Design",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

CUSTOM_CSS = """
<style>
    /* Global Typography & Palette */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
    
    html, body, [class*="css"] {
        font-family: 'Inter', sans-serif;
    }
    
    code, pre {
        font-family: 'JetBrains Mono', monospace !important;
    }
    
    /* Header Gradient Banner */
    .hero-banner {
        background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 50%, #2563eb 100%);
        padding: 24px 30px;
        border-radius: 14px;
        color: white;
        margin-bottom: 25px;
        box-shadow: 0 4px 15px rgba(15, 23, 42, 0.15);
    }
    .hero-title {
        font-size: 26px;
        font-weight: 700;
        margin: 0;
        color: #ffffff;
    }
    .hero-subtitle {
        font-size: 14px;
        color: #93c5fd;
        margin-top: 6px;
        margin-bottom: 0;
    }
    
    /* Metric Cards */
    .metric-card {
        background: #ffffff;
        border: 1px solid #e2e8f0;
        border-radius: 10px;
        padding: 16px;
        text-align: center;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .metric-val {
        font-size: 24px;
        font-weight: 700;
        color: #1e3a8a;
    }
    .metric-lbl {
        font-size: 12px;
        color: #64748b;
        text-transform: uppercase;
        font-weight: 600;
        letter-spacing: 0.5px;
    }
    
    /* Pill Badges */
    .badge-exact {
        background-color: #dcfce7;
        color: #15803d;
        padding: 3px 8px;
        border-radius: 12px;
        font-size: 11px;
        font-weight: 600;
    }
    .badge-wfa {
        background-color: #dbeafe;
        color: #1d4ed8;
        padding: 3px 8px;
        border-radius: 12px;
        font-size: 11px;
        font-weight: 600;
    }
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)


# ============================================================================
# STATE INITIALIZATION
# ============================================================================

if "motif_bank" not in st.session_state:
    st.session_state.motif_bank = MotifBank()
if "genome_registry" not in st.session_state:
    st.session_state.genome_registry = GenomeRegistry()
if "selected_motif_name" not in st.session_state:
    st.session_state.selected_motif_name = "GATA1"
if "selected_genome_id" not in st.session_state:
    st.session_state.selected_genome_id = "GENOME_HUMAN_HBB_CHR11"
if "hits" not in st.session_state:
    st.session_state.hits = []
if "scan_time" not in st.session_state:
    st.session_state.scan_time = 0.0


# ============================================================================
# SIDEBAR NAVIGATION & QUICK CONTROLS
# ============================================================================

with st.sidebar:
    st.markdown("### 🧬 DNA Compiler Suite")
    st.caption("Compiler Design & Theory of Computation Case Study")
    
    app_mode = st.radio(
        "Navigation Menu:",
        [
            "🔬 1. Motif Bank & PSSM Models",
            "🌍 2. Genome Registry & DNA Editor",
            "⚙️ 3. Automata Graph & VM Tracer",
            "⚡ 4. High-Throughput Genomic Scanner",
            "📊 5. Evaluation & Scalability Benchmarks"
        ],
        index=0
    )
    
    st.divider()
    
    # Active Motif Selector
    motif_names = st.session_state.motif_bank.list_motif_names()
    default_m_idx = motif_names.index("GATA1") if "GATA1" in motif_names else 0
    selected_m = st.selectbox(
        "Active Transcription Factor (20 Motifs):",
        motif_names,
        index=default_m_idx,
        help="Choose from 20 biologically verified transcription factor binding profiles (JASPAR)."
    )
    st.session_state.selected_motif_name = selected_m
    active_motif = st.session_state.motif_bank.get(selected_m)
    
    # Active Genome Selector
    all_genomes = st.session_state.genome_registry.get_all_records()
    genome_labels = [f"{g.locus_name} ({g.organism})" for g in all_genomes]
    genome_id_map = {f"{g.locus_name} ({g.organism})": g.sequence_id for g in all_genomes}
    
    selected_g_label = st.selectbox(
        "Active Genomic Target (20 Genomes):",
        genome_labels,
        index=0,
        help="Choose from 20 preloaded reference genomic sequences spanning human and model organisms."
    )
    active_genome_id = genome_id_map[selected_g_label]
    st.session_state.selected_genome_id = active_genome_id
    active_genome = st.session_state.genome_registry.get(active_genome_id)
    
    st.divider()
    st.info(f"**Target:** `{active_motif.name}` ({active_motif.length} bp)\n\n**Locus:** `{active_genome.locus_name}` ({active_genome.length:,} bp)")


# ============================================================================
# TOP BANNER
# ============================================================================

st.markdown(f"""
<div class="hero-banner">
    <h1 class="hero-title">Automata-Based Probabilistic DNA Motif Detection</h1>
    <p class="hero-subtitle">Compiler Lexical Analysis with Weighted Finite Automata (WFA) & 2D Mismatch State Space q(i, k)</p>
</div>
""", unsafe_allow_html=True)


# ============================================================================
# TAB 1: MOTIF BANK & PSSM MODELS
# ============================================================================

if "1. Motif Bank" in app_mode:
    st.subheader(f"🧬 Motif Profile: {active_motif.name} (JASPAR ID: {active_motif.motif_id})")
    st.write(f"**Organism:** {active_motif.organism} | **Target Function:** {active_motif.biological_function}")
    
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.markdown(f"""<div class="metric-card"><div class="metric-val">{active_motif.length} bp</div><div class="metric-lbl">Motif Length</div></div>""", unsafe_allow_html=True)
    with c2:
        st.markdown(f"""<div class="metric-card"><div class="metric-val">{active_motif.consensus_sequence}</div><div class="metric-lbl">Strict Consensus</div></div>""", unsafe_allow_html=True)
    with c3:
        st.markdown(f"""<div class="metric-card"><div class="metric-val">{active_motif.iupac_consensus}</div><div class="metric-lbl">IUPAC Consensus</div></div>""", unsafe_allow_html=True)
    with c4:
        st.markdown(f"""<div class="metric-card"><div class="metric-val">{active_motif.total_information_content:.2f} bits</div><div class="metric-lbl">Total Info Content</div></div>""", unsafe_allow_html=True)
        
    st.markdown("#### Position-Specific Scoring Matrix (PSSM Log-Odds Weights in Bits)")
    st.caption("Positive values indicate preferential nucleotide binding relative to genomic background ($Q=0.25$).")
    pssm_df = pd.DataFrame(
        active_motif.pssm,
        index=["A", "C", "G", "T"],
        columns=[f"Pos {j+1}" for j in range(active_motif.length)]
    )
    # Render with safe rounding without external matplotlib dependency
    st.dataframe(pssm_df.round(3), use_container_width=True)
    
    # Information Content Bar Chart
    st.markdown("#### 📊 Information Content Distribution Across Motif ($R_{seq}$ bits)")
    info_chart_df = pd.DataFrame({
        "Position": [f"Pos {info.position + 1} ({info.dominant_base})" for info in active_motif.information_per_position],
        "Information Content (bits)": [round(info.information_content, 3) for info in active_motif.information_per_position]
    }).set_index("Position")
    st.bar_chart(info_chart_df)
    
    st.markdown("#### Position Breakdown & Shannon Entropy")
    pos_data = []
    for info in active_motif.information_per_position:
        pos_data.append({
            "Position": info.position + 1,
            "Dominant Base": info.dominant_base,
            "IUPAC Symbol": info.consensus_symbol,
            "Shannon Entropy (bits)": round(info.shannon_entropy, 3),
            "Info Content R_seq (bits)": round(info.information_content, 3),
            "Relative Entropy D_KL": round(info.relative_entropy, 3),
            "Prob(A)": round(info.probabilities["A"], 3),
            "Prob(C)": round(info.probabilities["C"], 3),
            "Prob(G)": round(info.probabilities["G"], 3),
            "Prob(T)": round(info.probabilities["T"], 3),
        })
    st.dataframe(pd.DataFrame(pos_data), use_container_width=True)
    
    # Staden DP Inversion
    st.markdown("#### 📐 Staden Dynamic Programming P-Value Distribution")
    scores, pvals = active_motif.calculate_pvalue_score_distribution(granularity=0.05)
    sample_indices = np.linspace(0, len(scores)-1, min(60, len(scores)), dtype=int)
    chart_df = pd.DataFrame({
        "Score Threshold": [round(float(scores[i]), 2) for i in sample_indices],
        "P-Value P(Score >= T)": [float(pvals[i]) for i in sample_indices]
    }).set_index("Score Threshold")
    st.line_chart(chart_df)
    
    t_opt = active_motif.pvalue_to_threshold(1e-3)
    st.caption(f"Staden DP analytical threshold for statistical significance P <= 0.001 is **{t_opt:.3f}**.")


# ============================================================================
# TAB 2: GENOME REGISTRY & CUSTOM SEQUENCE APPENDER
# ============================================================================

elif "2. Genome Registry" in app_mode:
    st.subheader(f"🌍 Genomic Target: {active_genome.locus_name}")
    st.write(f"**Description:** {active_genome.description}")
    
    g1, g2, g3, g4 = st.columns(4)
    with g1:
        st.markdown(f"""<div class="metric-card"><div class="metric-val">{active_genome.length:,} bp</div><div class="metric-lbl">Total Length</div></div>""", unsafe_allow_html=True)
    with g2:
        st.markdown(f"""<div class="metric-card"><div class="metric-val">{active_genome.get_gc_content():.2%}</div><div class="metric-lbl">GC Content</div></div>""", unsafe_allow_html=True)
    with g3:
        st.markdown(f"""<div class="metric-card"><div class="metric-val">{active_genome.start_coordinate:,}</div><div class="metric-lbl">Start Coordinate</div></div>""", unsafe_allow_html=True)
    with g4:
        st.markdown(f"""<div class="metric-card"><div class="metric-val">{active_genome.end_coordinate:,}</div><div class="metric-lbl">End Coordinate</div></div>""", unsafe_allow_html=True)
        
    st.markdown("#### Sequence Viewer (First 600 bp)")
    preview = active_genome.sequence[:600] + ("..." if active_genome.length > 600 else "")
    st.code(preview, language="text")
    
    st.divider()
    
    # Custom DNA Sequence Adder / Appender
    st.markdown("### ➕ Add or Append DNA Sequences")
    st.write("You can paste arbitrary nucleotide sequences to append them to the current genome or register a brand-new locus.")
    
    append_mode = st.radio(
        "Action Mode:",
        ["Append extra sequence to current selected genome", "Create a brand new custom genome"],
        horizontal=True
    )
    
    seq_input = st.text_area(
        "Enter or paste DNA sequence (A, C, G, T):",
        placeholder="e.g. AGATAAGGGGCGGGGGCACGGCACGGCACGGTTTTATAAAATTTTTCACGTG",
        height=120
    )
    
    c_btn1, c_btn2 = st.columns([1, 4])
    with c_btn1:
        if st.button("Apply Sequence", type="primary"):
            cleaned = clean_dna_sequence(seq_input)
            if not cleaned:
                st.error("Please enter valid nucleotide characters (A, C, G, T).")
            else:
                if "Append" in append_mode:
                    prev_len = active_genome.length
                    new_len = active_genome.append_sequence(cleaned)
                    st.success(f"Successfully appended {len(cleaned):,} bp to {active_genome.locus_name}! New length: {new_len:,} bp.")
                    st.rerun()
                else:
                    new_id = f"CUSTOM_GENOME_{int(time.time())}"
                    rec = st.session_state.genome_registry.add_custom_genome(
                        sequence_id=new_id,
                        sequence_data=cleaned,
                        description="User Added Custom Genomic Sequence"
                    )
                    st.success(f"Registered new custom genome {new_id} ({len(cleaned):,} bp).")
                    st.rerun()
    with c_btn2:
        if st.button("Load Regulatory Fragment Example"):
            sample = "AGATAAGGGGCGGGGGCACGGCACGGCACGGTTTTATAAAATTTTTCACGTG"
            st.session_state.sample_dna = sample
            st.info(f"Sample DNA copied: `{sample}`. Paste it above to apply.")


# ============================================================================
# TAB 3: AUTOMATA GRAPH & VM TRACER
# ============================================================================

elif "3. Automata Graph" in app_mode:
    st.subheader(f"⚙️ Compiler State Automata: {active_motif.name}")
    st.write("Formal Automata Formulation: $M = (Q, \\Sigma, \\delta, W, q_0, F)$ where $\\delta(q_i, b) = q_{i+1}$ with weight $W(b, i) = \\text{PSSM}[b][i]$.")
    
    wfa = LinearWFA(active_motif)
    dfa = ClassicalExactDFA(f"{active_motif.name}_Exact", active_motif.iupac_consensus)
    
    # Draw Graphviz Diagram
    st.markdown("#### Automaton State Transition Chain")
    dot_lines = [
        "digraph WFA {",
        "  rankdir=LR;",
        '  node [fontname="Helvetica", shape=circle, style=filled, fillcolor="#f8fafc", color="#64748b"];',
        '  edge [fontname="Helvetica", fontsize=9, color="#94a3b8"];',
        f'  "q0" [fillcolor="#dcfce7", color="#16a34a", penwidth=2];',
        f'  "q{wfa.length}" [shape=doublecircle, fillcolor="#dbeafe", color="#2563eb", penwidth=2];'
    ]
    for pos in range(wfa.length):
        dominant = active_motif.consensus_sequence[pos]
        trans = wfa.transitions[wfa.states[pos]][dominant]
        dot_lines.append(f'  "q{pos}" -> "q{pos+1}" [label="{dominant} / {trans.weight:+.1f}", color="#2563eb", penwidth=2];')
    dot_lines.append("}")
    
    st.graphviz_chart("\n".join(dot_lines))
    
    st.divider()
    
    # Diagnostic Single-Window Tracer
    st.markdown("#### 🔍 Step-by-Step Diagnostic Window Tracer")
    test_cand = st.text_input(
        f"Enter candidate window ({active_motif.length} bp):",
        value=active_motif.consensus_sequence,
        max_chars=active_motif.length
    ).upper().strip()
    
    if len(test_cand) == active_motif.length:
        score, traces = wfa.trace_window(test_cand)
        is_exact = dfa.match_window(test_cand)
        
        c_sc1, c_sc2, c_sc3 = st.columns(3)
        with c_sc1:
            st.metric("Total Traversal Score", f"{score:.3f}")
        with c_sc2:
            st.metric("Status", "ACCEPTED" if score != float("-inf") else "REJECTED")
        with c_sc3:
            st.metric("Exact DFA Token Match", "YES" if is_exact else "NO")
            
        trace_data = []
        for t in traces:
            trace_data.append({
                "Step": t.step,
                "Source State": t.source_state.name,
                "Base Read": t.symbol_read,
                "Target State": t.target_state.name,
                "Transition Weight": f"{t.step_weight:+.3f}",
                "Running Score": f"{t.accumulated_score:+.3f}"
            })
        st.dataframe(pd.DataFrame(trace_data), use_container_width=True)
    else:
        st.warning(f"Please enter exactly {active_motif.length} bases.")


# ============================================================================
# TAB 4: GENOMIC SCANNER & RESULTS TABLE
# ============================================================================

elif "4. High-Throughput" in app_mode:
    st.subheader(f"⚡ High-Throughput Dual-Strand Scanner")
    st.write(f"Scanning target **{active_genome.locus_name}** ({active_genome.length:,} bp) for **{active_motif.name}** ({active_motif.length} bp) binding affinity.")
    
    # Scanner Parameters Form
    with st.expander("⚙️ Scanner Settings & Pruning Options", expanded=True):
        sc_col1, sc_col2, sc_col3, sc_col4 = st.columns(4)
        with sc_col1:
            rec_thresh = max(2.5, round(active_motif.max_possible_score * 0.55, 2))
            thresh_val = st.number_input("Score Threshold:", value=float(rec_thresh), step=0.5)
        with sc_col2:
            k_mismatch = st.selectbox("Mismatch Budget K:", [0, 1, 2, 3], index=2)
        with sc_col3:
            scan_dual = st.checkbox("Scan Both Strands (+ and -)", value=True)
        with sc_col4:
            use_prune = st.checkbox("Branch-and-Bound Pruning", value=True)
            
    if st.button("🚀 Execute WFA Genomic Scan", type="primary"):
        with st.spinner("Traversing Weighted Automaton across genomic sequence..."):
            pipeline = DualStrandScannerPipeline(active_motif, max_mismatches=k_mismatch)
            t0 = time.perf_counter()
            hits = pipeline.scan_genomic_record(
                active_genome,
                threshold=thresh_val,
                scan_both_strands=scan_dual,
                use_pruning=use_prune
            )
            t_el = (time.perf_counter() - t0) * 1000.0
            st.session_state.hits = hits
            st.session_state.scan_time = t_el
            
    hits = st.session_state.hits
    scan_time = st.session_state.scan_time
    
    if hits:
        fwd_hits = sum(1 for h in hits if h.strand == "+")
        rev_hits = sum(1 for h in hits if h.strand == "-")
        exact_hits = sum(1 for h in hits if h.is_exact_match)
        
        m1, m2, m3, m4, m5 = st.columns(5)
        with m1:
            st.markdown(f"""<div class="metric-card"><div class="metric-val">{len(hits)}</div><div class="metric-lbl">Total Candidate Hits</div></div>""", unsafe_allow_html=True)
        with m2:
            st.markdown(f"""<div class="metric-card"><div class="metric-val">{fwd_hits}</div><div class="metric-lbl">Forward (+) Strand</div></div>""", unsafe_allow_html=True)
        with m3:
            st.markdown(f"""<div class="metric-card"><div class="metric-val">{rev_hits}</div><div class="metric-lbl">Reverse (-) Strand</div></div>""", unsafe_allow_html=True)
        with m4:
            st.markdown(f"""<div class="metric-card"><div class="metric-val">{exact_hits}</div><div class="metric-lbl">Exact DFA Matches</div></div>""", unsafe_allow_html=True)
        with m5:
            st.markdown(f"""<div class="metric-card"><div class="metric-val">{scan_time:.2f} ms</div><div class="metric-lbl">Scan Time</div></div>""", unsafe_allow_html=True)
            
        st.markdown("#### Candidate Transcription Factor Binding Sites")
        
        # Display Table
        hits_table_data = []
        for idx, h in enumerate(hits, start=1):
            hits_table_data.append({
                "#": idx,
                "Start": h.start,
                "End": h.end,
                "Strand": h.strand,
                "Matched Sequence": h.matched_sequence,
                "WFA Score": round(h.score, 3),
                "P-Value": f"{h.p_value:.2e}",
                "Mismatches": h.mismatches,
                "Match Type": "Exact Match" if h.is_exact_match else "Probabilistic WFA"
            })
        hits_df = pd.DataFrame(hits_table_data)
        st.dataframe(hits_df, use_container_width=True, height=350)
        
        # Flanking Genomic Context Viewer
        st.markdown("#### 🔬 Inspect Local Sequence Context")
        sel_idx = st.number_input("Select Hit Row #:", min_value=1, max_value=len(hits), value=1, step=1)
        target_hit = hits[sel_idx - 1]
        
        loc_start = target_hit.start - active_genome.start_coordinate
        loc_end = target_hit.end - active_genome.start_coordinate
        full_s = active_genome.sequence
        flank = 20
        l_flank = full_s[max(0, loc_start - flank):loc_start]
        hit_s = full_s[loc_start:loc_end]
        r_flank = full_s[loc_end:min(len(full_s), loc_end + flank)]
        
        st.markdown(f"**Chromosome Coordinates:** `{target_hit.start:,}` - `{target_hit.end:,}` ({target_hit.strand} strand)")
        st.code(f"...{l_flank} [{hit_s}] {r_flank}...", language="text")
        
        # Export Actions
        st.markdown("#### 📥 Export Detected Binding Sites")
        d1, d2, d3, d4 = st.columns(4)
        
        csv_data = hits_df.to_csv(index=False).encode('utf-8')
        with d1:
            st.download_button("Download CSV", data=csv_data, file_name=f"{active_motif.name}_hits.csv", mime="text/csv")
            
        bed_lines = ["#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n"] + [h.to_bed_line() + "\n" for h in hits]
        with d2:
            st.download_button("Download BED6", data="".join(bed_lines), file_name=f"{active_motif.name}_hits.bed", mime="text/plain")
            
        gff_lines = ["##gff-version 3\n"] + [h.to_gff3_line() + "\n" for h in hits]
        with d3:
            st.download_button("Download GFF3", data="".join(gff_lines), file_name=f"{active_motif.name}_hits.gff3", mime="text/plain")
            
        json_data = json.dumps([h.to_dict() for h in hits], indent=2)
        with d4:
            st.download_button("Download JSON", data=json_data, file_name=f"{active_motif.name}_hits.json", mime="application/json")
    else:
        st.info("Click 'Execute WFA Genomic Scan' above to scan the active genomic sequence.")


# ============================================================================
# TAB 5: EVALUATION, ROC CURVES & SCALABILITY BENCHMARKS
# ============================================================================

elif "5. Evaluation" in app_mode:
    st.subheader(f"📊 Statistical Evaluation & Scalability Suite")
    st.write("Generates synthetic genomic datasets with known implanted mutations to compute un-biased Precision, Recall, and ROC curves.")
    
    col_e1, col_e2, col_e3 = st.columns(3)
    with col_e1:
        syn_len = st.number_input("Benchmark Sequence Length (bp):", value=15000, step=5000)
    with col_e2:
        syn_sites = st.number_input("Implanted Known True Sites:", value=35, step=5)
    with col_e3:
        syn_mut = st.slider("Mutation Noise Rate:", min_value=0.0, max_value=0.30, value=0.10, step=0.05)
        
    if st.button("Run Statistical Ground-Truth Evaluation", type="primary"):
        with st.spinner("Generating ground truth and sweeping score thresholds..."):
            dataset = GroundTruthGenerator.generate_dataset(
                motif_profile=active_motif,
                sequence_length=syn_len,
                num_sites=syn_sites,
                mutation_rate=syn_mut,
                seed=42
            )
            pts = ThresholdAnalyzer.sweep_thresholds(dataset, active_motif, num_steps=8)
            opt_pt = ThresholdAnalyzer.find_optimal_threshold(pts)
            auc = ThresholdAnalyzer.compute_auc_roc(pts)
            
            e_m1, e_m2, e_m3, e_m4 = st.columns(4)
            with e_m1:
                st.metric("Optimal Threshold", f"{opt_pt.threshold:.2f}")
            with e_m2:
                st.metric("Optimal F1-Score", f"{opt_pt.confusion_matrix.f1_score:.3f}")
            with e_m3:
                st.metric("Optimal Precision", f"{opt_pt.confusion_matrix.precision:.3f}")
            with e_m4:
                st.metric("AUC-ROC", f"{auc:.3f}")
                
            st.markdown("#### Confusion Matrix across Operating Thresholds")
            sweep_data = []
            for p in pts:
                cm = p.confusion_matrix
                sweep_data.append({
                    "Threshold": round(p.threshold, 2),
                    "TP": cm.tp,
                    "FP": cm.fp,
                    "FN": cm.fn,
                    "Precision": round(cm.precision, 3),
                    "Recall (Sensitivity)": round(cm.recall, 3),
                    "F1-Score": round(cm.f1_score, 3),
                    "MCC": round(cm.mcc, 3),
                    "FPR": round(cm.fpr, 4)
                })
            sw_df = pd.DataFrame(sweep_data)
            st.dataframe(sw_df, use_container_width=True)
            
            # ROC Plot
            st.markdown("#### Receiver Operating Characteristic (ROC) Curve")
            roc_chart = sw_df[["FPR", "Recall (Sensitivity)"]].set_index("FPR")
            st.line_chart(roc_chart)
            
    st.divider()
    
    st.markdown("#### Multi-Engine Computational Scalability Profiler (5,000 bp)")
    bench_results = ScalabilityProfiler.run_multi_engine_scalability_suite(active_motif, lengths=[5_000])
    b_data = []
    for r in bench_results:
        b_data.append({
            "Compiler Engine": r.engine_name,
            "Sequence Length": f"{r.sequence_length:,} bp",
            "Hits Found": r.num_hits_found,
            "Execution Time (ms)": round(r.wall_time_ms, 3),
            "Peak Memory (KB)": round(r.peak_memory_kb, 2),
            "Throughput (Mbp/s)": round(r.throughput_mbps, 3)
        })
    st.dataframe(pd.DataFrame(b_data), use_container_width=True)

