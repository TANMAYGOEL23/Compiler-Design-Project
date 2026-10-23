"""
Script to build the updated Microsoft Word (.docx) report documenting the
20-Motif, 20-Genome, Custom Sequence Appending, and Interactive GUI Architecture.
Generates an authentic WordprocessingML (.docx) package using standard zipfile.
"""

import zipfile
import os
import xml.sax.saxutils as saxutils

def escape(text: str) -> str:
    return saxutils.escape(text)

def build_report_docx(output_path: str) -> None:
    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
</Types>"""

    package_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>"""

    doc_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>"""

    styles_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:docDefaults>
    <w:rPrDefault>
      <w:rPr>
        <w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:cs="Calibri"/>
        <w:sz w:val="22"/>
        <w:color w:val="1F2937"/>
      </w:rPr>
    </w:rPrDefault>
    <w:pPrDefault>
      <w:pPr>
        <w:spacing w:line="260" w:lineRule="auto" w:after="120"/>
      </w:pPr>
    </w:pPrDefault>
  </w:docDefaults>
  
  <w:style w:type="paragraph" w:styleId="ReportTitle">
    <w:name w:val="Report Title"/>
    <w:rPr>
      <w:rFonts w:ascii="Arial" w:hAnsi="Arial"/>
      <w:b/>
      <w:sz w:val="42"/>
      <w:color w:val="0F172A"/>
    </w:rPr>
    <w:pPr>
      <w:spacing w:before="120" w:after="80"/>
      <w:jc w:val="center"/>
    </w:pPr>
  </w:style>

  <w:style w:type="paragraph" w:styleId="ReportSubtitle">
    <w:name w:val="Report Subtitle"/>
    <w:rPr>
      <w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/>
      <w:b/>
      <w:sz w:val="24"/>
      <w:color w:val="2563EB"/>
    </w:rPr>
    <w:pPr>
      <w:spacing w:after="200"/>
      <w:jc w:val="center"/>
    </w:pPr>
  </w:style>

  <w:style w:type="paragraph" w:styleId="Heading1">
    <w:name w:val="heading 1"/>
    <w:rPr>
      <w:rFonts w:ascii="Arial" w:hAnsi="Arial"/>
      <w:b/>
      <w:sz w:val="28"/>
      <w:color w:val="1E3A8A"/>
    </w:rPr>
    <w:pPr>
      <w:spacing w:before="240" w:after="100"/>
    </w:pPr>
  </w:style>

  <w:style w:type="paragraph" w:styleId="Heading2">
    <w:name w:val="heading 2"/>
    <w:rPr>
      <w:rFonts w:ascii="Arial" w:hAnsi="Arial"/>
      <w:b/>
      <w:sz w:val="24"/>
      <w:color w:val="2563EB"/>
    </w:rPr>
    <w:pPr>
      <w:spacing w:before="180" w:after="60"/>
    </w:pPr>
  </w:style>
</w:styles>"""

    body_parts = []

    def p(text: str, bold: bool = False, italic: bool = False, color: str = None, sz: int = None, before: int = 0, after: int = 100, align: str = "left") -> str:
        r_pr = ""
        if bold: r_pr += "<w:b/>"
        if italic: r_pr += "<w:i/>"
        if color: r_pr += f'<w:color w:val="{color}"/>'
        if sz: r_pr += f'<w:sz w:val="{sz}"/>'
        p_pr = f'<w:spacing w:before="{before}" w:after="{after}"/><w:jc w:val="{align}"/>'
        return f'<w:p><w:pPr>{p_pr}</w:pPr><w:r><w:rPr>{r_pr}</w:rPr><w:t xml:space="preserve">{escape(text)}</w:t></w:r></w:p>'

    def heading(text: str, level: int = 1) -> str:
        style = f"Heading{level}"
        return f'<w:p><w:pPr><w:pStyle w:val="{style}"/></w:pPr><w:r><w:t>{escape(text)}</w:t></w:r></w:p>'

    def bullet(text: str, bold_prefix: str = "") -> str:
        r_content = ""
        if bold_prefix:
            r_content += f'<w:r><w:rPr><w:b/><w:color w:val="1E3A8A"/></w:rPr><w:t xml:space="preserve">{escape(bold_prefix)} </w:t></w:r>'
        r_content += f'<w:r><w:t xml:space="preserve">{escape(text)}</w:t></w:r>'
        return f'<w:p><w:pPr><w:ind w:left="280"/><w:spacing w:before="30" w:after="50"/></w:pPr><w:r><w:rPr><w:color w:val="2563EB"/><w:b/></w:rPr><w:t>• </w:t></w:r>{r_content}</w:p>'

    def screenshot_placeholder(title: str, hint: str) -> str:
        out = []
        out.append(f'<w:p><w:pPr><w:spacing w:before="120" w:after="40"/></w:pPr><w:r><w:rPr><w:b/><w:color w:val="0F172A"/><w:sz w:val="22"/></w:rPr><w:t>{escape(title)}</w:t></w:r></w:p>')
        
        box_xml = f'''<w:tbl>
  <w:tblPr>
    <w:tblW w:w="9400" w:type="dxa"/>
    <w:tblBorders>
      <w:top w:val="dashed" w:sz="12" w:space="0" w:color="94A3B8"/>
      <w:left w:val="dashed" w:sz="12" w:space="0" w:color="94A3B8"/>
      <w:bottom w:val="dashed" w:sz="12" w:space="0" w:color="94A3B8"/>
      <w:right w:val="dashed" w:sz="12" w:space="0" w:color="94A3B8"/>
    </w:tblBorders>
    <w:tblCellMar>
      <w:top w:w="240" w:type="dxa"/><w:bottom w:w="240" w:type="dxa"/>
      <w:left w:w="200" w:type="dxa"/><w:right w:w="200" w:type="dxa"/>
    </w:tblCellMar>
  </w:tblPr>
  <w:tr>
    <w:tc>
      <w:tcPr><w:tcW w:w="9400" w:type="dxa"/><w:shd w:val="clear" w:color="auto" w:fill="F8FAFC"/></w:tcPr>
      <w:p>
        <w:pPr><w:jc w:val="center"/><w:spacing w:before="300" w:after="80"/></w:pPr>
        <w:r><w:rPr><w:b/><w:color w:val="2563EB"/><w:sz w:val="24"/></w:rPr><w:t>📷 [ PASTE SCREENSHOT HERE ]</w:t></w:r>
      </w:p>
      <w:p>
        <w:pPr><w:jc w:val="center"/><w:spacing w:before="40" w:after="300"/></w:pPr>
        <w:r><w:rPr><w:i/><w:color w:val="64748B"/><w:sz w:val="20"/></w:rPr><w:t>{escape(hint)}</w:t></w:r>
      </w:p>
    </w:tc>
  </w:tr>
</w:tbl>'''
        out.append(box_xml)
        return "".join(out)

    def team_table() -> str:
        headers = ["S.No.", "Team Member Name", "Registration / Roll No.", "Project Role & Module Contribution"]
        widths = [800, 2600, 2400, 3600]
        rows = [
            ["1", "Tanmay Goel (Leader)", "[ Reg / Roll No ]", "Module 1: 20-Motif Modeling, Information Theory & Staden DP"],
            ["2", "[ Team Member 2 Name ]", "[ Reg / Roll No ]", "Module 2: Formal Automata Engines, 2D Mismatch & VM"],
            ["3", "[ Team Member 3 Name ]", "[ Reg / Roll No ]", "Module 3: 20-Genome Registry, DNA Appender & Dual-Strand Scanner"],
            ["4", "[ Team Member 4 Name ]", "[ Reg / Roll No ]", "Module 4: Ground-Truth Evaluation, ROC Curves & Profiler"]
        ]
        
        out = ['<w:tbl><w:tblPr><w:tblW w:w="9400" w:type="dxa"/><w:tblBorders><w:top w:val="single" w:sz="8" w:space="0" w:color="CBD5E1"/><w:left w:val="single" w:sz="8" w:space="0" w:color="CBD5E1"/><w:bottom w:val="single" w:sz="8" w:space="0" w:color="CBD5E1"/><w:right w:val="single" w:sz="8" w:space="0" w:color="CBD5E1"/><w:insideH w:val="single" w:sz="4" w:space="0" w:color="E2E8F0"/><w:insideV w:val="single" w:sz="4" w:space="0" w:color="E2E8F0"/></w:tblBorders><w:tblCellMar><w:top w:w="100" w:type="dxa"/><w:bottom w:w="100" w:type="dxa"/><w:left w:w="140" w:type="dxa"/><w:right w:w="140" w:type="dxa"/></w:tblCellMar></w:tblPr>']
        out.append('<w:tr><w:trPr><w:tblHeader/></w:trPr>')
        for idx, h in enumerate(headers):
            out.append(f'<w:tc><w:tcPr><w:tcW w:w="{widths[idx]}" w:type="dxa"/><w:shd w:val="clear" w:color="auto" w:fill="1E3A8A"/></w:tcPr><w:p><w:pPr><w:spacing w:before="60" w:after="60"/></w:pPr><w:r><w:rPr><w:b/><w:color w:val="FFFFFF"/><w:sz w:val="20"/></w:rPr><w:t>{escape(h)}</w:t></w:r></w:p></w:tc>')
        out.append('</w:tr>')

        for r_idx, row in enumerate(rows):
            bg = "F8FAFC" if r_idx % 2 == 1 else "FFFFFF"
            out.append('<w:tr>')
            for idx, cell in enumerate(row):
                out.append(f'<w:tc><w:tcPr><w:tcW w:w="{widths[idx]}" w:type="dxa"/><w:shd w:val="clear" w:color="auto" w:fill="{bg}"/></w:tcPr><w:p><w:pPr><w:spacing w:before="50" w:after="50"/></w:pPr><w:r><w:rPr><w:sz w:val="20"/><w:color w:val="334155"/></w:rPr><w:t>{escape(cell)}</w:t></w:r></w:p></w:tc>')
            out.append('</w:tr>')

        out.append('</w:tbl>')
        return "".join(out)

    # --- TITLE & SUBTITLE ---
    body_parts.append('<w:p><w:pPr><w:pStyle w:val="ReportTitle"/></w:pPr><w:r><w:t>Automata-Based Probabilistic DNA Motif Detection</w:t></w:r></w:p>')
    body_parts.append('<w:p><w:pPr><w:pStyle w:val="ReportSubtitle"/></w:pPr><w:r><w:t>Compiler Design &amp; Theory of Computation Case Study</w:t></w:r></w:p>')

    # --- TEAM MEMBERS BLOCK ---
    body_parts.append(heading("Project Team Details", 2))
    body_parts.append(team_table())
    body_parts.append(p(""))

    # --- SECTION 1: HOW THIS IS ABOUT COMPILER DESIGN ---
    body_parts.append(heading("1. Compiler Design Connection (Theory & Concepts)", 1))
    body_parts.append(p("The fundamental goal of a compiler's Lexical Analyzer (Lexer) is pattern recognition over a finite alphabet. In a standard compiler, regular expressions are compiled into Deterministic Finite Automata (DFAs) that consume stream characters and classify tokens through state transitions."))
    body_parts.append(p("This project extends the classical compiler theory into bioinformatics by recognizing variable DNA patterns using Weighted Finite Automata (WFAs) and Mismatch-Aware 2D Automata:"))
    
    body_parts.append(bullet("DNA nucleotide alphabet Σ = {A, C, G, T} mapped to compiler input tokens.", "Finite Alphabet (Σ):"))
    body_parts.append(bullet("Target biological motifs (e.g. GATA1, p53, CTCF binding sites) treated as lexical tokens.", "Lexeme / Token: "))
    body_parts.append(bullet("States q0 → q1 → ... → qL, where each state represents consuming one position of the motif window.", "Automaton State Chain:"))
    body_parts.append(bullet("Instead of binary accept/reject (G → next), each transition carries a log-odds weight derived from position probability matrices: G / weight → next state.", "Weighted Transitions (W):"))
    body_parts.append(bullet("Tracks position i and consumed mismatch budget k in state space q(i, k), enforcing error tolerance constraints (0 ≤ k ≤ K_max).", "2D State Space q(i, k):"))
    body_parts.append(bullet("Pruning search branches when remaining maximum suffix weight cannot reach the threshold, mirroring compiler dead-code elimination.", "Lookahead Pruning:"))

    # --- SECTION 2: WHAT THE PROJECT ACHIEVES ---
    body_parts.append(heading("2. What the Project Achieves", 1))
    body_parts.append(bullet("Traditional regex/exact string matching fails on biological motifs because transcription factors tolerate mutations. The WFA scores sequence affinity quantitatively.", "Solves Exact Match Rigidity:"))
    body_parts.append(bullet("Expanded from a single motif to 20 biologically verified transcription factors (JASPAR profiles) spanning 6 bp to 20 bp.", "20-Motif Bank Registry:"))
    body_parts.append(bullet("Provides 20 preloaded genomic loci across human chromosomes and model organisms (Mouse, Zebrafish, Yeast, Viruses, Bacteria).", "20-Genome Reference Registry:"))
    body_parts.append(bullet("Allows users to enter new custom DNA sequences or append additional sequence bases directly to any selected genome.", "Custom DNA Sequence Appending:"))
    body_parts.append(bullet("Native Tkinter desktop application (gui_app.py) featuring interactive motif selection, genome viewers, canvas automaton graph drawing, and result tables.", "Modern Interactive Python UI:"))

    # --- SECTION 3: SYSTEM RESULTS & SCREENSHOT WALKTHROUGH ---
    body_parts.append(heading("3. Execution Results & Screenshot Walkthrough", 1))

    # --- SCREENSHOT 1 ---
    body_parts.append(screenshot_placeholder(
        "Screenshot 1: 20-Motif Explorer & PSSM Calculation (Module 1 / GUI Tab 1)",
        "Insert screenshot showing Motif Selector (GATA1 / p53 / CTCF), PFM matrix, PSSM Log-Odds, Information Content (bits), and Staden DP P-value inversion."
    ))
    body_parts.append(p("Minimal Result Explanation (Screenshot 1):", bold=True, color="1E3A8A"))
    body_parts.append(bullet("Converts raw JASPAR counts into smoothed probabilities using pseudocount α = 0.1.", "Probability Conversion:"))
    body_parts.append(bullet("Position 4 in GATA1 shows maximum conservation (1.55 bits, strong T preference), yielding Total Information Content = 5.904 bits.", "Information Content:"))
    body_parts.append(bullet("Staden DP discrete convolution calculates exact analytical P-values and inverts them to score thresholds.", "P-Value Threshold Inversion:"))
    body_parts.append(p(""))

    # --- SCREENSHOT 2 ---
    body_parts.append(screenshot_placeholder(
        "Screenshot 2: Interactive Automata Canvas Visualizer & Step Tracer (Module 2 / GUI Tab 2)",
        "Insert screenshot showing Canvas state graph q0->q1->...->qL, Step-by-step diagnostic trace for 'AGATAG', and Bytecode VM disassembly."
    ))
    body_parts.append(p("Minimal Result Explanation (Screenshot 2):", bold=True, color="1E3A8A"))
    body_parts.append(bullet("Visualizes states q0 through qL with color coding (Start in green, Accept in double-circle blue, active trace states in yellow).", "Automaton Graph Visualizer:"))
    body_parts.append(bullet("Demonstrates state walk q0→q1→...→q6 with weights (+0.89, +1.74, +1.73, +1.90, +1.64, -2.68) yielding Total Score = +5.221.", "1D-WFA Execution Trace:"))
    body_parts.append(bullet("Demonstrates linear opcode instructions (READ_BASE, ADD_WEIGHT, ACCEPT_HIT) generated for the Automaton Virtual Machine.", "Bytecode VM Disassembler:"))
    body_parts.append(p(""))

    # --- SCREENSHOT 3 ---
    body_parts.append(screenshot_placeholder(
        "Screenshot 3: 20-Genome Registry & Custom Sequence Appender (Module 3 / GUI Tab 1)",
        "Insert screenshot showing Genome dropdown, Sequence Viewer, and 'Add / Append DNA Sequence' action with success banner."
    ))
    body_parts.append(p("Minimal Result Explanation (Screenshot 3):", bold=True, color="1E3A8A"))
    body_parts.append(bullet("User can select from 20 authentic genomic loci (Human Chr11 HBB, Chr17 TP53, Chr8 MYC, Mouse Hbb, Zebrafish gata1a, etc.).", "20-Genome Registry:"))
    body_parts.append(bullet("Appends additional DNA bases (e.g. 50 bp regulatory fragment) directly to the active genome, instantly updating coordinates and length.", "Live DNA Sequence Appending:"))
    body_parts.append(bullet("Validates and cleans input DNA, stripping invalid characters and computing updated regional GC content.", "Sequence Sanitization:"))
    body_parts.append(p(""))

    # --- SCREENSHOT 4 ---
    body_parts.append(screenshot_placeholder(
        "Screenshot 4: Dual-Strand Genomic Scanner & Results Table (Module 3 / GUI Tab 3)",
        "Insert screenshot showing Dual-Strand (+ and -) results Treeview, surrounding genomic context preview, and BED/GFF3 export buttons."
    ))
    body_parts.append(p("Minimal Result Explanation (Screenshot 4):", bold=True, color="1E3A8A"))
    body_parts.append(bullet("Scans both forward (+) and reverse-complement (-) strands with exact chromosome coordinate mapping.", "Dual-Strand Coordinate Pipeline:"))
    body_parts.append(bullet("Interactive Treeview displays hit start/end, strand, matched sequence, score, p-value, and mismatch count.", "Results Table Inspector:"))
    body_parts.append(bullet("Clicking any hit row extracts and highlights the surrounding ±20 bp genomic context.", "Genomic Context Viewer:"))
    body_parts.append(p(""))

    # --- SCREENSHOT 5 ---
    body_parts.append(screenshot_placeholder(
        "Screenshot 5: Ground-Truth Evaluation & Scalability Profiling (Module 4 / GUI Tab 4)",
        "Insert screenshot showing Confusion Matrix (TP, FP, FN, TN), Optimal Operating Threshold, and Multi-Engine Scalability comparison table."
    ))
    body_parts.append(p("Minimal Result Explanation (Screenshot 5):", bold=True, color="1E3A8A"))
    body_parts.append(bullet("Synthetic sequence with known implanted mutated motifs evaluates Precision, Recall, F1, and MCC without circular bias.", "Ground-Truth Benchmark:"))
    body_parts.append(bullet("Systematic threshold sweep identifies optimal score cutoff maximizing F1-score.", "Optimal Threshold Detector:"))
    body_parts.append(bullet("Direct PSSM and Pruned WFA achieve high throughput (~0.78 Mbp/s) with minimal memory footprint (< 55 KB).", "Scalability & Profiler:"))

    # --- SECTION 4: SUMMARY & VIVA TAKEAWAYS ---
    body_parts.append(heading("4. Summary & Key Takeaways for Review / Viva", 1))
    body_parts.append(bullet("A standard compiler asks 'Is this token valid?' (Boolean accept/reject). Our Weighted Finite Automaton asks 'How strongly does this input match the pattern?' (Probabilistic log-odds scoring).", "Why Automata?"))
    body_parts.append(bullet("Exact regex matches only a few sites. The probabilistic WFA identifies hundreds of candidate regulatory sites, capturing biological mutation tolerance.", "Why WFA over DFA?"))
    body_parts.append(bullet("The 2D state space q(i, k) explicitly tracks position and mutation budget, bridging compiler state-space tracking with biological error tolerance.", "Why 2D Mismatch Automaton?"))
    body_parts.append(bullet("Preloaded 20 transcription factors and 20 reference genomes, with an interactive desktop GUI and sequence appender.", "Rich Biological Breadth:"))

    body_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:body>
    {''.join(body_parts)}
    <w:sectPr>
      <w:pgSz w:w="12240" w:h="15840"/>
      <w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440"/>
    </w:sectPr>
  </w:body>
</w:document>"""

    with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as docx:
        docx.writestr('[Content_Types].xml', content_types)
        docx.writestr('_rels/.rels', package_rels)
        docx.writestr('word/_rels/document.xml.rels', doc_rels)
        docx.writestr('word/styles.xml', styles_xml)
        docx.writestr('word/document.xml', body_xml)

    print(f"Successfully generated Report Word Document: {output_path}")

if __name__ == "__main__":
    out_file = "/Users/tanmaygoel/Documents/PROJECTS/Compiler Design/Compiler_Design_DNA_Motif_Detection_Report.docx"
    build_report_docx(out_file)
