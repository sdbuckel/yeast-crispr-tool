import streamlit as st
import requests
import re
from Bio.Seq import Seq

# 1. Global Codon Table
CODON_TABLE = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    '*': ['TAA', 'TAG', 'TGA'],
}

# --- CSS for Pixel-Perfect Alignment ---
st.markdown("""
    <style>
    .align-table { 
        font-family: 'Courier New', monospace; 
        border-collapse: collapse; 
        margin: 20px 0;
        background-color: #ffffff;
    }
    .align-table td { 
        padding: 0px 1px; 
        text-align: center; 
        width: 14px; 
        font-size: 1.1em;
    }
    .label-cell { 
        text-align: right !important; 
        width: 100px !important; 
        color: #666; 
        font-size: 0.8em !important; 
        padding-right: 15px !important;
        border-right: 1px solid #eee;
    }
    .pam-site { background-color: #d1ffbd; color: #006400; font-weight: bold; }
    .mut-site { background-color: #ffcccc; color: #cc0000; font-weight: bold; }
    </style>
    """, unsafe_allow_html=True)

# [Functions resolve_gene_name and get_gene_sequence remain same]
def resolve_gene_name(yeast_gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{yeast_gene}?content-type=application/json"
    response = requests.get(url)
    return response.json()[0]['id'] if response.status_code == 200 and len(response.json()) > 0 else None

def get_gene_sequence(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    response = requests.get(url)
    return response.text.strip() if response.status_code == 200 else None

def get_mutated_template(original_seq, mut_idx, target_aa):
    codons = CODON_TABLE.get(target_aa, ["???"])
    new_codon = codons[0]
    seq_list = list(original_seq)
    seq_list[mut_idx:mut_idx+3] = list(new_codon)
    return "".join(seq_list)

# --- UI Setup ---
st.title("🧬 Yeast AutoOligo: Grid Proofreader")

with st.sidebar:
    gene_input = st.text_input("Gene", "PHO13").strip()
    residue = st.number_input("Residue #", value=123, min_value=1)
    mutation_aa = st.text_input("Mutation", "R").upper().strip()
    run = st.button("Generate & Align", type="primary")

if run:
    gene_id = resolve_gene_name(gene_input)
    if not gene_id: st.error("Gene not found.")
    else:
        full_seq = get_gene_sequence(gene_id)
        mut_idx_gene = (residue - 1) * 3
        # Adjust start to be on a codon boundary for cleaner translation
        start_idx = max(0, mut_idx_gene - 24)
        end_idx = min(len(full_seq), mut_idx_gene + 27)
        
        wt_dna = full_seq[start_idx:end_idx]
        mut_dna = get_mutated_template(wt_dna, mut_idx_gene - start_idx, mutation_aa)
        
        wt_aa = [str(Seq(wt_dna[i:i+3]).translate()) for i in range(0, len(wt_dna), 3)]
        mut_aa = [str(Seq(mut_dna[i:i+3]).translate()) for i in range(0, len(mut_dna), 3)]
        mut_rel = mut_idx_gene - start_idx

        # --- Construct Grid Alignment ---
        html = '<table class="align-table">'
        
        # 1. WT Amino Acid Row
        html += '<tr><td class="label-cell">WT AA</td>'
        for aa in wt_aa: html += f'<td colspan="3" style="color:#888">{aa}</td>'
        html += '</tr>'

        # 2. WT DNA Row
        html += '<tr><td class="label-cell">WT DNA</td>'
        for i, b in enumerate(wt_dna):
            style = ' class="mut-site"' if i in range(mut_rel, mut_rel+3) else ''
            html += f'<td{style}>{b}</td>'
        html += '</tr>'

        # 3. Mutant DNA Row
        html += '<tr><td class="label-cell">MUT DNA</td>'
        for i, b in enumerate(mut_dna):
            style = ' class="mut-site"' if i in range(mut_rel, mut_rel+3) else ''
            html += f'<td{style}>{b}</td>'
        html += '</tr>'

        # 4. Mutant Amino Acid Row
        html += '<tr><td class="label-cell">MUT AA</td>'
        for aa in mut_aa: html += f'<td colspan="3" style="color:#000; font-weight:bold;">{aa}</td>'
        html += '</tr>'
        
        html += '</table>'
        
        st.markdown("### 🔍 Pixel-Perfect Alignment")
        st.markdown(html, unsafe_allow_html=True)

        # Copy Blocks
        st.subheader("📋 Sequences")
        st.code(f"Sense: {mut_dna}\nComplement: {str(Seq(mut_dna).complement())}", language="text")
