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

# --- Styling ---
st.markdown("""
    <style>
    .dna-seq { font-family: 'Courier New', monospace; font-size: 1.2em; line-height: 1.6; word-break: break-all; background: #f9f9f9; padding: 15px; border-radius: 8px; }
    .pam { background-color: #d1ffbd; color: #006400; padding: 2px 4px; border-radius: 3px; font-weight: bold; border-bottom: 2px solid #006400; }
    .mut { background-color: #ffcccc; color: #cc0000; padding: 2px 4px; border-radius: 3px; font-weight: bold; border-bottom: 2px solid #cc0000; }
    </style>
    """, unsafe_allow_html=True)

def resolve_gene_name(yeast_gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{yeast_gene}?content-type=application/json"
    response = requests.get(url)
    if response.status_code == 200 and len(response.json()) > 0: return response.json()[0]['id']
    return None

def get_gene_sequence(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    response = requests.get(url)
    return response.text.strip() if response.status_code == 200 else None

def find_cas9_sites(gene_sequence, aa_pos):
    nuc_pos = (aa_pos - 1) * 3
    window = 105
    start, end = max(0, nuc_pos - 50), min(len(gene_sequence), nuc_pos + 50)
    region = gene_sequence[start:end]
    sites = []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'reverse'})
    return sorted(sites, key=lambda x: abs(x['pos'] - nuc_pos))

# --- UI Setup ---
st.title("🧬 Yeast AutoOligo Helper")

with st.sidebar:
    gene_input = st.text_input("Gene", "PHO13")
    residue = st.number_input("Residue #", value=123, min_value=1)
    mutation_aa = st.text_input("Mutation", "R").upper()
    run = st.button("Generate Designs", type="primary")

if run:
    gene_id = resolve_gene_name(gene_input)
    if not gene_id: st.error("Gene not found.")
    else:
        full_seq = get_gene_sequence(gene_id)
        sites = find_cas9_sites(full_seq, residue)
        
        for i, site in enumerate(sites[:3]):
            with st.expander(f"Design {i+1} (Strand: {site['strand']})"):
                # 1. Oligo Design
                guide_20 = site['seq'][:-3] if site['strand'] == 'forward' else str(Seq(site['seq'][3:]).reverse_complement())
                oligo_a = f"gatc{guide_20}gttttagagctag"
                oligo_b = f"ctagctctaaaac{str(Seq(guide_20).reverse_complement())}"
                
                # 2. Repair Template Generation
                mut_idx = (residue - 1) * 3
                h_start, h_end = max(0, mut_idx - 40), min(len(full_seq), mut_idx + 43)
                repair_base = full_seq[h_start:h_end]
                
                # 3. Visual Highlighting Logic
                pam_start = site['pos'] - h_start
                mut_pos = mut_idx - h_start
                
                # Build HTML string for display
                html_output = '<div class="dna-seq">'
                for idx, char in enumerate(repair_base):
                    if idx in range(mut_pos, mut_pos+3):
                        html_output += f'<span class="mut">{char}</span>'
                    elif idx in range(pam_start + (20 if site['strand']=='forward' else 0), pam_start + (23 if site['strand']=='forward' else 3)):
                        html_output += f'<span class="pam">{char}</span>'
                    else:
                        html_output += char
                html_output += '</div>'
                
                st.markdown("### 🔬 Sequence Verification")
                st.markdown(html_output, unsafe_allow_html=True)
                st.caption("Green = PAM | Red = Mutation Site")
                
                # 4. Master Copy Block (Everything combined)
                st.markdown("### 📋 Copy All Sequences")
                master_text = f"--- DESIGN {i+1} ---\nOligo A: {oligo_a}\nOligo B: {oligo_b}\nRepair Template: {repair_base}"
                st.code(master_text, language="text")
