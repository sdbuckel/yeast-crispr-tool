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

# --- CSS for Highlighting ---
st.markdown("""
    <style>
    .align-table { font-family: 'Courier New', monospace; border-collapse: collapse; margin: 20px 0; background-color: #ffffff; }
    .align-table td { padding: 2px 1px; text-align: center; width: 18px; font-size: 1.1em; border: 1px solid #f0f0f0; }
    .label-cell { text-align: right !important; width: 100px !important; color: #666; font-size: 0.8em !important; padding-right: 15px !important; border-right: 2px solid #ddd !important; }
    .pam-site { background-color: #d1ffbd; color: #006400; font-weight: bold; }
    .mut-site { background-color: #ffcccc; color: #cc0000; font-weight: bold; }
    .overlap-site { background-color: #ffe5b4; color: #8b4513; font-weight: bold; } /* Orange for Mut+PAM overlap */
    </style>
    """, unsafe_allow_html=True)

# [Helper functions resolve_gene_name, get_gene_sequence, find_cas9_sites]
def resolve_gene_name(yeast_gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{yeast_gene}?content-type=application/json"
    response = requests.get(url)
    return response.json()[0]['id'] if response.status_code == 200 and len(response.json()) > 0 else None

def get_gene_sequence(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    response = requests.get(url)
    return response.text.strip() if response.status_code == 200 else None

def find_cas9_sites(gene_sequence, aa_pos):
    nuc_pos = (aa_pos - 1) * 3
    start, end = max(0, nuc_pos - 60), min(len(gene_sequence), nuc_pos + 60)
    region = gene_sequence[start:end]
    sites = []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'reverse'})
    return sorted(sites, key=lambda x: abs(x['pos'] - nuc_pos))

def get_mutated_template(original_seq, mut_idx, target_aa):
    codons = CODON_TABLE.get(target_aa, ["???"])
    seq_list = list(original_seq)
    seq_list[mut_idx:mut_idx+3] = list(codons[0])
    return "".join(seq_list)

# --- UI ---
st.title("🧬 Yeast AutoOligo: Proofreader Mode")

with st.sidebar:
    gene_input = st.text_input("Gene", "PHO13")
    residue = st.number_input("Residue #", value=123, min_value=1)
    mutation_aa = st.text_input("Mutation", "R").upper()
    run = st.button("Generate & Align", type="primary")

if run:
    gene_id = resolve_gene_name(gene_input)
    if gene_id:
        full_seq = get_gene_sequence(gene_id)
        mut_idx_gene = (residue - 1) * 3
        sites = find_cas9_sites(full_seq, residue)
        
        if sites:
            # We'll use the first (closest) site for the proofreader
            site = sites[0]
            start_idx = max(0, min(mut_idx_gene, site['pos']) - 15)
            # Align start_idx to codon boundary
            start_idx = (start_idx // 3) * 3
            end_idx = min(len(full_seq), max(mut_idx_gene + 3, site['pos'] + 23) + 15)
            
            wt_dna = full_seq[start_idx:end_idx]
            mut_dna = get_mutated_template(wt_dna, mut_idx_gene - start_idx, mutation_aa)
            
            # Identify PAM indices in this window
            pam_rel = site['pos'] - start_idx
            pam_indices = range(pam_rel + 20, pam_rel + 23) if site['strand'] == 'forward' else range(pam_rel, pam_rel + 3)
            mut_indices = range(mut_idx_gene - start_idx, mut_idx_gene - start_idx + 3)

            # Build HTML Table
            html = '<table class="align-table">'
            # AA Row
            html += '<tr><td class="label-cell">WT AA</td>'
            for i in range(0, len(wt_dna), 3): html += f'<td colspan="3" style="color:#888">{str(Seq(wt_dna[i:i+3]).translate())}</td>'
            html += '</tr>'
            
            # DNA Rows
            for label, dna in [("WT DNA", wt_dna), ("MUT DNA", mut_dna)]:
                html += f'<tr><td class="label-cell">{label}</td>'
                for i, b in enumerate(dna):
                    cls = ""
                    if i in mut_indices and i in pam_indices: cls = ' class="overlap-site"'
                    elif i in mut_indices: cls = ' class="mut-site"'
                    elif i in pam_indices: cls = ' class="pam-site"'
                    html += f'<td{cls}>{b}</td>'
                html += '</tr>'
            
            # Mut AA Row
            html += '<tr><td class="label-cell">MUT AA</td>'
            for i in range(0, len(mut_dna), 3): html += f'<td colspan="3" style="font-weight:bold;">{str(Seq(mut_dna[i:i+3]).translate())}</td>'
            html += '</tr></table>'
            
            st.markdown(html, unsafe_allow_html=True)
            st.markdown("🟢 **Green**: PAM Site | 🔴 **Red**: Mutation | 🟠 **Orange**: Overlap")
            
            st.subheader("📋 Ordering Sequences")
            st.code(f"Sense: {mut_dna}\nComplement: {str(Seq(mut_dna).complement())}", language="text")
