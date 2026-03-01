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

# --- CSS for Highlighting & Alignment ---
st.markdown("""
    <style>
    .align-table { font-family: 'Courier New', monospace; border-collapse: collapse; margin: 10px 0; background-color: #ffffff; width: 100%; }
    .align-table td { padding: 4px 1px; text-align: center; width: 20px; font-size: 1.1em; border: 1px solid #f0f0f0; }
    .label-cell { text-align: right !important; width: 120px !important; color: #666; font-size: 0.85em !important; padding-right: 15px !important; border-right: 3px solid #eee !important; font-weight: bold; }
    .pam-site { background-color: #d1ffbd; color: #006400; font-weight: bold; }
    .mut-site { background-color: #ffcccc; color: #cc0000; font-weight: bold; }
    .overlap-site { background-color: #ffe5b4; color: #8b4513; font-weight: bold; }
    .design-card { border: 2px solid #f0f2f6; padding: 20px; border-radius: 15px; margin-bottom: 25px; background-color: white; box-shadow: 2px 2px 10px rgba(0,0,0,0.05); }
    </style>
    """, unsafe_allow_html=True)

# --- Logic Functions ---
def resolve_gene_name(yeast_gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{yeast_gene}?content-type=application/json"
    r = requests.get(url)
    return r.json()[0]['id'] if r.status_code == 200 and r.json() else None

def get_gene_sequence(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    r = requests.get(url)
    return r.text.strip() if r.status_code == 200 else None

def find_cas9_sites(gene_sequence, aa_pos):
    nuc_pos = (aa_pos - 1) * 3
    window = 100
    start, end = max(0, nuc_pos - window), min(len(gene_sequence), nuc_pos + window)
    region = gene_sequence[start:end]
    sites = []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'reverse'})
    return sorted(sites, key=lambda x: abs(x['pos'] - nuc_pos))

# --- UI ---
st.title("🧬 Multi-Option CRISPR Designer")

with st.sidebar:
    st.header("Input Details")
    gene_input = st.text_input("Yeast Gene", "PHO13")
    residue = st.number_input("Residue #", value=123, min_value=1)
    mutation_aa = st.text_input("Mutation AA", "R").upper()
    run = st.button("Generate All Designs", type="primary", use_container_width=True)

if run:
    gene_id = resolve_gene_name(gene_input)
    if not gene_id:
        st.error("Gene not found.")
    else:
        full_seq = get_gene_sequence(gene_id)
        mut_idx_gene = (residue - 1) * 3
        sites = find_cas9_sites(full_seq, residue)
        
        if not sites:
            st.warning("No Cas9 sites found near this residue.")
        else:
            st.success(f"Found {len(sites)} potential sgRNA sites. Designs are ordered by proximity to residue {residue}.")
            
            for i, site in enumerate(sites[:5]): # Show top 5
                with st.container():
                    st.markdown(f'<div class="design-card">', unsafe_allow_html=True)
                    st.subheader(f"Option {i+1}: sgRNA on {site['strand'].capitalize()} Strand")
                    
                    # 1. Cloning Oligos
                    guide_20 = site['seq'][:-3] if site['strand'] == 'forward' else str(Seq(site['seq'][3:]).reverse_complement())
                    oligo_a = f"gatc{guide_20.lower()}gttttagagctag"
                    oligo_b = f"ctagctctaaaac{str(Seq(guide_20).reverse_complement()).lower()}"
                    
                    # 2. Repair Data
                    start_idx = (max(0, min(mut_idx_gene, site['pos']) - 12) // 3) * 3
                    end_idx = min(len(full_seq), max(mut_idx_gene + 3, site['pos'] + 23) + 12)
                    wt_dna = full_seq[start_idx:end_idx]
                    
                    mut_codon = CODON_TABLE.get(mutation_aa, ["???"])[0]
                    mut_dna_list = list(wt_dna)
                    rel_mut = mut_idx_gene - start_idx
                    mut_dna_list[rel_mut:rel_mut+3] = list(mut_codon)
                    mut_dna = "".join(mut_dna_list)
                    
                    # 3. Highlighting Logic
                    pam_rel = site['pos'] - start_idx
                    pam_idx = range(pam_rel+20, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+3)
                    mut_idx = range(rel_mut, rel_mut+3)

                    # 4. Display Logic
                    col1, col2 = st.columns([1, 1])
                    with col1:
                        st.markdown("**pML104 Cloning Oligos**")
                        st.code(f"A: {oligo_a}\nB: {oligo_b}", language="text")
                    with col2:
                        st.markdown("**Repair Template (80-90bp)**")
                        st.code(f"Sense: {mut_dna}\nComp : {str(Seq(mut_dna).complement())}", language="text")

                    # 5. The Proofreader Grid
                    html = '<table class="align-table">'
                    html += '<tr><td class="label-cell">WT PROTEIN</td>'
                    for j in range(0, len(wt_dna), 3): html += f'<td colspan="3" style="color:#888">{str(Seq(wt_dna[j:j+3]).translate())}</td>'
                    html += '</tr>'
                    
                    for label, dna in [("WT DNA", wt_dna), ("MUT DNA", mut_dna)]:
                        html += f'<tr><td class="label-cell">{label}</td>'
                        for idx, char in enumerate(dna):
                            cls = ""
                            if idx in mut_idx and idx in pam_idx: cls = ' class="overlap-site"'
                            elif idx in mut_idx: cls = ' class="mut-site"'
                            elif idx in pam_idx: cls = ' class="pam-site"'
                            html += f'<td{cls}>{char}</td>'
                        html += '</tr>'
                    
                    html += '<tr><td class="label-cell">MUT PROTEIN</td>'
                    for j in range(0, len(mut_dna), 3): html += f'<td colspan="3" style="font-weight:bold;">{str(Seq(mut_dna[j:j+3]).translate())}</td>'
                    html += '</tr></table>'
                    
                    st.markdown(html, unsafe_allow_html=True)
                    st.markdown('</div>', unsafe_allow_html=True)
