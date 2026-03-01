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

# --- CSS Styling (Clean Highlights) ---
st.markdown("""
    <style>
    .dna-seq { 
        font-family: 'Courier New', monospace; 
        font-size: 1.2em; 
        line-height: 1.8; 
        word-break: break-all; 
        background: #ffffff; 
        padding: 15px; 
        border: 1px solid #eee;
        border-radius: 8px; 
    }
    .pam-site { background-color: #d1ffbd; color: #006400; padding: 0px 2px; border-radius: 2px; }
    .mut-site { background-color: #ffcccc; color: #cc0000; padding: 0px 2px; border-radius: 2px; }
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
    start, end = max(0, nuc_pos - 60), min(len(gene_sequence), nuc_pos + 60)
    region = gene_sequence[start:end]
    sites = []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'reverse'})
    return sorted(sites, key=lambda x: abs(x['pos'] - nuc_pos))

# --- Main Interface ---
st.title("🧬 Yeast AutoOligo Designer")

with st.sidebar:
    st.header("Parameters")
    gene_input = st.text_input("Gene", "PHO13").strip()
    residue = st.number_input("Residue #", value=123, min_value=1)
    mutation_aa = st.text_input("Mutation", "R").upper().strip()
    st.divider()
    st.markdown("**Legend:**")
    st.markdown('<span style="background:#ffcccc; color:#cc0000; padding:2px 5px; border-radius:3px">Mutation Site</span>', unsafe_allow_html=True)
    st.markdown('<span style="background:#d1ffbd; color:#006400; padding:2px 5px; border-radius:3px">PAM Site</span>', unsafe_allow_html=True)
    run = st.button("Generate Designs", type="primary", use_container_width=True)

if run:
    gene_id = resolve_gene_name(gene_input)
    if not gene_id:
        st.error(f"Gene '{gene_input}' not found.")
    else:
        full_seq = get_gene_sequence(gene_id)
        sites = find_cas9_sites(full_seq, residue)
        
        if not sites:
            st.warning("No CRISPR sites found nearby.")
        else:
            for i, site in enumerate(sites[:3]):
                with st.expander(f"Design Option {i+1} ({site['strand'].capitalize()})"):
                    
                    # 1. Oligo Design
                    guide_20 = site['seq'][:-3] if site['strand'] == 'forward' else str(Seq(site['seq'][3:]).reverse_complement())
                    oligo_a = f"gatc{guide_20.lower()}gttttagagctag"
                    oligo_b = f"ctagctctaaaac{str(Seq(guide_20).reverse_complement()).lower()}"
                    
                    # 2. Repair Template (Sense and Complement)
                    mut_idx = (residue - 1) * 3
                    h_start, h_end = max(0, mut_idx - 40), min(len(full_seq), mut_idx + 43)
                    repair_sense = full_seq[h_start:h_end]
                    # Generate the complementary strand (not reverse complement, just complement)
                    repair_comp = str(Seq(repair_sense).complement())
                    
                    # Map positions for visual highlighting
                    pam_rel = site['pos'] - h_start
                    mut_rel = mut_idx - h_start
                    pam_range = range(pam_rel + 20, pam_rel + 23) if site['strand'] == 'forward' else range(pam_rel, pam_rel + 3)

                    # Build Visual HTML
                    html_dna = '<div class="dna-seq">'
                    for idx, char in enumerate(repair_sense):
                        if idx in range(mut_rel, mut_rel + 3):
                            html_dna += f'<span class="mut-site">{char}</span>'
                        elif idx in range(pam_range.start, pam_range.stop):
                            html_dna += f'<span class="pam-site">{char}</span>'
                        else:
                            html_dna += char
                    html_dna += '</div>'
                    
                    st.markdown("### 🔬 Visual Map (Sense Strand)")
                    st.markdown(html_dna, unsafe_allow_html=True)
                    
                    st.markdown("### 📋 Combined Sequences")
                    master_text = (
                        f"--- DESIGN {i+1} ---\n"
                        f"Oligo A: {oligo_a}\n"
                        f"Oligo B: {oligo_b}\n"
                        f"Repair Template (Sense): {repair_sense}\n"
                        f"Repair Template (Complement): {repair_comp}"
                    )
                    st.code(master_text, language="text")
