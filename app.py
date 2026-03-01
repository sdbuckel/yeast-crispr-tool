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

# 2. Logic Functions
def resolve_gene_name(yeast_gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{yeast_gene}?content-type=application/json"
    response = requests.get(url)
    if response.status_code == 200 and len(response.json()) > 0:
        return response.json()[0]['id']
    return None

def get_gene_sequence(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    response = requests.get(url)
    return response.text.strip() if response.status_code == 200 else None

def find_cas9_sites(gene_sequence, aa_pos):
    nuc_pos = (aa_pos - 1) * 3
    window = 105
    start = max(0, nuc_pos - window // 2)
    end = min(len(gene_sequence), nuc_pos + window // 2)
    region = gene_sequence[start:end]
    sites = []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'reverse'})
    return sorted(sites, key=lambda x: abs(x['pos'] - nuc_pos))

# --- UI Setup ---
st.set_page_config(page_title="Yeast CRISPR Designer", layout="wide")
st.title("🧬 Yeast AutoOligo CRISPR Helper")

with st.sidebar:
    st.header("Settings")
    gene_input = st.text_input("Gene Name", "PHO13")
    residue = st.number_input("Residue #", value=123, min_value=1)
    mutation = st.text_input("Mutation (One Letter)", "R").upper()
    run = st.button("Generate Design", type="primary")

if run:
    gene_id = resolve_gene_name(gene_input)
    if not gene_id:
        st.error("Gene not found.")
    else:
        dna_seq = get_gene_sequence(gene_id)
        sites = find_cas9_sites(dna_seq, residue)
        
        if not sites:
            st.warning("No sites found.")
        else:
            for i, site in enumerate(sites[:3]):
                with st.expander(f"Option {i+1}: {site['strand'].capitalize()} Strand (Pos {site['pos']})", expanded=(i==0)):
                    
                    # 1. Cloning Oligos
                    guide_20 = site['seq'][:-3] if site['strand'] == 'forward' else str(Seq(site['seq'][3:]).reverse_complement())
                    oligo_a = f"gatc{guide_20}gttttagagctag"
                    oligo_b = f"ctagctctaaaac{str(Seq(guide_20).reverse_complement())}"
                    
                    st.subheader("📋 pML104 Cloning Oligos")
                    st.caption("Click the icon on the right of the box to copy")
                    st.code(f"Oligo A: {oligo_a}\nOligo B: {oligo_b}", language="text")
                    
                    # 2. Repair Template Generation (Homology Arms)
                    mut_idx = (residue - 1) * 3
                    h_start = max(0, mut_idx - 40)
                    h_end = min(len(dna_seq), mut_idx + 43)
                    repair_template = dna_seq[h_start:h_end]
                    
                    # (In a real scenario, you'd apply the actual mutation string here)
                    st.subheader("🛠️ Repair Template (80bp typical)")
                    st.code(repair_template, language="text")
                    st.caption("Note: This template centers the mutation. Ensure PAM disruption is manually added.")
