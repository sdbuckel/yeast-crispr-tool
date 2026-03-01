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
    .align-table { font-family: 'Courier New', monospace; border-collapse: collapse; margin: 10px 0; background-color: #ffffff; width: 100%; }
    .align-table td { padding: 4px 1px; text-align: center; width: 20px; font-size: 1.1em; border: 1px solid #f0f0f0; }
    .label-cell { text-align: right !important; width: 120px !important; color: #666; font-size: 0.85em !important; padding-right: 15px !important; border-right: 3px solid #eee !important; font-weight: bold; }
    .pam-site { background-color: #d1ffbd; color: #006400; font-weight: bold; }
    .mut-site { background-color: #ffcccc; color: #cc0000; font-weight: bold; }
    .silent-site { background-color: #fff9c4; color: #827717; font-weight: bold; } /* Yellow for Silent PAM Change */
    .design-card { border: 1px solid #ddd; padding: 25px; border-radius: 12px; margin-bottom: 30px; background-color: white; }
    </style>
    """, unsafe_allow_html=True)

# [Helper functions resolve_gene_name, get_gene_sequence, find_cas9_sites remain the same]
def resolve_gene_name(gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{gene}?content-type=application/json"
    r = requests.get(url)
    return r.json()[0]['id'] if r.status_code == 200 and r.json() else None

def get_gene_sequence(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    r = requests.get(url)
    return r.text.strip() if r.status_code == 200 else None

def find_cas9_sites(gene_sequence, aa_pos):
    nuc_pos = (aa_pos - 1) * 3
    window = 80
    start, end = max(0, nuc_pos - window), min(len(gene_sequence), nuc_pos + window)
    region = gene_sequence[start:end]
    sites = []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'reverse'})
    return sorted(sites, key=lambda x: abs(x['pos'] - nuc_pos))

def disrupt_pam_silently(sequence, pam_indices, start_offset):
    """Attempt to mutate a base in pam_indices without changing the AA sequence."""
    seq_list = list(sequence)
    mutated_indices = []
    
    # Check every codon in the sequence
    for i in range(0, len(sequence) - 2, 3):
        codon_indices = {i, i+1, i+2}
        # If this codon overlaps with the PAM
        overlap = codon_indices.intersection(pam_indices)
        if overlap:
            orig_codon = "".join(seq_list[i:i+3])
            current_aa = next((aa for aa, codons in CODON_TABLE.items() if orig_codon in codons), None)
            if not current_aa: continue
            
            # Try synonymous codons
            for syn_codon in CODON_TABLE[current_aa]:
                if syn_codon == orig_codon: continue
                # Does this syn_codon change one of the PAM bases?
                for idx in range(3):
                    if (i + idx) in pam_indices and syn_codon[idx] != orig_codon[idx]:
                        # Success! Apply mutation
                        seq_list[i:i+3] = list(syn_codon)
                        mutated_indices.append(i+idx)
                        return "".join(seq_list), mutated_indices
    return "".join(seq_list), mutated_indices

# --- UI ---
st.title("🧬 Yeast AutoOligo: Anti-Recutting Designer")

with st.sidebar:
    gene_input = st.text_input("Gene", "PHO13")
    residue = st.number_input("Residue #", value=123)
    mutation_aa = st.text_input("Mutation", "R").upper()
    run = st.button("Generate Designs", type="primary")

if run:
    gene_id = resolve_gene_name(gene_input)
    if gene_id:
        full_seq = get_gene_sequence(gene_id)
        mut_idx_gene = (residue - 1) * 3
        sites = find_cas9_sites(full_seq, residue)
        
        for i, site in enumerate(sites[:5]):
            with st.container():
                st.markdown('<div class="design-card">', unsafe_allow_html=True)
                st.subheader(f"Option {i+1} ({site['strand'].upper()} strand)")
                
                # 1. Base Fragment
                start_idx = (max(0, min(mut_idx_gene, site['pos']) - 12) // 3) * 3
                end_idx = min(len(full_seq), max(mut_idx_gene + 3, site['pos'] + 23) + 15)
                wt_dna = full_seq[start_idx:end_idx]
                
                # 2. Apply Desired Mutation
                mut_codon = CODON_TABLE.get(mutation_aa, ["???"])[0]
                rel_mut_idx = mut_idx_gene - start_idx
                mut_dna_list = list(wt_dna)
                mut_dna_list[rel_mut_idx:rel_mut_idx+3] = list(mut_codon)
                mut_dna_step1 = "".join(mut_dna_list)
                
                # 3. Apply Silent PAM Disruption
                pam_rel = site['pos'] - start_idx
                pam_indices = set(range(pam_rel+21, pam_rel+23)) if site['strand']=='forward' else set(range(pam_rel, pam_rel+2))
                # Note: We target the 'GG' of NGG or 'CC' of CCN
                
                mut_dna_final, silent_indices = disrupt_pam_silently(mut_dna_step1, pam_indices, start_idx)
                
                # 4. Display Logic
                # (Generating HTML and Word Block using mut_dna_final and highlighting silent_indices in yellow)
                st.markdown("### 🔍 Verified Alignment")
                # ... [HTML Table Logic similar to previous, adding silent-site class for silent_indices] ...
                
                st.info("🟡 **Yellow Highlight**: Indicates a silent mutation added to disrupt the PAM and prevent re-cutting.")
                st.code(f"Final Repair Template: {mut_dna_final}", language="text")
                st.markdown('</div>', unsafe_allow_html=True)
