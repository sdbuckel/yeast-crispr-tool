import streamlit as st
import requests
import re
from Bio.Seq import Seq

# 1. Global Codon Table
CODON_TABLE = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'C': ['TGT', 'TGC'], 'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'], 'F': ['TTT', 'TTC'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'], 'I': ['ATT', 'ATC', 'ATA'], 'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M': ['ATG'], 'N': ['AAT', 'AAC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'W': ['TGG'], 'Y': ['TAT', 'TAC'], '*': ['TAA', 'TAG', 'TGA']
}

# --- CSS Styling ---
st.markdown("""
    <style>
    .align-table { font-family: 'Courier New', monospace; border-collapse: collapse; margin: 10px 0; background-color: #ffffff; width: 100%; }
    .align-table td { padding: 4px 1px; text-align: center; width: 22px; font-size: 1.1em; border: 1px solid #f0f0f0; }
    .label-cell { text-align: right !important; width: 120px !important; color: #666; font-size: 0.85em !important; padding-right: 15px !important; border-right: 3px solid #eee !important; font-weight: bold; }
    .pam-site { background-color: #d1ffbd; color: #006400; font-weight: bold; }
    .mut-site { background-color: #ffcccc; color: #cc0000; font-weight: bold; }
    .silent-site { background-color: #fff9c4; color: #827717; font-weight: bold; }
    .success-badge { background-color: #28a745; color: white; padding: 5px 10px; border-radius: 5px; font-weight: bold; display: inline-block; margin-bottom: 10px; }
    .design-card { border: 1px solid #ddd; padding: 25px; border-radius: 12px; margin-bottom: 35px; background-color: white; box-shadow: 3px 3px 10px rgba(0,0,0,0.05); }
    </style>
    """, unsafe_allow_html=True)

def resolve_gene_name(gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{gene}?content-type=application/json"
    r = requests.get(url)
    return r.json()[0]['id'] if r.status_code == 200 and r.json() else None

def get_gene_sequence_with_flanks(ensembl_id):
    # Expanded to include 50bp of 5' and 3' flanking regions
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain;expand_5prime=50;expand_3prime=50"
    r = requests.get(url)
    return r.text.strip() if r.status_code == 200 else None

def find_unique_cas9_sites(gene_sequence, mut_pos_in_seq):
    # Search window around the mutation
    window = 100 
    start = max(0, mut_pos_in_seq - window)
    end = min(len(gene_sequence), mut_pos_in_seq + window)
    region = gene_sequence[start:end]
    sites = []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'reverse'})
    return sorted(sites, key=lambda x: abs(x['pos'] - mut_pos_in_seq))

# --- UI ---
st.title("🧬 Yeast CRISPR Designer (UTR-Aware)")
with st.sidebar:
    gene_input = st.text_input("Gene", "PHO13").strip()
    residue = st.number_input("Residue #", value=1, min_value=1)
    mutation_aa = st.text_input("Mutation", "A").upper().strip()
    run = st.button("Generate Designs", type="primary")

if run:
    gene_id = resolve_gene_name(gene_input)
    if gene_id:
        full_seq_with_flanks = get_gene_sequence_with_flanks(gene_id)
        # Because we added 50bp to the 5' end, the mutation position shifts
        offset = 50
        mut_idx_in_seq = offset + (residue - 1) * 3
        
        # CDS is from index 50 to (len - 50)
        cds_end = len(full_seq_with_flanks) - 50
        
        if mut_idx_in_seq + 3 > cds_end:
            st.error(f"Residue {residue} is out of bounds for the coding region of {gene_input}.")
        elif mutation_aa not in CODON_TABLE:
            st.error(f"'{mutation_aa}' is not a valid amino acid.")
        else:
            wt_codon = full_seq_with_flanks[mut_idx_in_seq:mut_idx_in_seq+3]
            wt_aa = str(Seq(wt_codon).translate())
            
            if mutation_aa == wt_aa:
                st.error(f"Position {residue} is already {wt_aa}. Please choose a different mutation.")
            else:
                all_sites = find_unique_cas9_sites(full_seq_with_flanks, mut_idx_in_seq)
                
                for i, site in enumerate(all_sites[:5]):
                    with st.container():
                        st.markdown('<div class="design-card">', unsafe_allow_html=True)
                        # Determine if PAM is in UTR
                        location_note = ""
                        if site['pos'] < offset: location_note = " (5' UTR)"
                        elif site['pos'] > cds_end: location_note = " (3' UTR)"
                        
                        st.subheader(f"Option {i+1}: PAM Site{location_note}")

                        # Window for visualization
                        v_start = (max(0, min(mut_idx_in_seq, site['pos']) - 12) // 3) * 3
                        v_end = min(len(full_seq_with_flanks), max(mut_idx_in_seq + 3, site['pos'] + 23) + 15)
                        wt_dna = full_seq_with_flanks[v_start:v_end]
                        
                        # Apply mutation
                        rel_mut = mut_idx_in_seq - v_start
                        mut_dna_list = list(wt_dna)
                        mut_dna_list[rel_mut:rel_mut+3] = list(CODON_TABLE[mutation_aa][0])
                        mut_dna_final = "".join(mut_dna_list)
                        
                        # Silent PAM disruption (only if PAM is in CDS)
                        # Note: If PAM is in UTR, we generally don't need silent mutation 
                        # because UTRs don't code for AAs, but for safety we check it anyway.
                        pam_rel = site['pos'] - v_start
                        crit_idx = range(pam_rel+21, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+2)
                        
                        # Logic for disruption markers
                        is_broken = any(mut_dna_final[p] != wt_dna[p] for p in crit_idx)
                        if is_broken:
                            st.markdown('<div class="success-badge">✨ PAM Disrupted by Mutation or UTR position</div>', unsafe_allow_html=True)
                        
                        # Oligos
                        if site['strand'] == 'forward':
                            guide_20 = site['seq'][:-3].upper()
                        else:
                            guide_20 = str(Seq(site['seq'][3:]).reverse_complement()).upper()
                        
                        ol_a = f"GATC{guide_20}GTTTTAGAGCTAG"
                        ol_b = f"CTAGCTCTAAAAC{str(Seq(guide_20).reverse_complement()).upper()}"
                        
                        # Formatting for proofreader
                        # (Translation only valid within CDS bounds)
                        st.code(f"Oligo A: {ol_a}\nOligo B: {ol_b}\nRepair: {mut_dna_final}", language="text")
                        st.markdown('</div>', unsafe_allow_html=True)
    else:
        st.error(f"Gene {gene_input} not found.")
