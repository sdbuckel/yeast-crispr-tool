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
    .design-card { border: 1px solid #ddd; padding: 25px; border-radius: 12px; margin-bottom: 35px; background-color: white; box-shadow: 3px 3px 10px rgba(0,0,0,0.05); }
    </style>
    """, unsafe_allow_html=True)

def resolve_gene_name(gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{gene}?content-type=application/json"
    r = requests.get(url)
    return r.json()[0]['id'] if r.status_code == 200 and r.json() else None

def get_gene_sequence_with_flanks(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain;expand_5prime=50;expand_3prime=50"
    r = requests.get(url)
    return r.text.strip() if r.status_code == 200 else None

def find_unique_cas9_sites(gene_sequence, mut_pos_in_seq):
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
st.title("🧬 Yeast CRISPR Lab Designer")
with st.sidebar:
    gene_input = st.text_input("Gene", "PHO13").strip()
    residue = st.number_input("Residue #", value=1, min_value=1)
    mutation_aa = st.text_input("Mutation", "A").upper().strip()
    run = st.button("Generate Designs", type="primary")

if run:
    gene_id = resolve_gene_name(gene_input)
    if gene_id:
        full_seq = get_gene_sequence_with_flanks(gene_id)
        offset = 50
        mut_idx = offset + (residue - 1) * 3
        cds_end = len(full_seq) - 50
        
        if mut_idx + 3 > cds_end:
            st.error(f"Residue {residue} is out of bounds.")
        elif mutation_aa not in CODON_TABLE:
            st.error(f"'{mutation_aa}' is invalid.")
        else:
            wt_aa = str(Seq(full_seq[mut_idx:mut_idx+3]).translate())
            if mutation_aa == wt_aa:
                st.error(f"Position {residue} is already {wt_aa}.")
            else:
                sites = find_unique_cas9_sites(full_seq, mut_idx)
                for i, site in enumerate(sites[:5]):
                    with st.container():
                        st.markdown('<div class="design-card">', unsafe_allow_html=True)
                        st.subheader(f"Option {i+1}: {'Forward' if site['strand']=='forward' else 'Reverse'} Site")

                        # Frame-aware window
                        raw_start = max(0, min(mut_idx, site['pos']) - 12)
                        v_start = offset + ((raw_start - offset) // 3) * 3
                        v_end = min(len(full_seq), max(mut_idx + 3, site['pos'] + 23) + 15)
                        wt_dna_seg = full_seq[v_start:v_end]
                        
                        rel_mut = mut_idx - v_start
                        mut_dna_list = list(wt_dna_seg)
                        mut_dna_list[rel_mut:rel_mut+3] = list(CODON_TABLE[mutation_aa][0])
                        mut_dna_step1 = "".join(mut_dna_list)
                        
                        pam_rel = site['pos'] - v_start
                        crit = range(pam_rel+21, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+2)
                        is_broken = any(mut_dna_step1[p] != wt_dna_seg[p] for p in crit)
                        
                        mut_dna_final = mut_dna_step1
                        changed_indices = list(range(rel_mut, rel_mut+3))
                        if not is_broken:
                            for p in crit:
                                c_start = (p // 3) * 3
                                if c_start < 0 or c_start + 3 > len(mut_dna_step1): continue
                                o_cod = mut_dna_step1[c_start:c_start+3]
                                aa_orig = next((a for a, c in CODON_TABLE.items() if o_cod in c), None)
                                if aa_orig:
                                    for syn in CODON_TABLE[aa_orig]:
                                        if syn != o_cod and syn[p % 3] != o_cod[p % 3]:
                                            tmp = list(mut_dna_step1)
                                            tmp[c_start:c_start+3] = list(syn)
                                            mut_dna_final = "".join(tmp)
                                            changed_indices.extend(range(c_start, c_start+3))
                                            is_broken = True; break
                                if is_broken: break

                        final_display_dna = ""
                        for idx, char in enumerate(mut_dna_final):
                            if idx in changed_indices and mut_dna_final[idx] != wt_dna_seg[idx]:
                                final_display_dna += char.lower()
                            else:
                                final_display_dna += char.upper()

                        # Formatting the text alignment manually (1 char per DNA base)
                        wt_aa_line = ""
                        mu_aa_line = ""
                        for j in range(0, len(wt_dna_seg), 3):
                            gp = v_start + j
                            in_cds = gp >= offset and gp + 3 <= cds_end
                            
                            # Place AA at first char of codon, then 2 spaces
                            if in_cds:
                                w_a = str(Seq(wt_dna_seg[j:j+3]).translate())
                                m_a = str(Seq(mut_dna_final[j:j+3]).translate())
                                wt_aa_line += w_a + "  "
                                mu_aa_line += m_a + "  "
                            else:
                                wt_aa_line += "---"
                                mu_aa_line += "---"

                        # Ensure DNA and AA lines are trimmed/padded equally
                        wt_aa_line = wt_aa_line[:len(wt_dna_seg)]
                        mu_aa_line = mu_aa_line[:len(wt
