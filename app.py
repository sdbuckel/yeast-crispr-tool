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
    .warning-badge { background-color: #ff4b4b; color: white; padding: 5px 10px; border-radius: 5px; font-weight: bold; display: inline-block; margin-bottom: 10px; }
    .design-card { border: 1px solid #ddd; padding: 25px; border-radius: 12px; margin-bottom: 35px; background-color: white; box-shadow: 3px 3px 10px rgba(0,0,0,0.05); }
    </style>
    """, unsafe_allow_html=True)

def resolve_gene_name(gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{gene}?content-type=application/json"
    r = requests.get(url)
    return r.json()[0]['id'] if r.status_code == 200 and r.json() else None

def get_gene_sequence(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    r = requests.get(url)
    return r.text.strip() if r.status_code == 200 else None

def find_unique_cas9_sites(gene_sequence, aa_pos):
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
st.title("🧬 Yeast CRISPR Multi-Site Designer")
with st.sidebar:
    gene_input = st.text_input("Gene", "PHO13").strip()
    residue = st.number_input("Residue #", value=134, min_value=1)
    mutation_aa = st.text_input("Mutation", "C").upper().strip()
    run = st.button("Generate Independent Designs", type="primary")

if run:
    gene_id = resolve_gene_name(gene_input)
    if gene_id:
        full_seq = get_gene_sequence(gene_id)
        mut_idx_gene = (residue - 1) * 3
        all_sites = find_unique_cas9_sites(full_seq, residue)
        
        for i, site in enumerate(all_sites[:5]):
            with st.container():
                st.markdown(f'<div class="design-card">', unsafe_allow_html=True)
                st.subheader(f"Option {i+1}: Independent PAM Site")

                start_idx = (max(0, min(mut_idx_gene, site['pos']) - 12) // 3) * 3
                end_idx = min(len(full_seq), max(mut_idx_gene + 3, site['pos'] + 23) + 15)
                wt_dna = full_seq[start_idx:end_idx]
                rel_mut = mut_idx_gene - start_idx
                
                # Apply Primary Mutation
                mut_dna_list = list(wt_dna)
                new_codon = CODON_TABLE.get(mutation_aa, ["???"])[0]
                mut_dna_list[rel_mut:rel_mut+3] = list(new_codon)
                mut_dna_step1 = "".join(mut_dna_list)
                
                # PAM Logic: Critical Indices (GG or CC)
                pam_rel = site['pos'] - start_idx
                if site['strand'] == 'forward':
                    crit_idx = range(pam_rel+21, pam_rel+23)
                else:
                    crit_idx = range(pam_rel, pam_rel+2)
                
                # Check for "Bonus" Self-Disruption
                is_broken = any(mut_dna_step1[p] != wt_dna[p] for p in crit_idx)
                
                mut_dna_final = mut_dna_step1
                silent_idx = []
                
                if is_broken:
                    st.markdown('<div class="success-badge">✨ Mutation Self-Disrupts PAM</div>', unsafe_allow_html=True)
                else:
                    # Attempt Silent disruption if the mutation didn't do it
                    for p_idx in crit_idx:
                        codon_start = (p_idx // 3) * 3
                        orig_codon = mut_dna_step1[codon_start:codon_start+3]
                        current_aa = next((aa for aa, codons in CODON_TABLE.items() if orig_codon in codons), None)
                        if current_aa:
                            for syn in CODON_TABLE.get(current_aa, []):
                                if syn != orig_codon and syn[p_idx % 3] != orig_codon[p_idx % 3]:
                                    temp_list = list(mut_dna_step1)
                                    temp_list[codon_start:codon_start+3] = list(syn)
                                    mut_dna_final = "".join(temp_list)
                                    silent_idx = [p_idx]
                                    is_broken = True
                                    break
                        if is_broken: break
                
                if not is_broken:
                    st.markdown('<div class="warning-badge">⚠️ PAM LOCKED</div>', unsafe_allow_html=True)

                # Pre-calculate Oligo and Complement strings to avoid SyntaxErrors
                if site['strand'] == 'forward':
                    guide_20 = site['seq'][:-3].upper()
                else:
                    guide_20 = str(Seq(site['seq'][3:]).reverse_complement()).upper()
                
                rev_guide_comp = str(Seq(guide_20).reverse_complement()).upper()
                oligo_a = f"GATC{guide_20}GTTTTAGAGCTAG"
                oligo_b = f"CTAGCTCTAAAAC{rev_guide_comp}"
                repair_comp = str(Seq(mut_dna_final).complement())
                
                # Table Data
                aa_line = [str(Seq(wt_dna[j:j+3]).translate()) for j in range(0, len(wt_dna), 3)]
                ma_line = [str(Seq(mut_dna_final[j:j+3]).translate()) for j in range(0, len(mut_dna_final), 3)]
                
                html = '<table class="align-table"><tr><td class="label-cell">WT PROTEIN</td>'
                for aa in aa_line: html += f'<td colspan="3" style="color:#888">{aa}</td>'
                html += '</tr><tr><td class="label-cell">WT DNA</td>'
                full_pam = range(pam_rel+20, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+3)
                for idx, char in enumerate(wt_dna):
                    cls = ' class="pam-site"' if idx in full_pam else ''
                    html += f'<td{cls}>{char}</td>'
                html += '</tr><tr><td class="label-cell">MUT DNA</td>'
                for idx, char in enumerate(mut_dna_final):
                    cls = ' class="mut-site"' if idx in range(rel_mut, rel_mut+3) else (' class="silent-site"' if idx in silent_idx else '')
                    html += f'<td{cls}>{char}</td>'
                html += '</tr><tr><td class="label-cell">MUT PROTEIN</td>'
                for aa in ma_line: html += f'<td colspan="3" style="font-weight:bold;">{aa}</td>'
                html += '</tr></table>'
                st.markdown(html, unsafe_allow_html=True)

                # Final Word Block
                word_text = f"OPTION {i+1}\nOligo A: {oligo_a}\nOligo B: {oligo_b}\n"
                word_text += f"Repair Sense: {mut_dna_final}\nRepair Comp:  {repair_comp}"
                st.code(word_text, language="text")
                st.markdown('</div>', unsafe_allow_html=True)
