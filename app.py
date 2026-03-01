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
        offset = 50 # 5' UTR length
        mut_idx = offset + (residue - 1) * 3
        cds_end = len(full_seq) - 50
        
        if mut_idx + 3 > cds_end:
            st.error(f"Residue {residue} is out of bounds.")
        elif mutation_aa not in CODON_TABLE:
            st.error(f"'{mutation_aa}' is an invalid amino acid.")
        else:
            wt_aa = str(Seq(full_seq[mut_idx:mut_idx+3]).translate())
            if mutation_aa == wt_aa:
                st.error(f"**Entry Error:** Position {residue} is already {wt_aa}.")
            else:
                sites = find_unique_cas9_sites(full_seq, mut_idx)
                for i, site in enumerate(sites[:5]):
                    with st.container():
                        st.markdown('<div class="design-card">', unsafe_allow_html=True)
                        loc = " (5' UTR)" if site['pos'] < offset else (" (3' UTR)" if site['pos'] > cds_end else "")
                        st.subheader(f"Option {i+1}: PAM Site{loc}")

                        # --- CRITICAL FIX: FRAME-AWARE VISUALIZATION ---
                        raw_start = max(0, min(mut_idx, site['pos']) - 12)
                        # Ensure v_start relative to offset is always a multiple of 3
                        dist_from_start = raw_start - offset
                        v_start = offset + (dist_from_start // 3) * 3
                        
                        v_end = min(len(full_seq), max(mut_idx + 3, site['pos'] + 23) + 15)
                        wt_dna = full_seq[v_start:v_end]
                        
                        # Apply Primary Mutation
                        rel_mut = mut_idx - v_start
                        mut_dna_list = list(wt_dna)
                        mut_dna_list[rel_mut:rel_mut+3] = list(CODON_TABLE[mutation_aa][0])
                        mut_dna_step1 = "".join(mut_dna_list)
                        
                        # Silent PAM Disruption
                        pam_rel = site['pos'] - v_start
                        crit = range(pam_rel+21, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+2)
                        is_broken = any(mut_dna_step1[p] != wt_dna[p] for p in crit)
                        
                        mut_dna_final, sil_idx = mut_dna_step1, []
                        if is_broken:
                            st.markdown('<div class="success-badge">✨ PAM Disrupted</div>', unsafe_allow_html=True)
                        else:
                            for p in crit:
                                # Codon must align with v_start frame
                                c_start = (p // 3) * 3
                                if c_start < 0 or c_start + 3 > len(mut_dna_step1): continue
                                o_cod = mut_dna_step1[c_start:c_start+3]
                                aa = next((a for a, c in CODON_TABLE.items() if o_cod in c), None)
                                if aa:
                                    for syn in CODON_TABLE[aa]:
                                        if syn != o_cod and syn[p % 3] != o_cod[p % 3]:
                                            tmp = list(mut_dna_step1)
                                            tmp[c_start:c_start+3] = list(syn)
                                            mut_dna_final, sil_idx, is_broken = "".join(tmp), [p], True
                                            break
                                if is_broken: break

                        # Oligo Calculation
                        g20 = site['seq'][:-3].upper() if site['strand']=='forward' else str(Seq(site['seq'][3:]).reverse_complement()).upper()
                        ol_a = f"GATC{g20}GTTTTAGAGCTAG"
                        ol_b = f"CTAGCTCTAAAAC{str(Seq(g20).reverse_complement()).upper()}"
                        
                        # Proofreader Table with Frame Checking
                        aa_wt, aa_mu = [], []
                        for j in range(0, len(wt_dna), 3):
                            global_pos = v_start + j
                            if global_pos >= offset and global_pos + 3 <= cds_end:
                                aa_wt.append(str(Seq(wt_dna[j:j+3]).translate()))
                                aa_mu.append(str(Seq(mut_dna_final[j:j+3]).translate()))
                            else:
                                aa_wt.append("UTR")
                                aa_mu.append("UTR")
                        
                        html = '<table class="align-table"><tr><td class="label-cell">WT PROTEIN</td>'
                        for a in aa_wt: html += f'<td colspan="3">{a}</td>'
                        html += '</tr><tr><td class="label-cell">WT DNA</td>'
                        f_pam = range(pam_rel+20, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+3)
                        for idx, char in enumerate(wt_dna):
                            cls = ' class="pam-site"' if idx in f_pam else ''
                            html += f'<td{cls}>{char}</td>'
                        html += '</tr><tr><td class="label-cell">MUT DNA</td>'
                        for idx, char in enumerate(mut_dna_final):
                            cls = ' class="mut-site"' if idx in range(rel_mut, rel_mut+3) else (' class="silent-site"' if idx in sil_idx else '')
                            html += f'<td{cls}>{char}</td>'
                        html += '</tr><tr><td class="label-cell">MUT PROTEIN</td>'
                        for a in aa_mu: html += f'<td colspan="3" style="font-weight:bold;">{a}</td>'
                        html += '</tr></table>'
                        st.markdown(html, unsafe_allow_html=True)

                        # Word Export Block
                        word_text = f"OPTION {i+1} - {site['strand'].upper()} STRAND\n" + "="*40 + "\n"
                        word_text += f"Oligo A (Top):    {ol_a}\nOligo B (Bottom): {ol_b}\n\n"
                        word_text += f"Repair (Sense):   {mut_dna_final}\nRepair (Comp):    {str(Seq(mut_dna_final).complement())}\n\n"
                        word_text += f"WT DNA Segment:  {wt_dna}\nMUT DNA Segment: {mut_dna_final}"
                        st.code(word_text, language="text")
                        st.markdown('</div>', unsafe_allow_html=True)
