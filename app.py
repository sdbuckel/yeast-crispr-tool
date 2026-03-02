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
    .instruction-box { background-color: #f0f7ff; border-left: 5px solid #007bff; padding: 15px; border-radius: 5px; margin-bottom: 25px; color: #004085; }
    .align-table { font-family: 'Courier New', monospace; border-collapse: collapse; margin: 20px 0; background-color: #ffffff; width: 100%; }
    .align-table td { padding: 5px 2px; text-align: center; width: 24px; font-size: 1.1em; border: 1px solid #eee; }
    .label-cell { text-align: right !important; width: 140px !important; color: #444; font-size: 0.9em !important; padding-right: 15px !important; border-right: 3px solid #ddd !important; font-weight: bold; background-color: #f9f9f9; }
    .pam-site { background-color: #d1ffbd; color: #006400; font-weight: bold; }
    .mut-site { background-color: #ffcccc; color: #cc0000; font-weight: bold; }
    .silent-site { background-color: #fff9c4; color: #827717; font-weight: bold; }
    .design-card { border: 1px solid #ddd; padding: 25px; border-radius: 12px; margin-bottom: 40px; background-color: white; box-shadow: 2px 2px 8px rgba(0,0,0,0.05); }
    </style>
    """, unsafe_allow_html=True)

# [Helper functions: resolve_gene_name, get_gene_sequence_with_flanks, find_unique_cas9_sites, format_alignment remain same]
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

def format_alignment(dna_seq, v_start, offset, cds_end):
    aa_str = ""
    for j in range(0, len(dna_seq), 3):
        gp = v_start + j
        if gp >= offset and gp + 3 <= cds_end:
            aa_str += str(Seq(dna_seq[j:j+3]).translate()) + "  "
        else:
            aa_str += "---"
    return aa_str[:len(dna_seq)]

# --- UI Header ---
st.title("🧬 Yeast CRISPR Lab Designer")

st.markdown("""
    <div class="instruction-box">
    <strong>How to use:</strong> Enter the systematic or common yeast gene name (e.g., PHO13). 
    Add the <strong>Amino Acid number</strong> you want to change in the <em>Residue #</em> box 
    and add the <strong>Amino Acid</strong> you want to change it to using 
    <strong>single-letter abbreviations</strong> (e.g., 'C' for Cysteine).
    </div>
    """, unsafe_allow_html=True)



with st.sidebar:
    st.header("Design Parameters")
    gene_input = st.text_input("Gene Name", "PHO13").strip()
    residue = st.number_input("Residue # (Position)", value=1, min_value=1)
    mutation_aa = st.text_input("New Amino Acid (1-letter)", "A").upper().strip()
    
    with st.expander("Amino Acid Cheat Sheet"):
        st.caption("A: Ala | C: Cys | D: Asp | E: Glu | F: Phe")
        st.caption("G: Gly | H: His | I: Ile | K: Lys | L: Leu")
        st.caption("M: Met | N: Asn | P: Pro | Q: Gln | R: Arg")
        st.caption("S: Ser | T: Thr | V: Val | W: Trp | Y: Tyr")
    
    run = st.button("Generate Designs", type="primary")

if run:
    gene_id = resolve_gene_name(gene_input)
    if gene_id:
        full_seq = get_gene_sequence_with_flanks(gene_id)
        offset, cds_end = 50, len(full_seq) - 50
        mut_idx = offset + (residue - 1) * 3
        
        if mut_idx + 3 > cds_end:
            st.error(f"Residue {residue} is out of bounds for {gene_input}.")
        elif mutation_aa not in CODON_TABLE:
            st.error(f"'{mutation_aa}' is not a valid 1-letter amino acid code.")
        else:
            wt_aa = str(Seq(full_seq[mut_idx:mut_idx+3]).translate())
            if mutation_aa == wt_aa:
                st.error(f"Residue {residue} is already {wt_aa}. Please choose a different target mutation.")
            else:
                sites = find_unique_cas9_sites(full_seq, mut_idx)
                if not sites:
                    st.warning("No PAM sites found nearby.")
                
                for i, site in enumerate(sites[:5]):
                    with st.container():
                        st.markdown('<div class="design-card">', unsafe_allow_html=True)
                        st.subheader(f"Option {i+1}: {'Forward' if site['strand']=='forward' else 'Reverse'} Strand Site")
                        
                        v_start = offset + (((max(0, min(mut_idx, site['pos']) - 12)) - offset) // 3) * 3
                        v_end = min(len(full_seq), max(mut_idx + 3, site['pos'] + 23) + 15)
                        wt_dna_seg = full_seq[v_start:v_end]
                        
                        rel_mut = mut_idx - v_start
                        mut_dna_list = list(wt_dna_seg)
                        mut_dna_list[rel_mut:rel_mut+3] = list(CODON_TABLE[mutation_aa][0])
                        
                        pam_rel = site['pos'] - v_start
                        crit = range(pam_rel+21, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+2)
                        mut_dna_step1 = "".join(mut_dna_list)
                        is_broken = any(mut_dna_step1[p] != wt_dna_seg[p] for p in crit)
                        
                        mut_dna_final, changed_indices = mut_dna_step1, list(range(rel_mut, rel_mut+3))
                        if not is_broken:
                            for p in crit:
                                c_start = (p // 3) * 3
                                if 0 <= c_start <= len(mut_dna_step1)-3:
                                    o_cod = mut_dna_step1[c_start:c_start+3]
                                    aa_orig = next((a for a, c in CODON_TABLE.items() if o_cod in c), None)
                                    if aa_orig:
                                        for syn in CODON_TABLE[aa_orig]:
                                            if syn != o_cod and syn[p % 3] != o_cod[p % 3]:
                                                tmp = list(mut_dna_step1)
                                                tmp[c_start:c_start+3] = list(syn)
                                                mut_dna_final, is_broken = "".join(tmp), True
                                                changed_indices.extend(range(c_start, c_start+3))
                                                break
                                if is_broken: break

                        display_mut_dna = "".join([c.lower() if (idx in changed_indices and c != wt_dna_seg[idx]) else c.upper() for idx, c in enumerate(mut_dna_final)])
                        
                        # Visual Table
                        aa_wt_row = [str(Seq(wt_dna_seg[j:j+3]).translate()) if (v_start+j >= offset and v_start+j+3 <= cds_end) else "---" for j in range(0, len(wt_dna_seg), 3)]
                        aa_mu_row = [str(Seq(mut_dna_final[j:j+3]).translate()) if (v_start+j >= offset and v_start+j+3 <= cds_end) else "---" for j in range(0, len(mut_dna_final), 3)]
                        
                        html = '<table class="align-table"><tr><td class="label-cell">WT PROTEIN</td>'
                        for a in aa_wt_row: html += f'<td colspan="3" style="color:#777">{a}</td>'
                        html += '</tr><tr><td class="label-cell">WT DNA</td>'
                        f_pam = range(pam_rel+20, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+3)
                        for idx, char in enumerate(wt_dna_seg):
                            cls = ' class="pam-site"' if idx in f_pam else ''
                            html += f'<td{cls}>{char}</td>'
                        html += '</tr><tr><td class="label-cell">MUT DNA</td>'
                        for idx, char in enumerate(display_mut_dna):
                            cls = ' class="mut-site"' if idx in range(rel_mut, rel_mut+3) else (' class="silent-site"' if idx in changed_indices else '')
                            html += f'<td{cls}>{char}</td>'
                        html += '</tr><tr><td class="label-cell">MUT PROTEIN</td>'
                        for a in aa_mu_row: html += f'<td colspan="3" style="font-weight:bold;">{a}</td>'
                        html += '</tr></table>'
                        st.markdown(html, unsafe_allow_html=True)

                        # Word Block
                        wt_prot_line = format_alignment(wt_dna_seg, v_start, offset, cds_end)
                        mu_prot_line = format_alignment(mut_dna_final, v_start, offset, cds_end)
                        g20 = site['seq'][:-3].upper() if site['strand']=='forward' else str(Seq(site['seq'][3:]).reverse_complement()).upper()
                        ol_a, ol_b = f"GATC{g20}GTTTTAGAGCTAG", f"CTAGCTCTAAAAC{str(Seq(g20).reverse_complement()).upper()}"
                        
                        word_output = f"--- OPTION {i+1} LAB DATA ---\n"
                        word_output += f"Oligo A: {ol_a}\nOligo B: {ol_b}\n\n"
                        word_output += f"Repair (Sense): {display_mut_dna}\n"
                        word_output += f"Repair (Comp):  {str(Seq(mut_dna_final).complement()).upper()}\n\n"
                        word_output += "ALIGNMENT (Set font to Courier New):\n"
                        word_output += f"WT PROT: {wt_prot_line}\n"
                        word_output += f"WT DNA:  {wt_dna_seg.upper()}\n"
                        word_output += f"MUT DNA: {display_mut_dna}\n"
                        word_output += f"MUT PROT:{mu_prot_line}\n"
                        st.code(word_output, language="text")
                        
                        st.markdown('</div>', unsafe_allow_html=True)
    else:
        st.error(f"Gene '{gene_input}' not found in Ensembl database.")
