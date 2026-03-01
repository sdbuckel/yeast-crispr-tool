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
    .warning-badge { background-color: #ff4b4b; color: white; padding: 5px 10px; border-radius: 5px; font-weight: bold; display: inline-block; margin-bottom: 10px; }
    .design-card { border: 1px solid #ddd; padding: 25px; border-radius: 12px; margin-bottom: 35px; background-color: white; box-shadow: 3px 3px 10px rgba(0,0,0,0.05); }
    </style>
    """, unsafe_allow_html=True)

# [Helper functions resolve_gene_name, get_gene_sequence, find_cas9_sites remain same]
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
    window = 60
    start, end = max(0, nuc_pos - window), min(len(gene_sequence), nuc_pos + window)
    region = gene_sequence[start:end]
    sites = []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'reverse'})
    return sorted(sites, key=lambda x: abs(x['pos'] - nuc_pos))

def attempt_strict_pam_disruption(sequence, pam_critical_indices):
    """Targets ONLY the GG or CC bases to break the PAM."""
    seq_list = list(sequence)
    for p_idx in pam_critical_indices:
        codon_start = (p_idx // 3) * 3
        orig_codon = "".join(seq_list[codon_start:codon_start+3])
        current_aa = next((aa for aa, codons in CODON_TABLE.items() if orig_codon in codons), None)
        
        if current_aa:
            for syn in CODON_TABLE[current_aa]:
                if syn != orig_codon:
                    pos_in_codon = p_idx % 3
                    if syn[pos_in_codon] != orig_codon[pos_in_codon]:
                        seq_list[codon_start:codon_start+3] = list(syn)
                        return "".join(seq_list), [p_idx], True
    return "".join(seq_list), [], False

# --- UI Layout ---
st.title("🧬 Yeast AutoOligo: Lab Edition")
with st.sidebar:
    gene_input = st.text_input("Gene", "PHO13").strip()
    residue = st.number_input("Residue #", value=123, min_value=1)
    mutation_aa = st.text_input("Mutation", "R").upper().strip()
    run = st.button("Generate Designs", type="primary", use_container_width=True)

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
                
                # 1. Oligo Design
                guide_20 = site['seq'][:-3] if site['strand'] == 'forward' else str(Seq(site['seq'][3:]).reverse_complement())
                oligo_a = f"gatc{guide_20.lower()}gttttagagctag"
                oligo_b = f"ctagctctaaaac{str(Seq(guide_20).reverse_complement()).lower()}"
                
                # 2. Base Sequences
                start_idx = (max(0, min(mut_idx_gene, site['pos']) - 12) // 3) * 3
                end_idx = min(len(full_seq), max(mut_idx_gene + 3, site['pos'] + 23) + 15)
                wt_dna = full_seq[start_idx:end_idx]
                
                # 3. Apply Primary Mutation
                rel_mut = mut_idx_gene - start_idx
                mut_dna_list = list(wt_dna)
                mut_dna_list[rel_mut:rel_mut+3] = list(CODON_TABLE.get(mutation_aa, ["???"])[0])
                mut_dna_base = "".join(mut_dna_list)
                
                # 4. Strict PAM Disruption (Only GG or CC)
                pam_rel = site['pos'] - start_idx
                # Critical bases are the last two for forward (NGG) or first two for reverse (CCN)
                critical_indices = range(pam_rel+21, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+2)
                mut_dna_final, silent_idx, success = attempt_strict_pam_disruption(mut_dna_base, critical_indices)
                
                if not success:
                    st.markdown('<div class="warning-badge">⚠️ PAM LOCKED: No silent mutation possible for the GG/CC motif</div>', unsafe_allow_html=True)

                # 5. Alignment Table
                aa_line_vals = [str(Seq(wt_dna[j:j+3]).translate()) for j in range(0, len(wt_dna), 3)]
                ma_line_vals = [str(Seq(mut_dna_final[j:j+3]).translate()) for j in range(0, len(mut_dna_final), 3)]
                
                html = '<table class="align-table"><tr><td class="label-cell">WT PROTEIN</td>'
                for aa in aa_line_vals: html += f'<td colspan="3" style="color:#888">{aa}</td>'
                html += '</tr><tr><td class="label-cell">WT DNA</td>'
                full_pam_range = range(pam_rel+20, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+3)
                for idx, char in enumerate(wt_dna):
                    cls = ' class="pam-site"' if idx in full_pam_range else ''
                    html += f'<td{cls}>{char}</td>'
                html += '</tr><tr><td class="label-cell">MUT DNA</td>'
                for idx, char in enumerate(mut_dna_final):
                    cls = ' class="mut-site"' if idx in range(rel_mut, rel_mut+3) else (' class="silent-site"' if idx in silent_idx else '')
                    html += f'<td{cls}>{char}</td>'
                html += '</tr><tr><td class="label-cell">MUT PROTEIN</td>'
                for aa in ma_line_vals: html += f'<td colspan="3" style="font-weight:bold;">{aa}</td>'
                html += '</tr></table>'
                st.markdown(html, unsafe_allow_html=True)

                # 6. Word Export
                st.markdown("### 📋 Copy for Microsoft Word")
                word_text = (
                    f"DESIGN OPTION {i+1} ({site['strand'].upper()})\n"
                    f"OLIGO A: {oligo_a}\n"
                    f"OLIGO B: {oligo_b}\n"
                    f"REPAIR SENSE: {mut_dna_final}\n"
                    f"REPAIR COMP : {str(Seq(mut_dna_final).complement())}\n\n"
                    f"ALIGNMENT (Use Courier New Font):\n"
                    f"WT AA:  " + "  ".join(aa_line_vals) + "\n"
                    f"WT DNA: {wt_dna}\n"
                    f"MT DNA: {mut_dna_final}\n"
                    f"MT AA:  " + "  ".join(ma_line_vals)
                )
                st.code(word_text, language="text")
                st.markdown('</div>', unsafe_allow_html=True)
