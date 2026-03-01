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
    .overlap-site { background-color: #ffe5b4; color: #8b4513; font-weight: bold; }
    .design-card { border: 1px solid #ddd; padding: 25px; border-radius: 12px; margin-bottom: 30px; background-color: white; }
    </style>
    """, unsafe_allow_html=True)

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
    window = 80
    start, end = max(0, nuc_pos - window), min(len(gene_sequence), nuc_pos + window)
    region = gene_sequence[start:end]
    sites = []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', region):
        sites.append({'pos': start + m.start(), 'seq': m.group(1), 'strand': 'reverse'})
    return sorted(sites, key=lambda x: abs(x['pos'] - nuc_pos))

# --- UI ---
st.title("🧬 Yeast AutoOligo Designer")

with st.sidebar:
    st.header("Input")
    gene_input = st.text_input("Gene", "PHO13")
    residue = st.number_input("Residue #", value=123, min_value=1)
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
                
                # Oligos
                guide_20 = site['seq'][:-3] if site['strand'] == 'forward' else str(Seq(site['seq'][3:]).reverse_complement())
                oligo_a = f"gatc{guide_20.lower()}gttttagagctag"
                oligo_b = f"ctagctctaaaac{str(Seq(guide_20).reverse_complement()).lower()}"
                
                # Sequence Logic
                start_idx = (max(0, min(mut_idx_gene, site['pos']) - 12) // 3) * 3
                end_idx = min(len(full_seq), max(mut_idx_gene + 3, site['pos'] + 23) + 15)
                wt_dna = full_seq[start_idx:end_idx]
                mut_codon = CODON_TABLE.get(mutation_aa, ["???"])[0]
                mut_dna_list = list(wt_dna)
                rel_mut = mut_idx_gene - start_idx
                mut_dna_list[rel_mut:rel_mut+3] = list(mut_codon)
                mut_dna = "".join(mut_dna_list)
                
                # HTML Proofreader (UI only)
                pam_rel = site['pos'] - start_idx
                pam_idx = range(pam_rel+20, pam_rel+23) if site['strand']=='forward' else range(pam_rel, pam_rel+3)
                mut_idx = range(rel_mut, rel_mut+3)

                html = '<table class="align-table">'
                html += '<tr><td class="label-cell">WT PROTEIN</td>'
                for j in range(0, len(wt_dna), 3): html += f'<td colspan="3">{str(Seq(wt_dna[j:j+3]).translate())}</td>'
                html += '</tr><tr><td class="label-cell">WT DNA</td>'
                for idx, char in enumerate(wt_dna):
                    cls = ' class="pam-site"' if idx in pam_idx else ''
                    html += f'<td{cls}>{char}</td>'
                html += '</tr><tr><td class="label-cell">MUT DNA</td>'
                for idx, char in enumerate(mut_dna):
                    cls = ' class="mut-site"' if idx in mut_idx else ''
                    html += f'<td{cls}>{char}</td>'
                html += '</tr><tr><td class="label-cell">MUT PROTEIN</td>'
                for j in range(0, len(mut_dna), 3): html += f'<td colspan="3" style="font-weight:bold;">{str(Seq(mut_dna[j:j+3]).translate())}</td>'
                html += '</tr></table>'
                st.markdown(html, unsafe_allow_html=True)

                # --- WORD FRIENDLY COPY BLOCK ---
                st.markdown("### 📝 Copy to Word")
                
                # Generate a text-based alignment
                aa_line = "WT AA:  " + "  ".join([str(Seq(wt_dna[j:j+3]).translate()) for j in range(0, len(wt_dna), 3)])
                wt_line = "WT DNA: " + wt_dna
                mt_line = "MUT DNA:" + mut_dna
                ma_line = "MUT AA: " + "  ".join([str(Seq(mut_dna[j:j+3]).translate()) for j in range(0, len(mut_dna), 3)])
                
                word_export = (
                    f"DESIGN OPTION {i+1}\n"
                    f"Gene: {gene_input} | Residue: {residue} | Mutation: {mutation_aa}\n"
                    f"{'='*40}\n"
                    f"CLONING OLIGOS (pML104)\n"
                    f"Oligo A: {oligo_a}\n"
                    f"Oligo B: {oligo_b}\n\n"
                    f"REPAIR TEMPLATE\n"
                    f"Sense:      {mut_dna}\n"
                    f"Complement: {str(Seq(mut_dna).complement())}\n\n"
                    f"PROOFREADING ALIGNMENT\n"
                    f"{aa_line}\n"
                    f"{wt_line}\n"
                    f"{mt_line}\n"
                    f"{ma_line}\n"
                    f"{'-'*40}\n"
                )
                
                st.code(word_export, language="text")
                st.caption("Copy the box above. In Word, use a monospace font like 'Courier New'.")
                st.markdown('</div>', unsafe_allow_html=True)
