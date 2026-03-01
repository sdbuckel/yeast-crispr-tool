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
    """Creates a protein string where AA sits on the 1st base of each codon."""
    aa_str = ""
    for j in range(0, len(dna_seq), 3):
        gp = v_start + j
        if gp >= offset and gp + 3 <= cds_end:
            codon = dna_seq[j:j+3]
            amino_acid = str(Seq(codon).translate())
            aa_str += amino_acid + "  " # AA + 2 spaces = 3 chars
        else:
            aa_str += "---"
    return aa_str[:len(dna_seq)]

# --- UI Setup ---
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
        offset, cds_end = 50, len(full_seq) - 50
        mut_idx = offset + (residue - 1) * 3
        
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
                        st.subheader(f"Option {i+1}: {'Forward' if site['strand']=='forward' else 'Reverse'} Strand Site")
                        
                        # Visualization bounds
                        v_start = offset + (((max(0, min(mut_idx, site['pos']) - 12)) - offset) // 3) * 3
                        v_end = min(len(full_seq), max(mut_idx + 3, site['pos'] + 23) + 15)
                        wt_dna_seg = full_seq[v_start:v_end]
                        
                        # Apply mutation
                        rel_mut = mut_idx - v_start
                        mut_dna_list = list(wt_dna_seg)
                        mut_dna_list[rel_mut:rel_mut+3] = list(CODON_TABLE[mutation_aa][0])
                        
                        # Silent PAM disruption
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

                        # Format for display (lowercase mutations)
                        display_mut_dna = "".join([c.lower() if (idx in changed_indices and c != wt_dna_seg[idx]) else c.upper() for idx, c in enumerate(mut_dna_final)])
                        
                        # Create aligned protein lines
                        wt_prot_line = format_alignment(wt_dna_seg, v_start, offset, cds_end)
                        mu_prot_line = format_alignment(mut_dna_final, v_start, offset, cds_end)
                        
                        # Calculate Oligos
                        g20 = site['seq'][:-3].upper() if site['strand']=='forward' else str(Seq(site['seq'][3:]).reverse_complement()).upper()
                        ol_a, ol_b = f"GATC{g20}GTTTTAGAGCTAG", f"CTAGCTCTAAAAC{str(Seq(g20).reverse_complement()).upper()}"
                        
                        # The "Word Block"
                        word_output = f"OPTION {i+1} LAB REPORT DATA\n" + "="*50 + "\n"
                        word_output += f"Oligo A: {ol_a}\nOligo B: {ol_b}\n\n"
                        word_output += f"Repair Template (Sense): {display_mut_dna}\n"
                        word_output += f"Repair Template (Comp):  {str(Seq(mut_dna_final).complement()).upper()}\n\n"
                        word_output += "PROOFREADING ALIGNMENT (Set font to Courier New in Word):\n"
                        word_output += f"WT PROTEIN: {wt_prot_line}\n"
                        word_output += f"WT DNA:     {wt_dna_seg.upper()}\n"
                        word_output += f"MUT DNA:    {display_mut_dna}\n"
                        word_output += f"MUT PROTEIN:{mu_prot_line}\n"
                        
                        st.code(word_output, language="text")
                        st.markdown("---")
    else:
        st.error(f"Gene {gene_input} not found.")
