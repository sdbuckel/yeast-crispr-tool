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
    .citation-text { font-size: 0.9em; color: #555; font-style: italic; margin-bottom: 15px; line-height: 1.4; }
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
    for j in range(0, len(dna_seq),
