import streamlit as st
import requests
import re
from Bio.Seq import Seq

CODON_TABLE = {'A':['GCT','GCC','GCA','GCG'],'C':['TGT','TGC'],'D':['GAT','GAC'],'E':['GAA','GAG'],'F':['TTT','TTC'],'G':['GGT','GGC','GGA','GGG'],'H':['CAT','CAC'],'I':['ATT','ATC','ATA'],'K':['AAA','AAG'],'L':['TTA','TTG','CTT','CTC','CTA','CTG'],'M':['ATG'],'N':['AAT','AAC'],'P':['CCT','CCC','CCA','CCG'],'Q':['CAA','CAG'],'R':['CGT','CGC','CGA','CGG','AGA','AGG'],'S':['TCT','TCC','TCA','TCG','AGT','AGC'],'T':['ACT','ACC','ACA','ACG'],'V':['GTT','GTC','GTA','GTG'],'W':['TGG'],'Y':['TAT','TAC'],'*':['TAA','TAG','TGA']}

st.markdown("<style>.citation-text{font-size:0.9em;color:#555;font-style:italic;margin-bottom:15px;}.instruction-box{background-color:#f0f7ff;border-left:5px solid #007bff;padding:15px;border-radius:5px;margin-bottom:25px;color:#004085;}.align-table{font-family:'Courier New',monospace;border-collapse:collapse;width:100%;margin:20px 0;}.align-table td{padding:5px 2px;text-align:center;width:24px;border:1px solid #eee;}.label-cell{text-align:right!important;width:140px!important;font-weight:bold;background-color:#f9f9f9;border-right:3px solid #ddd!important;}.pam-site{background-color:#d1ffbd;color:#006400;font-weight:bold;}.mut-site{background-color:#ffcccc;color:#cc0000;font-weight:bold;}.silent-site{background-color:#fff9c4;color:#827717;font-weight:bold;}.cut-mark{color:red;font-size:1.2em;font-weight:bold;}.design-card{border:1px solid #ddd;padding:25px;border-radius:12px;margin-bottom:40px;background-color:white;}</style>", unsafe_allow_html=True)

def resolve_gene_name(gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{gene}?content-type=application/json"
    r = requests.get(url); return r.json()[0]['id'] if r.status_code == 200 and r.json() else None

def get_gene_seq(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain;expand_5prime=50;expand_3prime=50"
    r = requests.get(url); return r.text.strip() if r.status_code == 200 else None

def find_sites(seq, m_pos):
    win = 100; s, e = max(0, m_pos-win), min(len(seq), m_pos+win)
    reg, res_list = seq[s:e], []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', reg): res_list.append({'pos':s+m.start(),'seq':m.group(1),'strand':'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', reg): res_list.append({'pos':s+m.start(),'seq':m.group(1),'strand':'reverse'})
    return sorted(res_list, key=lambda x: abs(x['pos']-m_pos))

st.title("🧬 Yeast CRISPR Oligo Designer")
st.markdown('<div class="citation-text">This tool will assist you in making single point mutations in various yeast genes. The references include:<br>• Laughery et al (2015) <em>Yeast</em> volume 32 pages 711-720<br>• Laughery and Wryck (2019) <em>Current Protocols in Molecular Biology</em> Volume 129 pages e110.</div>', unsafe_allow_html=True)
st.markdown
