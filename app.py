import streamlit as st
import requests
import re
from Bio.Seq import Seq

# 1. Global Codon Table
CODON_TABLE = {'A':['GCT','GCC','GCA','GCG'],'C':['TGT','TGC'],'D':['GAT','GAC'],'E':['GAA','GAG'],'F':['TTT','TTC'],'G':['GGT','GGC','GGA','GGG'],'H':['CAT','CAC'],'I':['ATT','ATC','ATA'],'K':['AAA','AAG'],'L':['TTA','TTG','CTT','CTC','CTA','CTG'],'M':['ATG'],'N':['AAT','AAC'],'P':['CCT','CCC','CCA','CCG'],'Q':['CAA','CAG'],'R':['CGT','CGC','CGA','CGG','AGA','AGG'],'S':['TCT','TCC','TCA','TCG','AGT','AGC'],'T':['ACT','ACC','ACA','ACG'],'V':['GTT','GTC','GTA','GTG'],'W':['TGG'],'Y':['TAT','TAC'],'*':['TAA','TAG','TGA']}

# --- CSS Styling ---
st.markdown("<style>.citation-text{font-size:0.9em;color:#555;font-style:italic;margin-bottom:15px;}.instruction-box{background-color:#f0f7ff;border-left:5px solid #007bff;padding:15px;border-radius:5px;margin-bottom:25px;color:#004085;}.align-table{font-family:'Courier New',monospace;border-collapse:collapse;width:100%;margin:20px 0;}.align-table td{padding:5px 2px;text-align:center;width:24px;border:1px solid #eee;}.label-cell{text-align:right!important;width:140px!important;font-weight:bold;background-color:#f9f9f9;border-right:3px solid #ddd!important;}.pam-site{background-color:#d1ffbd;color:#006400;font-weight:bold;}.mut-site{background-color:#ffcccc;color:#cc0000;font-weight:bold;}.silent-site{background-color:#fff9c4;color:#827717;font-weight:bold;}.design-card{border:1px solid #ddd;padding:25px;border-radius:12px;margin-bottom:40px;background-color:white;}</style>", unsafe_allow_html=True)

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

# --- UI Header ---
st.title("🧬 Yeast CRISPR Oligo Designer")

st.markdown('<div class="citation-text">This tool will assist you in making single point mutations in various yeast genes. The references include:<br>• Laughery et al (2015) <em>Yeast</em> volume 32 pages 711-720<br>• Laughery and Wryck (2019) <em>Current Protocols in Molecular Biology</em> Volume 129 pages e110.</div>', unsafe_allow_html=True)
st.markdown('<div class="instruction-box"><strong>How to use:</strong> Enter the systematic or common yeast gene name. Add the <strong>amino acid number</strong> to change in the <em>Residue #</em> box and the <strong>new amino acid</strong> using 1-letter abbreviations.</div>', unsafe_allow_html=True)

with st.sidebar:
    g_in = st.text_input("Gene Name", "PHO13").strip()
    res = st.number_input("Residue #", value=1, min_value=1)
    m_aa = st.text_input("New AA", "A").upper().strip()
    run = st.button("Generate Designs", type="primary")

if run:
    g_id = resolve_gene_name(g_in)
    if g_id:
        f_seq = get_gene_seq(g_id)
        off, c_end = 50, len(f_seq)-50
        m_idx = off + (res-1)*3
        if m_idx+3 > c_end: st.error("Residue out of bounds.")
        elif m_aa not in CODON_TABLE: st.error("Invalid AA code.")
        else:
            sites = find_sites(f_seq, m_idx)
            for i, site in enumerate(sites[:5]):
                with st.container():
                    st.markdown('<div class="design-card">', unsafe_allow_html=True)
                    st.subheader(f"Option {i+1}: {site['strand'].title()} Site")
                    v_s = off + (((max(0, min(m_idx, site['pos'])-12)) - off)//3)*3
                    v_e = min(len(f_seq), max(m_idx+3, site['pos']+23)+15)
                    wt_dna = f_seq[v_s:v_e]
                    r_mut = m_idx - v_s
                    m_dna_l = list(wt_dna); m_dna_l[r_mut:r_mut+3] = list(CODON_TABLE[m_aa][0])
                    p_rel = site['pos'] - v_s
                    crit = range(p_rel+21, p_rel+23) if site['strand']=='forward' else range(p_rel, p_rel+2)
                    m_dna_s1 = "".join(m_dna_l); is_bk = any(m_dna_s1[p] != wt_dna[p] for p in crit)
                    m_dna_f, c_idx = m_dna_s1, list(range(r_mut, r_mut+3))
                    if not is_bk:
                        for p in crit:
                            cs = (p//3)*3
                            if 0 <= cs <= len(m_dna_s1)-3:
                                oc = m_dna_s1[cs:cs+3]; aa_o = next((a for a, c in CODON_TABLE.items() if oc in c), None)
                                if aa_o:
                                    for syn in CODON_TABLE[aa_o]:
                                        if syn != oc and syn[p%3] != oc[p%3]:
                                            tmp = list(m_dna_s1); tmp[cs:cs+3] = list(syn); m_dna_f, is_bk = "".join(tmp), True
                                            c_idx.extend(range(cs, cs+3)); break
                            if is_bk: break
                    dis_dna = "".join([c.lower() if (idx in c_idx and c != wt_dna[idx]) else c.upper() for idx, c in enumerate(m_dna_f)])
                    aa_wt = [str(Seq(wt_dna[j:j+3]).translate()) if (v_s+j>=off and v_s+j+3<=c_end) else "---" for j in range(0, len(wt_dna), 3)]
                    aa_mu = [str(Seq(m_dna_f[j:j+3]).translate()) if (v_s+j>=off and v_s+j+3<=c_end) else "---" for j in range(0, len(m_dna_f), 3)]
                    
                    # Visual HTML Table
                    h = '<table class="align-table"><tr><td class="label-cell">WT PROT</td>'
                    for a in aa_wt: h += f'<td colspan="3" style="color:#777">{a}</td>'
                    h += '</tr><tr><td class="label-cell">WT DNA</td>'
                    fp = range(p_rel+20, p_rel+23) if site['strand']=='forward' else range(p_rel, p_rel+3)
                    for idx, char in enumerate(wt_dna): h += f'<td{" class=pam-site" if idx in fp else ""}>{char}</td>'
                    h += '</tr><tr><td class="label-cell">MUT DNA</td>'
                    for idx, char in enumerate(dis_dna):
                        cl = ' class="mut-site"' if idx in range(r_mut, r_mut+3) else (' class="silent-site"' if idx in c_idx else '')
                        h += f'<td{cl}>{char}</td>'
                    h += '</tr><tr><td class="label-cell">MUT PROT</td>'
                    for a in aa_mu: h += f'<td colspan="3" style="font-weight:bold;">{a}</td>'
                    h += '</tr></table>'
                    st.markdown(h, unsafe_allow_html=True)

                    # Export Block
                    g20 = site['seq'][:-3].upper() if site['strand']=='forward' else str(Seq(site['seq'][3:]).reverse_complement()).upper()
                    wt_p_l = "".join([a + "  " for a in aa_wt])[:len(wt_dna)]
                    mu_p_l = "".join([a + "  " for a in aa_mu])[:len(wt_dna)]
                    txt = f"Oligo A: GATC{g20}GTTTTAGAGCTAG\nOligo B: CTAGCTCTAAAAC{str(Seq(g20).reverse_complement()).upper()}\n\n"
                    txt += f"Repair (Sense): {dis_dna}\nRepair (Comp):  {str(Seq(m_dna_f).complement()).upper()}\n\n"
                    txt += f"WT PROT: {wt_p_l}\nWT DNA:  {wt_dna.upper()}\nMUT DNA: {dis_dna}\nMUT PROT:{mu_p_l}"
                    st.code(txt, language="text")
                    st.markdown('</div>', unsafe_allow_html=True)
    else: st.error("Gene not found.")
