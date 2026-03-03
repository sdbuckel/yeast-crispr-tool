import streamlit as st
import requests
import re
from Bio.Seq import Seq

CODON_TABLE = {'A':['GCT','GCC','GCA','GCG'],'C':['TGT','TGC'],'D':['GAT','GAC'],'E':['GAA','GAG'],'F':['TTT','TTC'],'G':['GGT','GGC','GGA','GGG'],'H':['CAT','CAC'],'I':['ATT','ATC','ATA'],'K':['AAA','AAG'],'L':['TTA','TTG','CTT','CTC','CTA','CTG'],'M':['ATG'],'N':['AAT','AAC'],'P':['CCT','CCC','CCA','CCG'],'Q':['CAA','CAG'],'R':['CGT','CGC','CGA','CGG','AGA','AGG'],'S':['TCT','TCC','TCA','TCG','AGT','AGC'],'T':['ACT','ACC','ACA','ACG'],'V':['GTT','GTC','GTA','GTG'],'W':['TGG'],'Y':['TAT','TAC'],'*':['TAA','TAG','TGA']}

# --- CSS Styling ---
st.markdown("""
<style>
.citation-text{font-size:0.9em;color:#555;font-style:italic;margin-bottom:15px;}
.instruction-box{background-color:#f0f7ff;border-left:5px solid #007bff;padding:15px;border-radius:5px;margin-bottom:25px;color:#004085;}
.align-table{font-family:'Courier New',monospace;border-collapse:collapse;width:100%;margin:20px 0;}
.align-table td{padding:5px 2px;text-align:center;width:24px;border:1px solid #eee;}
.label-cell{text-align:right!important;width:140px!important;font-weight:bold;background-color:#f9f9f9;border-right:3px solid #ddd!important;}
.pam-site{background-color:#d1ffbd;color:#006400;font-weight:bold;}
.mut-site{background-color:#ffcccc;color:#cc0000;font-weight:bold;}
.silent-site{background-color:#fff9c4;color:#827717;font-weight:bold;}
.cut-mark{color:red;font-size:1.2em;font-weight:bold;}
.design-card{border:1px solid #ddd;padding:25px;border-radius:12px;margin-bottom:40px;background-color:white;}
.legend-box{font-size:0.85em;color:#444;margin:10px 0;padding:10px;background:#f9f9f9;border:1px solid #eee;}
</style>
""", unsafe_allow_html=True)

def resolve_gene_name(gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{gene}?content-type=application/json"
    r = requests.get(url)
    return r.json()[0]['id'] if r.status_code == 200 and r.json() else None

def get_gene_seq(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain;expand_5prime=100;expand_3prime=100"
    r = requests.get(url)
    return r.text.strip() if r.status_code == 200 else None

def find_sites(seq, m_pos):
    win = 100; s, e = max(0, m_pos-win), min(len(seq), m_pos+win)
    reg, res_list = seq[s:e], []
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', reg): res_list.append({'pos':s+m.start(),'seq':m.group(1),'strand':'forward'})
    for m in re.finditer(r'(?=(CC[ACGT][ACGT]{20}))', reg): res_list.append({'pos':s+m.start(),'seq':m.group(1),'strand':'reverse'})
    return sorted(res_list, key=lambda x: abs(x['pos']-m_pos))

# --- UI Header ---
st.title("🧬 Yeast CRISPR Oligo Designer")
st.markdown('<div class="citation-text">Laughery et al (2015) <i>Yeast</i> 32:711-720; Laughery & Wryck (2019) <i>CPMB</i> 129:e110.</div>', unsafe_allow_html=True)

st.markdown("""
<div class="instruction-box">
<strong>How to use:</strong><br>
1. Enter the systematic or common <strong>Gene Name</strong> (e.g., ADE2, PHO13).<br>
2. Enter the <strong>Residue #</strong> (Amino Acid position) you wish to mutate.<br>
3. Provide the 1-letter code for the <strong>New Amino Acid</strong>.<br>
4. Adjust the <strong>Slider</strong> to set the total length of the repair template (60-100bp).<br>
5. Review the visual alignment below to ensure the <strong>PAM site</strong> is broken and the <strong>Cut Site (▲)</strong> is near your mutation.
</div>
""", unsafe_allow_html=True)

with st.sidebar:
    st.header("Parameters")
    g_in = st.text_input("Gene Name", "PHO13").strip()
    res = st.number_input("Residue #", value=1, min_value=1)
    m_aa = st.text_input("New AA", "A").upper().strip()
    t_len = st.slider("Repair Template Size (bp)", 60, 100, 90)
    run = st.button("Generate Designs", type="primary")

if run:
    g_id = resolve_gene_name(g_in)
    if g_id:
        f_seq = get_gene_seq(g_id)
        off, c_end = 100, len(f_seq)-100
        m_idx = off + (res-1)*3
        if m_idx+3 > c_end: st.error("Residue out of bounds.")
        elif m_aa not in CODON_TABLE: st.error("Invalid AA code.")
        else:
            sites = find_sites(f_seq, m_idx)
            for i, site in enumerate(sites[:5]):
                with st.container():
                    st.markdown('<div class="design-card">', unsafe_allow_html=True)
                    st.subheader(f"Option {i+1}: {site['strand'].title()} Site")
                    half = t_len // 2
                    v_s = off + (((max(0, m_idx - half) - off)//3)*3)
                    v_e = v_s + t_len
                    wt_dna = f_seq[v_s:v_e]
                    r_mut = m_idx - v_s
                    m_dna_l = list(wt_dna); m_dna_l[r_mut:r_mut+3] = list(CODON_TABLE[m_aa][0])
                    p_rel = site['pos'] - v_s
                    crit = range(p_rel+21, p_rel+23) if site['strand']=='forward' else range(p_rel, p_rel+2)
                    m_dna_f, c_idx = "".join(m_dna_l), list(range(r_mut, r_mut+3))
                    
                    is_bk = False
                    if 0 <= min(crit) and max(crit) < len(wt_dna):
                        is_bk = any(m_dna_f[p] != wt_dna[p] for p in crit)
                    
                    if not is_bk:
                        scan = range(p_rel+18, p_rel+23) if site['strand']=='forward' else range(p_rel, p_rel+5)
                        for p in reversed(scan):
                            cs = ((v_s + p - off)//3)*3 + off - v_s
                            if 0 <= cs <= len(m_dna_f)-3:
                                codon = m_dna_f[cs:cs+3]; amino = str(Seq(codon).translate())
                                for syn in CODON_TABLE.get(amino, []):
                                    if syn != codon:
                                        test_l = list(m_dna_f); test_l[cs:cs+3] = list(syn); test_s = "".join(test_l)
                                        if (site['strand']=='forward' and test_s[p_rel+21:p_rel+23] != "GG") or (site['strand']=='reverse' and test_s[p_rel:p_rel+2] != "CC"):
                                            m_dna_f, is_bk = test_s, True
                                            c_idx.extend(range(cs, cs+3)); break
                            if is_bk: break
                    
                    mut_indices = [idx for idx, (m, w) in enumerate(zip(m_dna_f, wt_dna)) if m != w]
                    dis_dna_sense = "".join([c.lower() if idx in mut_indices else c.upper() for idx, c in enumerate(m_dna_f)])
                    dis_dna_comp = "".join([c.lower() if idx in mut_indices else c.upper() for idx, c in enumerate(str(Seq(m_dna_f).complement()))])

                    aa_wt = [str(Seq(wt_dna[j:j+3]).translate()) if (v_s+j>=off and v_s+j+3<=c_end) else "---" for j in range(0, len(wt_dna), 3)]
                    aa_mu = [str(Seq(m_dna_f[j:j+3]).translate()) if (v_s+j>=off and v_s+j+3<=c_end) else "---" for j in range(0, len(m_dna_f), 3)]
                    cut_idx = p_rel + 17
