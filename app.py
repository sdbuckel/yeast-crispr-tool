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

# --- CSS for Alignment View ---
st.markdown("""
    <style>
    .alignment-box { 
        font-family: 'Courier New', monospace; 
        white-space: pre; 
        background-color: #f8f9fa; 
        padding: 20px; 
        border-radius: 10px; 
        line-height: 1.2;
        overflow-x: auto;
        border: 1px solid #ddd;
    }
    .label { color: #888; font-size: 0.8em; }
    .pam-site { background-color: #d1ffbd; color: #006400; font-weight: bold; }
    .mut-site { background-color: #ffcccc; color: #cc0000; font-weight: bold; }
    </style>
    """, unsafe_allow_html=True)

def resolve_gene_name(yeast_gene):
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{yeast_gene}?content-type=application/json"
    response = requests.get(url)
    if response.status_code == 200 and len(response.json()) > 0: return response.json()[0]['id']
    return None

def get_gene_sequence(ensembl_id):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    response = requests.get(url)
    return response.text.strip() if response.status_code == 200 else None

def get_mutated_template(original_seq, mut_idx, target_aa):
    """Replaces the codon at mut_idx with the first available codon for target_aa."""
    codons = CODON_TABLE.get(target_aa, ["???"])
    new_codon = codons[0]
    seq_list = list(original_seq)
    seq_list[mut_idx:mut_idx+3] = list(new_codon)
    return "".join(seq_list)

# --- UI Setup ---
st.title("🧬 Yeast AutoOligo: Proofreader Mode")

with st.sidebar:
    st.header("Parameters")
    gene_input = st.text_input("Gene", "PHO13").strip()
    residue = st.number_input("Residue #", value=123, min_value=1)
    mutation_aa = st.text_input("Mutation", "R").upper().strip()
    run = st.button("Generate & Proofread", type="primary", use_container_width=True)

if run:
    gene_id = resolve_gene_name(gene_input)
    if not gene_id: st.error("Gene not found.")
    else:
        full_seq = get_gene_sequence(gene_id)
        # 1. Define region (30bp left, 30bp right of mutation for clarity)
        mut_idx_in_gene = (residue - 1) * 3
        start_idx = max(0, mut_idx_in_gene - 30)
        end_idx = min(len(full_seq), mut_idx_in_gene + 33)
        
        wt_dna = full_seq[start_idx:end_idx]
        mut_dna = get_mutated_template(wt_dna, mut_idx_in_gene - start_idx, mutation_aa)
        
        # 2. Translate (Handling frame)
        wt_aa = str(Seq(wt_dna).translate())
        mut_aa = str(Seq(mut_dna).translate())
        
        # 3. Highlight PAM (find sites in this specific window)
        # (Simplified: highlighting the central mutation codon)
        mut_rel = mut_idx_in_gene - start_idx
        
        # 4. Construct Stacked Visual
        def format_dna(dna, mut_pos):
            res = ""
            for i, base in enumerate(dna):
                if i in range(mut_pos, mut_pos+3): res += f'<span class="mut-site">{base}</span>'
                else: res += base
            return res

        def format_aa(aa):
            # Spread AA letters to align with 3-nucleotide codons
            return "  ".join(list(aa))

        st.markdown("### 🔍 Verification Alignment")
        st.markdown(f"""
        <div class="alignment-box">
<span class="label">WT AA :</span>  {format_aa(wt_aa)}
<span class="label">WT DNA:</span>  {format_dna(wt_dna, mut_rel)}
<span class="label">      </span>  {" " * mut_rel}|||
<span class="label">MUT DNA:</span>  {format_dna(mut_dna, mut_rel)}
<span class="label">MUT AA:</span>  {format_aa(mut_aa)}
        </div>
        """, unsafe_allow_html=True)

        st.info("💡 **Instructions:** Verify that only the highlighted codon changes the Amino Acid above/below it.")

        # Combined Copy Block
        st.subheader("📋 Final Sequences for Ordering")
        st.code(f"WT Fragment: {wt_dna}\nMutant Fragment: {mut_dna}\nComplement: {str(Seq(mut_dna).complement())}", language="text")
