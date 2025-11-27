Got it — you’d like a **ready-to-run Python script** that fetches the protein sequence for **UniProt ID P16043** (which is *CD44 antigen, human*) and then runs **ProtParam analysis** with plots. Here’s a complete solution:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein analysis for UniProt ID P16043 (CD44 antigen, human)
- Fetch sequence from UniProt
- Use Biopython ProtParam for statistics
- Plot amino acid composition and hydropathy profile
"""

import requests
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# -----------------------------
# Fetch sequence from UniProt
# -----------------------------
uniprot_id = "P16043"
url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
response = requests.get(url)

if not response.ok:
    raise RuntimeError(f"Failed to fetch sequence for {uniprot_id}")

# Parse FASTA
lines = response.text.strip().split("\n")
sequence = "".join(line.strip() for line in lines if not line.startswith(">"))

print(f"Fetched sequence for {uniprot_id}, length {len(sequence)} aa")

# -----------------------------
# ProtParam analysis
# -----------------------------
analysis = ProteinAnalysis(sequence)

length = len(sequence)
mw = analysis.molecular_weight()
pI = analysis.isoelectric_point()
instability = analysis.instability_index()
aromaticity = analysis.aromaticity()
gravy = analysis.gravy()
aa_comp = analysis.count_amino_acids()
aa_percent = analysis.get_amino_acids_percent()

# -----------------------------
# Print summary
# -----------------------------
print("=== Protein summary ===")
print(f"Length: {length} aa")
print(f"Molecular weight: {mw/1000:.2f} kDa")
print(f"Isoelectric point (pI): {pI:.2f}")
print(f"Instability index: {instability:.2f}")
print(f"Aromaticity: {aromaticity:.3f}")
print(f"GRAVY (hydropathy): {gravy:.3f}")

# -----------------------------
# Amino acid composition plot
# -----------------------------
df_comp = pd.DataFrame({
    "AA": list(aa_comp.keys()),
    "Count": list(aa_comp.values()),
    "Fraction": [aa_percent[aa] for aa in aa_comp.keys()]
}).sort_values("AA")

sns.set(style="whitegrid", context="talk")
fig1, ax1 = plt.subplots(figsize=(12, 6))
sns.barplot(data=df_comp, x="AA", y="Fraction", color="#4C78A8", ax=ax1)
ax1.set_title("Amino acid composition (fraction)")
ax1.set_ylabel("Fraction")
ax1.set_xlabel("Residue")
for i, row in df_comp.iterrows():
    ax1.text(i, row["Fraction"] + 0.005, f'{row["Count"]}', ha="center", va="bottom", fontsize=9)
plt.tight_layout()
plt.show()

# -----------------------------
# Hydropathy profile (Kyte-Doolittle)
# -----------------------------
kd_values = analysis.protein_scale(
    param_dict={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,
                'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,
                'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3},
    window=9
)

fig2, ax2 = plt.subplots(figsize=(12, 6))
ax2.plot(range(1, len(kd_values)+1), kd_values, color="#F58518", lw=2)
ax2.axhline(0, color="black", lw=1, ls="--")
ax2.set_title("Kyte-Doolittle hydropathy (window=9)")
ax2.set_xlabel("Residue index")
ax2.set_ylabel("Hydropathy")
plt.tight_layout()
plt.show()
```

---

### How it works
- **Fetches sequence** directly from UniProt FASTA endpoint.
- Runs **ProtParam analysis** (Biopython) for MW, pI, instability, aromaticity, GRAVY.
- Plots:
  - **Amino acid composition** (bar chart with counts and fractions).
  - **Kyte-Doolittle hydropathy profile** (window size = 9).

---

Would you like me to extend this script to also **save the plots as PNG files** (so you can reuse them in reports), or keep it interactive with `plt.show()` only?
