Perfect — let’s refine the solution to **explicitly use Biopython’s ProtParam** for the protein statistics and graphs. ProtParam provides accurate calculations for molecular weight, isoelectric point, instability index, aromaticity, and GRAVY, so we’ll integrate it directly.

Here’s a clean Python script:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein analysis using Biopython ProtParam
- Computes statistics (MW, pI, instability, aromaticity, GRAVY)
- Plots amino acid composition and hydropathy profile
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# -----------------------------
# Input: protein sequence (one-letter amino acid codes)
# -----------------------------
sequence = "MKWVTFISLLFLFSSAYSLSGVRGSHHHHHHGATPQVSTPTLVEVSRNLGKVVSTQGNHQLENYCN"

# Initialize ProtParam
analysis = ProteinAnalysis(sequence)

# -----------------------------
# Compute statistics
# -----------------------------
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
kd_values = [analysis.protein_scale(window=1, param_dict={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,
                                                         'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                                                         'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,
                                                         'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3}, 
                                   window=9)[i] for i in range(len(sequence)-8)]

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

### What this script does
- Uses **ProtParam** for:
  - Molecular weight
  - Isoelectric point (pI)
  - Instability index
  - Aromaticity
  - GRAVY (hydropathy)
- Plots:
  - **Amino acid composition** (bar chart with counts and fractions)
  - **Kyte-Doolittle hydropathy profile** (line plot with sliding window)

---

👉 If you give me a **specific protein sequence or UniProt ID**, I can adapt this script to automatically fetch the sequence and run ProtParam analysis on it. Would you like me to extend this so it can pull sequences directly from UniProt and generate the graphs for any protein you choose?
