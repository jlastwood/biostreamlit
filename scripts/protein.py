#Protein analysis using Biopython ProtParam
#- Computes statistics (MW, pI, instability, aromaticity, GRAVY)
#- Plots amino acid composition and hydropathy profile

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import streamlit as st
from collections import Counter
import requests

def proteinanal(sequence, seqid):
  # if the protein is blank we will use an ID
  # Example UniProt accession ID
  protein_id = seqid  # Human Hemoglobin subunit alpha
  if seqid != "ID":

# Query UniProt REST API
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
    response = requests.get(url)

    if response.status_code == 200:
      data = response.json()
    # Extract some fields
      protein_name = data["proteinDescription"]["recommendedName"]["fullName"]["value"]
      organism = data["organism"]["scientificName"]
      sequence = data["sequence"]["value"]

      st.write("Protein Name:", protein_name)
      st.write("Organism:", organism)
      st.write("Sequence Length:", len(sequence))
      st.write("First 50 residues:", sequence[:50])
    else:
      st.write("Error fetching data:", response.status_code)
# -----------------------------
# Input: protein sequence (one-letter amino acid codes)
# -----------------------------
#  sequence = "MKWVTFISLLFLFSSAYSLSGVRGSHHHHHHGATPQVSTPTLVEVSRNLGKVVSTQGNHQLENYCN"

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
  st.write("=== Protein summary ===")
  st.write(f"Length: {length} aa")
  st.write(f"Molecular weight: {mw/1000:.2f} kDa")
  st.write(f"Isoelectric point (pI): {pI:.2f}")
  st.write(f"Instability index: {instability:.2f}")
  st.write(f"Aromaticity: {aromaticity:.3f}")
  st.write(f"GRAVY (hydropathy): {gravy:.3f}")

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
  #plt.show()
  st.pyplot(fig1)

  hydrophobic = set("AILMFWV")
  hydrophilic = set("RNDQEKHSTYCGP")

  hydrophobic_count = sum(sequence.count(aa) for aa in hydrophobic)
  hydrophilic_count = sum(sequence.count(aa) for aa in hydrophilic)

  fig2, ax1 = plt.subplots()
  plt.pie([hydrophobic_count, hydrophilic_count],
        labels=["Hydrophobic", "Hydrophilic"],
        autopct="%1.1f%%",
        colors=["gold", "lightgreen"])
  plt.title("Hydrophobic vs Hydrophilic Composition")
#plt.show()
  st.pyplot(fig2)

  # this function can handle a list but we have only one
  sequences = [sequence]

  weights = [ProteinAnalysis(seq).molecular_weight() for seq in sequences]

  fig3, ax1 = plt.subplots()
  plt.hist(weights, bins=10, color="orchid", edgecolor="black")
  plt.xlabel("Molecular Weight (Da)")
  plt.ylabel("Count")
  plt.title("Protein Molecular Weight Distribution")
#plt.show()
  st.pyplot(fig3)

# Count amino acids
  counts = Counter(sequence)

# Plot
  fig4, ax1 = plt.subplots()
  plt.bar(counts.keys(), counts.values(), color="skyblue")
  plt.xlabel("Amino Acid")
  plt.ylabel("Frequency")
  plt.title("Amino Acid Frequency in Protein Sequence")
# plt.show()
  st.pyplot(fig4)

  hydrophobic = set("AILMFWV")
  hydrophilic = set("RNDQEKHSTYCGP")
