#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import streamlit as st
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from itertools import islice

# Path to your FASTA file

def readfasta(fasta_gz): 

#fasta_gz = "files/uniprot.gz"
#fasta_file = "files/uniproto.fa"

# Open the gzip file in text mode
 with gzip.open(fasta_gz, "rt") as handle:
    # Parse all sequences
    records = list(SeqIO.parse(handle, "fasta"))

# Dictionary keyed by sequence IDs
 seq_dict = {}
# for record in SeqIO.parse(handle, "fasta"):
 for record in records:
    seq_dict[record.id] = str(record.seq)
    # st.write(record)

 #diroutput = dir(record)
 #st.write(diroutput)

#for index, record in enumerate(records, "fasta"):
#    print(
#        "index %i, ID = %s, length %i, with %i features"
#        % (index, record.id, len(record.seq), len(record.features))
#    )

# Concatenate all sequences into one string
 dna_seq = "".join(seq_dict.values())

# format of seqrecord object to create a file
#newrecord = SeqRecord(
#    Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
#    id="YP_025292.1",
#    name="HokC",
#    description="toxic membrane protein, small",
#)
#print(newrecord)

#print("=== Fasta record format ===")
 st.write(record.format("fasta"))

# Store sequences in a dictionary keyed by ID
#  seq_dict = {record.id: str(record.seq) for record in records}
# --- Outputs ---
# st.write("=== Individual sequences first 5 ===")
# st.write(f"Length: {len(seq_dict)}")
# for sid, seq in islice(seq_dict.items(),5):
#    st.write(f"{sid}: length {len(seq)}")

# print("\n=== Combined sequence print first 50 ===")
# print(f"Length: {len(dna_seq)}")
# print(dna_seq)
# (dna_seq[:50] + "...")  # print first 200 bases for preview

 return seq_dict
