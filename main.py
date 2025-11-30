from PIL import Image
import numpy as np
import io
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import neatbio as nt
import neatbio.sequtils as utils
from collections import Counter
from Bio import SeqIO
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from fastaframes import to_df
from scripts.expression import expression
from scripts.readfasta import readfasta
from itertools import islice
from scripts.protein import proteinanal

def seqvalid(seq):

  #seq = 'actgtactgzzoootcga'

  seq = set(seq.lower())
  IUPAC = ['a','c','g','t','r','y','s','w','k','m','b','d','h','v','n','.','-','i','f','q','l','p','e']

  for char in IUPAC:
   if char in seq:
      seq.remove(char)
  if len(seq) > 0:
   st.write( 'sequence contains non IUPAC characters')
   st.write(seq)
  else:
   st.write( 'sequence passes check')

def replacer(m):
   return f"{m.group(1)}**{m.group(2)}**{m.group(3)}"

def delta(x, y):
    return 0 if x == y else 1

def M(seq1, seq2, i, j, k):
    return sum(delta(x, y) for x, y in zip(seq1[i:i + k], seq2[j:j + k]))


def makeMatrix(seq1, seq2, k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1, seq2, i, j, k) for j in range(m - k + 1)]
            for i in range(n - k + 1)]


def plotMatrix(M, t, seq1, seq2, nonblank=chr(0x25A0), blank=' '):
    print(' |' + seq2)
    print('-' * (2 + len(seq2)))
    for label, row in zip(seq1, M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)


def dotplot(seq1, seq2, k=1, t=1):
    M = makeMatrix(seq1, seq2, k)
    plotMatrix(M, t, seq1, seq2)  # experiment with character choice


# Convert to Fxn
def dotplotx(seq1, seq2):
    plt.imshow(np.array(makeMatrix(seq1, seq2, 1)))
    # on x-axis list all sequences of seq 2
    xt = plt.xticks(np.arange(len(list(seq2))), list(seq2))
    # on y-axis list all sequences of seq 1
    yt = plt.yticks(np.arange(len(list(seq1))), list(seq1))
    plt.show()


def gc_content(seq):
    result = float(str(seq).count('G') + str(seq).count('C')) / len(seq) * 100
    return result


def at_content(seq):
    result = float(str(seq).count('A') + str(seq).count('T')) / len(seq) * 100
    return result


def main():
    """A Simple Streamlit App """
    st.title("BioInformatics App")
    if 'seq_length' not in st.session_state:
      st.session_state['seq_length'] = 0
    if 'seq_dict' not in st.session_state:
      st.session_state['seq_dict'] = {}
    activity = ['Intro', 'Fasta', 'Expression', 'Protein', 'DotPlot', "About"]
    choice = st.sidebar.selectbox("Select Activity", activity)
    if choice == 'Intro':
        st.subheader("Intro")
        st.write("An app using the Bioinformatics library to read Fasta files and search for patterns")
    elif choice == 'Fasta':
        st.subheader("Upload a Fasta gz File")
        st.write("Upload a fasta file and read into a session variable for processing.  Remember this is a gz file")
        fasta = st.file_uploader("Upload FASTA File", type=[".gz"])
        if fasta is not None:
            st.session_state.seq_dict = readfasta(fasta)
            st.session_state.seq_length = len(st.session_state.seq_dict)
            st.write(f"Length: {len(st.session_state.seq_dict)}")
            st.write("The first five sequences with their length")
            for sid, seq in islice(st.session_state.seq_dict.items(),5):
               st.write(f"{sid}: length {len(seq)}")

            #st.write(f"Length: {len(dna_seq)}")
            #st.write(dna_seq[:50] + "...")  # print first 200 bases for preview

    elif choice == 'Expression':
        st.subheader("Expression Search")
        st.write("Enter the expression to search.  You can use a simple format like ATCGNWYS.  You can also use [] in the expression, such as [A|C]TC[AC]WYS.  After entering the search press enter, the results of the expression are displayed. Spances and special characters are ignored. ") 
        exp_text = st.text_input("enter Expression to search")
        #  do the work of creating the regex expressions
        (regex, explength)  = expression(exp_text)
        st.write("This is the regular expression for search", ", ".join(regex), "has ", len(regex), "iterations and is a length of ", explength)
        st.write("The number of sequences is ", st.session_state.seq_length, " If the sequences are 0, go to Upload Fasta.  To perform a search, click Analyse Sequence")
        distance = st.slider("select the distance", max_value = 1000, value=100)
        stopseq = st.slider("limit sequence search", max_value = 1000000, value=1000000)
        if st.session_state.seq_length > 0:
         gobutton = st.button ("Analyse Sequence")
         matchcount = 0
         matchmore = 0
         if gobutton:
          for idx, regex_part in enumerate(regex):
           st.write("search", regex_part)
           posfirst = idx
           for sid, seq in islice(st.session_state.seq_dict.items(),stopseq): 
             match = re.search(regex_part, seq)
             if match:
              start = match.start() +1
              end = match.end()
              group = match.group()
              # distance between groups
              distanceseq = end - start - explength +1
              endfirst = start + posfirst 
              startsend = end - explength +2  
              if  distanceseq > distance:
                matchmore=matchmore+1
              else:
                st.write("Found:", sid, " start ", match.start(), " span ", match.span(), " distance ", distanceseq, " 1st ", start, endfirst, " 2nd ", startsend, end)
              highlighted = f"{seq[:start]}**{seq[start:end]}**{seq[end:]}"
              if  distanceseq < distance and distanceseq > 0:
                matchcount=matchcount+1
                # st.code (seq, language=None, wrap_lines=True)
                st.markdown (highlighted)
              if "*" not in regex_part and distanceseq == 0:
                matchcount=matchcount+1
                st.markdown (highlighted)
         st.write("Done", matchcount, " matches found and ", matchmore, " found greater than distance ")
 
    elif choice == "Protein":
        st.subheader("Protein Analysis")
        st.write("enter a sequence ID from uniproto OR a sequence in the text box")
        seqid = st.text_input ("Enter ID Sequence", value = "ID")
        seq = st.text_area ("Enter DNA Sequence")
        seqvalid(seq)
        dna_seq = SeqRecord(
           Seq(seq),
           id=seqid,
           name="seq name",
           description="seq desc",
        )
        if len(seq) > 5 or seqid != "ID":
          proteinanal(seq, seqid)

            #  st.write(dna_seq) 
            #  streamlit sample blog.jcharistech.com
            #  details = st.radio("Details", ("Description", "Sequence"))
            # Top Most Common Amino

    elif choice == "DotPlot":
        st.subheader("Generate Dot Plot For Two Sequences")

        seq1 = st.text_area("Enter sequence 1")
        seqvalid (seq1)
        seq2 = st.text_area("Enter sequence 2")
        seqvalid (seq2)
        
        cus_limit = st.slider(
                "Select Max number of Nucleotide", value=50, max_value=200, min_value=10)

        if st.button("Dot Plot"):
                st.write(
                    "Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                fig, ax = plt.subplots()
                ax = dotplotx(seq1[0:cus_limit], seq2[0:cus_limit])
                st.pyplot(fig)

    elif choice == "About":
        st.subheader("About")


if __name__ == '__main__':
    main()
