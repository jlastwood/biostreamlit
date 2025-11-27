from PIL import Image
import numpy as np
import io
from Bio.Seq import Seq
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
    activity = ['Intro', 'Upload Fasta', 'Expression', 'DNA', 'DotPlot', "About"]
    choice = st.sidebar.selectbox("Select Activity", activity)
    if choice == 'Intro':
        st.subheader("Intro")

    elif choice == 'Upload Fasta':
        st.subheader("Upload")
        fasta = st.file_uploader("Upload FASTA File", type=[".gz"])
        if fasta is not None:
            st.session_state.seq_dict = readfasta(fasta)
            st.session_state.seq_length = len(st.session_state.seq_dict)
            st.write(f"Length: {len(st.session_state.seq_dict)}")
            for sid, seq in islice(st.session_state.seq_dict.items(),5):
               st.write(f"{sid}: length {len(seq)}")

            #st.write(f"Length: {len(dna_seq)}")
            #st.write(dna_seq[:50] + "...")  # print first 200 bases for preview

    elif choice == 'Expression':
        st.subheader("Expression Search")
        st.write("Enter the expression to search.  You can use a simple format like ATCGNWYS.  You can also use [] in the expression, such as [A|C]TC[AC]WYS.  After entering the search press enter, the results of the expression are displayed. ") 
        exp_text = st.text_input("enter Expression to search")
        regex = expression(exp_text)
        st.write("This is the regular expression for search", regex)
        st.write("The number of sequences is ", st.session_state.seq_length, " If the sequences are 0, go to Upload Fasta.  To perform the search, click Analyse Sequence")
        if st.session_state.seq_length > 0:
         gobutton = st.button ("Analyse Sequence")
         matchcount = 0
         if gobutton:
           for sid, seq in islice(st.session_state.seq_dict.items(),1000000): 
             match = regex.search(seq)
             if match:
              st.write("Found:", sid, match.group(), match.start(), match.span())
              matchcount=matchcount+1
           st.write("Done", matchcount, " matches found")
 
    elif choice == "DNA":
        st.subheader("DNA Sequence Analysis")
        fasta = st.file_uploader("Upload FASTA File", type=[".fasta", ".fa"])

        if fasta is not None:
            dna_df = to_df(fasta)
            st.write(dna_df) 

            stringio0 = io.StringIO(fasta.getvalue().decode("utf-8"))
            dna_record = SeqIO.parse(stringio0, "fasta")
            
            details = st.radio("Details", ("Description", "Sequence"))

            for seq1 in dna_record:
               dna_seq = seq1.seq

            # Nucleotide Frequencies
            st.subheader("Nucleotide Frequency")
            dna_freq = Counter(dna_seq)
            st.write(dna_freq)
            adenine_color = st.color_picker("Adenine Color")
            thymine_color = st.color_picker("thymine Color")
            guanine_color = st.color_picker("Guanine Color")
            cytosil_color = st.color_picker("cytosil Color")

            if st.button("Plot Freq"):
                fig, barlist = plt.subplots()
                barlist = plt.bar(dna_freq.keys(), dna_freq.values())
                barlist[2].set_color(adenine_color)
                barlist[3].set_color(thymine_color)
                barlist[1].set_color(guanine_color)
                barlist[0].set_color(cytosil_color)

                st.pyplot(fig)

            st.subheader("DNA Composition")
            gc_score = utils.gc_content(str(dna_seq))
            at_score = utils.at_content(str(dna_seq))
            st.json({"GC Content": gc_score, "AT Content": at_score})

            # Nucleotide Count
            nt_count = st.text_input(
                "Enter Nucleotide Here",
                "Type Nucleotide Alphabet")
            st.write("Number of {} Nucleotide is ::{}".format(
                (nt_count), str(dna_seq).count(nt_count)))

            # Protein Synthesis
            st.subheader("Protein Synthesis")
            p1 = dna_seq.translate()
            aa_freq = Counter(str(p1))

            if st.checkbox("Transcription"):
                st.write(dna_seq.transcribe())

            elif st.checkbox("Translation"):
                st.write(dna_seq.translate())

            elif st.checkbox("Complement"):
                st.write(dna_seq.complement())

            elif st.checkbox("AA Frequency"):
                st.write(aa_freq)

            elif st.checkbox("Plot AA Frequency"):
                aa_color = st.color_picker("Pick An Amino Acid Color")
                # barlist = plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                # barlist[2].set_color(aa_color)
                fig, ax = plt.subplots()
                ax = plt.bar(aa_freq.keys(), aa_freq.values(), color=aa_color)
                st.pyplot(fig)

            elif st.checkbox("Full Amino Acid Name"):
                aa_name = str(p1).replace("*", "")
                aa3 = utils.convert_1to3(aa_name)
                st.write(aa_name)
                st.write("=====================")
                st.write(aa3)

                st.write("=====================")
                st.write(utils.get_acid_name(aa3))

            # Top Most Common Amino

    elif choice == "DotPlot":
        st.subheader("Generate Dot Plot For Two Sequences")
        seq_file1 = st.file_uploader(
            "Upload 1st FASTA File", type=[
                "fasta", "fa"])
        seq_file2 = st.file_uploader(
            "Upload 2nd FASTA File", type=[
                "fasta", "fa"])

        if seq_file1 and seq_file2 is not None:
            stringio1 = io.StringIO(seq_file1.getvalue().decode("utf-8"))
            dna_record1 = SeqIO.parse(stringio1, "fasta")
            stringio2 = io.StringIO(seq_file2.getvalue().decode("utf-8"))
            dna_record2 = SeqIO.parse(stringio2, "fasta")
            # st.write(dna_record1.composition)
            for seq1 in dna_record1:
               dna_seq1 = seq1.seq
            for seq2 in dna_record2:
               dna_seq2 = seq2.seq
            details = st.radio("Details", ("Description", "Sequence"))
            if details == "Description":
              for seq1 in dna_record1:
                st.write(seq1.description)
                st.write("=====================")
              for seq2 in dna_record2:
                st.write(seq2.description)
                st.write("=====================")
            elif details == "Sequence":
              for seq1 in dna_record1:
                st.write(seq1.seq)
                st.write("=====================")
              for seq2 in dna_record2:
                st.write(seq2.seq)
                st.write("=====================")

            cus_limit = st.number_input(
                "Select Max number of Nucleotide", 10, 200, 50)
            if st.button("Dot Plot"):
                st.write(
                    "Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                fig, ax = plt.subplots()
                ax = dotplotx(dna_seq1[0:cus_limit], dna_seq2[0:cus_limit])
                st.pyplot(fig)

    elif choice == "About":
        st.subheader("About")


if __name__ == '__main__':
    main()
