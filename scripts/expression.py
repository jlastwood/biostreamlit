import re
from os.path import join
import streamlit as st

def expression(patternorig):

#  build the expressions for the search
 text = "This is a long string where ACTC     GNWYS appear."
# patternc = r"ACTCGNWYS"
# pattern = "[A|C]CTC[GA]NWYS"

# Example: search for "ABC...XYZ" where they can be separated by any characters

# Build a regex that allows the pattern to be split at any position
 regex_parts = []

 pattern = re.sub(r'[#*. ()0123456789]', '', patternorig)
 st.write("Clean search pattern", pattern, " original ", patternorig)

# Regex: match either [ ... ] groups OR single uppercase letters
 tokens = re.findall(r"\[[^\]]+\]|[A-Za-z]", pattern)
 strlength = len(tokens)
 # st.write(len(tokens), tokens)

#for i in range(1, len(pattern)):
#    left = str(pattern[:i])
#    right = str(pattern[i:])
#    regex_parts.append(f"{left}.*{right}")

 for i in range(1, len(tokens)):
    left = "".join(tokens[:i])
    right = "".join(tokens[i:])
    regex_parts.append(f"{left}.*{right}")
    # st.write (left,right)

# add the string for exact matches
 regex_parts.append(pattern)

# search for the pattern as a single string
 # st.write("there are ", len(regex_parts), " iterations and the string is ", strlength)
 # st.write(", ".join(regex_parts))
 # cannot use re.compile because it has a limit
 # regex = re.compile("|".join(regex_parts))
 # st.write(regex, ))
 return (regex_parts, strlength)

#  this is for testing
#match = regex.search(text)
#if match:
#    print("Found:", match.group(), match.start(), match.span())
#else:
#    print("Not found")
