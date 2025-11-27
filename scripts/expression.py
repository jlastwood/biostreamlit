import re
from os.path import join

def expression(pattern):

#  build the expressions for the search
 text = "This is a long string where ACTC     GNWYS appear."
# patternc = r"ACTCGNWYS"
# pattern = "[A|C]CTC[GA]NWYS"

# Example: search for "ABC...XYZ" where they can be separated by any characters

# Build a regex that allows the pattern to be split at any position
 regex_parts = []
# first add the string
 regex_parts.append(pattern)

# Regex: match either [ ... ] groups OR single uppercase letters
 tokens = re.findall(r"\[[^\]]+\]|[A-Za-z]", pattern)
 #print(len(tokens), tokens)

#for i in range(1, len(pattern)):
#    left = str(pattern[:i])
#    right = str(pattern[i:])
#    regex_parts.append(f"{left}.*{right}")

 for i in range(1, len(tokens)):
    left = "".join(tokens[:i])
    right = "".join(tokens[i:])
#   print (left,right)
    regex_parts.append(f"{left}.*{right}")

# search for the pattern as a single string
 regex = re.compile("|".join(regex_parts))
 return (regex)

#  this is for testing
#match = regex.search(text)
#if match:
#    print("Found:", match.group(), match.start(), match.span())
#else:
#    print("Not found")
