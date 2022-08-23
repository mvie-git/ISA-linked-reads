import re
import fileinput
import sys

# Function to replace the wrong chromosome names with the right ones
def replace_all(text, dic):
    text=text.upper()
    for key, value in dic.items():
            text = re.sub("SN:"+key.upper()+'\t', "SN:"+value+'\t', text)
    return text

d = {
         "chr1": "chr1", "chr2": "chr10", "chr3" :"chr11", "chr4": "chr12",
         "chr5": "chr13", "chr6": "chr14", "chr7" :"chr15", "chr8": "chr16",
         "chr9": "chr17", "chr10": "chr18", "chr11" :"chr19", "chr12": "chr2",
         "chr13": "chr3", "chr14": "chr4", "chr15" :"chr5", "chr16": "chr6",
         "chr17": "chr7", "chr18": "chr8", "chr19" :"chr9", "chr20": "chrX",
         "chr21": "chrY"
     }


#d = {
#        "chr1": "chr1", "chr2": "chr12", "chr3" :"chr13", "chr4": "chr14",
#        "chr5": "chr15", "chr6": "chr16", "chr7" :"chr17", "chr8": "chr18",
#        "chr9": "chr19", "chr10": "chr2", "chr11" :"chr3", "chr12": "chr4",
#        "chr13": "chr5", "chr14": "chr6", "chr15" :"chr7", "chr16": "chr8",
#        "chr17": "chr9", "chr18": "chr10", "chr19" :"chr11", "CHRX": "chrX",
#        "CHRY": "chrY"
#    }


# file = open("sample.txt", "r")
# lines = file.readlines()
# file.close()
#
# for i in range(1,len(lines)):
#     lines[i]=replace_all(lines[i], d)
#
# file = open("sample.txt", "w")
# file.writelines(lines)
# file.close()

i=0 # Adapter si il y a un entête

# fileinput() module: iterates over the lines of all files listed in sys.argv[1:]
for line in fileinput.input(inplace=True):
    # if i>0 and i<22:
    if line.startswith("@SQ"):
        line = replace_all(line, d)
    sys.stdout.write(line)
    i=i+1
#    print("Ligne numéro : " + str(i) + " (début de la numérotation à 1)")
