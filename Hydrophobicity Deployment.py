import math

import Bio.pairwise2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment


def distance(x1, y1, z1, x2, y2, z2):
    d = math.sqrt(math.pow(x2 - x1, 2) +
                  math.pow(y2 - y1, 2) +
                  math.pow(z2 - z1, 2) * 1.0)
    #print("Distance is ")
    return(d)

def main():
    #Set maximum distance threshold
    threshold = 20

    #Set a dictionary to convert three letter amino acid code to one letter
    aminoacid_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'AGL': 'Q', 'BGL': 'Q'}

    hydrophobicity_dict = {"D": -3.50, "E": -3.50, "R": -4.50, "H": -3.20, "K": -3.90, "G": -0.40, "P": -1.60,
                           "A": 1.80,
                           "N": -3.50, "C": 2.50, "Q": -3.50, "I": 4.50, "L": 3.80, "M": 1.90, "F": 2.80,
                           "S": -0.80, "T": -0.70, "W": -0.90, "Y": -1.30, "V": 4.20, "-": 0, "X": 0}

    #Set up dictionaries to store amino acid coordinates with amino acid id as key
    covid_data = {}
    ace2_data = {}
    spatiallycloseaminoacids = {}
    #aminoacid_id = "" #parsed value
    #covid_data[aminoacid_id] = []
    #covid_data["ALA_32"].append(float(x), float(y), float(z))
    #ace2_data[aminoacid_id] = []
    #if parsing for if key exists in dictionary
    #else just append

    #Open PDB file
    filename = "6lzg.pdb"
    with open(filename, "r") as file:
        line = file.readline()
        for line in file:
            splitline = line.split(" ")
            #look at key positions in the file
            #fixed length string
            #use line instead of splitline
            #line.substring[si, ei]
            #if splitline[0] == "ATOM" and splitline[2] == "C":

            #Select only lines that correspond to atoms
            if splitline[0] == "ATOM":
                splitline = list(filter(None, splitline))
                #print(splitline)
                aminoacid_id = splitline[3] + str(splitline[5])
                #print(aminoacid_id)

                #If atom belongs to chain B
                if splitline[4] == "B":
                    #print(splitline)
                    #print(aminoacid_id)
                    if aminoacid_id not in covid_data:
                        x = splitline[6]
                        y = splitline[7]
                        z = splitline[8]
                        #print(x)
                        #print(y)
                        #print(z)
                        aminoacid_distance = []
                        #aminoacid_distance.append(float(x), float(y), float(z))
                        aminoacid_distance.append(float(x))
                        aminoacid_distance.append(float(y))
                        aminoacid_distance.append(float(z))
                        covid_data[aminoacid_id] = aminoacid_distance
                else:
                    if aminoacid_id not in ace2_data:
                        x = splitline[6]
                        y = splitline[7]
                        z = splitline[8]
                        #ace2_data[aminoacid_id].append(float(x), float(y), float(z))
                        aminoacid_distance = []
                        #aminoacid_distance.append( float(x), float(y), float(z))
                        aminoacid_distance.append(float(x))
                        aminoacid_distance.append(float(y))
                        aminoacid_distance.append(float(z))
                        ace2_data[aminoacid_id] = aminoacid_distance
        print(covid_data)
        print(ace2_data)


        #covid_aminoacid = list(filter(None, covid_data.keys()))
        #ace2_aminoacid = list(filter(None, ace2_data.keys()))
        #print(covid_aminoacid)
        #print(ace2_aminoacid)
        p1 = []
        p2 = []
        for aminoacid in covid_data:
            #print("Covid Amino Acid: " + str(aminoacid) + " X: " + str(covid_data[aminoacid][0]) + " Y: " + str(
             #   covid_data[aminoacid][1]) + " Z: " + str(covid_data[aminoacid][2]))
            for aminoacidace in ace2_data:
            #    print("Ace2 Amino Acid: " + str(aminoacidace) + " X: " + str(ace2_data[aminoacidace][0]) + " Y: " + str(
             #       ace2_data[aminoacidace][1]) + " Z: " + str(ace2_data[aminoacidace][2]))
                dist = distance(covid_data[aminoacid][0], covid_data[aminoacid][1], covid_data[aminoacid][2], ace2_data[aminoacidace][0], ace2_data[aminoacidace][1], ace2_data[aminoacidace][2])
               # print("Comparing: " + str(aminoacid) + " to " + str(aminoacidace) + " distance: " + str(dist))
                if dist <= threshold:
                    spatiallycloseaminoacids[aminoacid] = [aminoacidace, dist]
    print(spatiallycloseaminoacids)
    spatiallycloseaminoacids_covid = list(spatiallycloseaminoacids.keys())
    print(spatiallycloseaminoacids_covid)
    covid_sequence_three = []
    for aminoacidcovid in spatiallycloseaminoacids_covid:
        covid_sequence_three.append(aminoacidcovid[0:3])
    print(covid_sequence_three)
    covid_sequence = []
    for aminoacid_covid in covid_sequence_three:
        covid_sequence.append(aminoacid_dict[aminoacid_covid])
    print(covid_sequence)
    #Turn list of amino acids into a string
    separator = ""
    covid_sequence_one = separator.join(covid_sequence)
    #print(covid_sequence_one)

    #Read variant fasta file
    filename_variant = "sequence BA2.fasta"
    variantsequence = SeqIO.read(filename_variant, "fasta")
    variantsequence = variantsequence.format("fasta")
    variantsequence = variantsequence.replace('\n','')
    variantsequence = variantsequence[82:]
    print(variantsequence)

    #Align reference and variant sequences
    alignment = Bio.pairwise2.align.globalxx(covid_sequence_one, variantsequence)
    alignment = format_alignment(*alignment[0])
    print(alignment)
    alignments = alignment.split('\n')
    reference_sequence = alignments[0]
    variant_sequence = alignments[2]
    print(reference_sequence)

    hydrophobicity_differencelist = []
    for i in range(len(reference_sequence)):
        variant_hydrophobicity = float(hydrophobicity_dict[variant_sequence[i-0]])
        reference_hydrophobicity = float(hydrophobicity_dict[reference_sequence[i-0]])
        hydrophobicity_difference = float(reference_hydrophobicity) - float(variant_hydrophobicity)
        hydrophobicity_differencelist.append(float(hydrophobicity_difference))
    print(hydrophobicity_differencelist)



if __name__ == "__main__":
    main()
    #Sort through relevant file lines
    #Calculate distance between each atom
    #Add atom to new file if it meets threshold