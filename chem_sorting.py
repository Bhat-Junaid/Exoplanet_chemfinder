import pubnetworks as pb


elements= {"H", "He", "O"}

####################################### PART I ###################################
##### CHEMICAL SCHEME GENERATION INCLUDING A SET OF CHOSEN ELEMENTS ##############
##### EXAMPLE: elements = {"H", "He", "O"}; H - He - O CHEMICAL NETWORK ##########
##### IMPORTANT: metadata = 'N' ; in the functions used in this part- always #####
#################################################################################

pb.vulcan(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Input_VULCAN_fullnet.txt",
    allowed_elements=elements)


pb.hu(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Input_HU2012.txt",
    allowed_elements=elements)

pb.moses(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/MOSES_2011_data",
    allowed_elements=elements)

pb.stand(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Input_STAND2020.dat",
    allowed_elements=elements)

pb.velliet_venot(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/VELLIET_VENOT_2025",
    allowed_elements=elements)

pb.agundez(
    path = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Input_Agundez2025.dat",
    allowed_elements = elements)



####################################### PART II ###############################

"""
THIS PART IS BROKEN DOWN INTO TWO SUB-PARTS:

A) CHEMICAL REACTION FILE REPROCESSING: To a homogenous nomeclature of the 
reaction networks. Standardization of all different reaction writing proc-
dures used by different authors.

INPUT: The CHEMICAL NETWORK FILES produced in PART I (without metadata!!!).
OUTPUT: Files with name "reprocess_{input_file_name}.dat".
These files contain the reactions in a one standardized uniform format.

IMP :: GENERAL CHARACTERISTICS OF THE FILES PRODUCED
- HV = Used to denote light (h\nu) across all the files.
- BRACKETS () = Used to write the excited states. 
                e.g: O(3P), O(1S), O(X2Pig), CH2(*) etc
- ^+ or ^- or ^-- = In general the caret ^ is used to denote ions in the file.

- FOR STAND network file processing only; OXYRANE (in input file) = C2H4O (its
chemical formula) in the output file
"""

pb.reprocess_velliet_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/velliet_venot_output_H_HE_O.dat")
pb.reprocess_agundez_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Agundez_output_H_HE_O.dat")
pb.reprocess_hu_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Hu_output_H_He_O.dat")
pb.reprocess_vulcan_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/VULCAN_output_H_He_O.dat")
pb.reprocess_moses_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/MOSES_output_H_He_O.dat")
pb.reprocess_stand_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/stand_output_H_He_O.dat")


