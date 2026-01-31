import pubnetworks as pb
import os
import re
from collections import OrderedDict

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

"""
pb.moses(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/MOSES_2011_data",
    allowed_elements=elements)
"""
pb.moses_vII(
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

###################################### PART II.A #############################

A) CHEMICAL REACTION FILE REPROCESSING: To a homogenous nomeclature of the 
reaction networks. Standardization of all different reaction writing proc-
dures used by different authors.

INPUT: The CHEMICAL NETWORK FILES produced in PART I (without metadata!!!).
OUTPUT: Files with name "reprocess_{input_file_name}.dat".
These files contain the reactions in a one standardized uniform format.

"""



pb.reprocess_velliet_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/velliet_output_H_HE_O.dat")
pb.reprocess_agundez_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Agundez_output_H_HE_O.dat")
pb.reprocess_hu_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Hu_output_H_He_O.dat")
pb.reprocess_vulcan_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/VULCAN_output_H_He_O.dat")
pb.reprocess_mosesvii_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/MOSESII_output_H_He_O.dat")
pb.reprocess_stand_file("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/stand_output_H_He_O.dat")






"""
###################################### PART II.B #############################

B) UNIQUE CHEMICAL REACTION PRODUCTION: Produces all the reactions that occur 
at least once in the six networks (in chosen element network e.g H-He-O).

INPUT: The CHEMICAL NETWORK FILES produced in PART II.A (reprocessed).
OUTPUT: 3 files: 
i) Uniq_reactions_E1_E2_E3_E4.dat :: contains all the unique reactions as 1 col
                                     format dat file.

ii) Uniq_reactions_metadata_E1_E2_E3_E4.dat 
                                  :: contains all the unique reactions in the 1
                                  st col and then in subsquent cols the line no.
                                  in which the reaction was found in the input 
                                  reprocessed line [Good for sanity checks and 
                                  cross linking.]

iii) Uniq_reactions_longtable_E1_E2_E3_E4.dat 
                                  :: contains 6 columns with header correspondi-
                                  ng to each network (Hu, Agundez, VULCAN, STAND
                                  ,Velliet_Venot, Moses) and under them the reac-
                                  tions present in the respective networks.
                                  All the similar/repeated reactions like in the 
                                  same row, just under diff columns.               
                                     
same with k = rate constant files as well but the
"""



path_vv = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_velliet_output_H_HE_O.dat"
path_agundez= "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_Agundez_output_H_HE_O.dat"
path_hu = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_Hu_output_H_He_O.dat"
path_vulcan = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_VULCAN_output_H_He_O.dat"
path_moses = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_MOSESII_output_H_He_O.dat"
path_stand = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_stand_output_H_HE_O.dat"
path_agm = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_agm2007_output_CO2.dat"

reprocessed_file_paths = [path_vv, path_agundez, path_hu, path_vulcan, path_moses, path_stand, path_agm]



pb.build_unique_reactions(reprocessed_file_paths, long_table="Y")


k_path_vv = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_rate_const_velliet_output_H_HE_O.dat"
k_path_agm = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_rate_const_agm2007_output_CO2.dat"
k_path_agundez = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_rate_const_Agundez_output_H_HE_O.dat"
k_path_hu = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_rate_const_Hu_output_H_He_O.dat"
k_path_vulcan = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_rate_const_VULCAN_output_H_He_O.dat"
k_path_moses = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_rate_const_MOSESII_output_H_He_O.dat"
k_path_stand = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reprocess_rate_const_stand_output_H_He_O.dat"

reprocessed_kfile_paths = [k_path_hu, k_path_agundez, k_path_vulcan, k_path_vv, k_path_moses, k_path_stand, k_path_agm]


#     Outputs:
#  1) uniq_reactions_rateconst_<TAG>.dat

pb.build_unique_rxns_ratek(reprocessed_kfile_paths)