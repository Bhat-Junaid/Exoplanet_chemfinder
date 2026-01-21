import inputchem_backend as icb


## Identifies the number of unique species in the final generated output file 
## filename: uniq_reactions_A_B_C.dat

path_agm_h_he_o = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/uniq_reactions_H_HE_O.dat"

#n, outp = icb.write_unique_species_list(path_agm_h_he_o)
#print(n, outp)

##  Sort a .dat reaction file by number of reactants (ascending).

##--------------------------------------- SORTING ACC TO THE NO. OF REACTANTS ----------------------------
## INPUT FILE: uniq_reactions_{elements}.dat, uniq_reactions_rateconst_metadata_{elements}.dat
## This input is the final output of the chem_reprocess.py
## Gives sorted reactions acc to number of reactants

rxn_netpath = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/uniq_reactions_H_HE_O.dat"
krxn_netpath = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/uniq_reactions_rateconst_metadata_H_HE_O.dat"

def sort_final_rxns(s = "N"):
    if s == "Y":
        icb.sort_rxns_reactant(
            input_path=rxn_netpath,     
            kpath=krxn_netpath)
    else:
        print("Give argument \"Y\" to the function")

sort_final_rxns()


sorted_rxnpath = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/sorted_uniq_reactions_H_HE_O.dat"
sorted_krxnpath = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/sorted_uniq_reactions_rateconst_metadata_H_HE_O.dat"

icb.rm_nonzero_agm(
    input_path = sorted_rxnpath,
    clean_metadata = 'Y',
    ksort = 'Y',
    kpath = sorted_krxnpath,
    kstack = 9)

