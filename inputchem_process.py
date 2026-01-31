import inputchem_backend as icb


## Identifies the number of unique species in the final generated output file 
## filename: uniq_reactions_A_B_C.dat

path_agm_h_he_o = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/uniq_reactions_rateconst_H_HE_O.dat"

#n, outp = icb.write_unique_species_list(path_agm_h_he_o)
#print(n, outp)

##  Sort a .dat reaction file by number of reactants (ascending).

##--------------------------------------- SORTING ACC TO THE NO. OF REACTANTS ----------------------------
## INPUT FILE: uniq_reactions_{elements}.dat, uniq_reactions_rateconst_metadata_{elements}.dat
## This input is the final output of the chem_xmatch.py
## Gives sorted reactions acc to number of reactants
## This sort function is important as well...don't leave it as something optional 
## IMP: important if you want to do sort on all reactions including the ones present in Antonios
## Otherwise if you want to sort on the net reactions that you will add then the VII is fine
## The reaction order sorting works only on the removal of antonio reactions
## OUTPUT FILE :: sorted_{inputfilename}.dat
##----------------------------------------------------------------------------------------------------------

rxn_netpath = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/uniq_reactions_metadata_H_HE_O.dat"
krxn_netpath = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/uniq_reactions_rateconst_metadata_H_HE_O.dat"

def sort_final_rxns(s = "N"):
    if s == "Y":
        icb.sort_rxns_reactant(
            input_path=rxn_netpath,     
            kpath=krxn_netpath)
    else:
        print("Give argument \"Y\" to the function")

sort_final_rxns()


## =========================================================================================================

#sorted_rxnpath = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/sorted_uniq_reactions_H_HE_O.dat"
#sorted_krxnpath = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/sorted_uniq_reactions_rateconst_metadata_H_HE_O.dat"

## =========================================================================================================

## ========================================== REMOVAL OF REACTIONS ALREADY PRESENT IN CHEMISTRY.DAT ==============
## INPUT ::  uniq_reactions_rateconst_metadata...dat , only works with metadata files
## OUTPUT :: netrxn_addition_uniq_reactions_rateconst/metadata....dat
## TWO MODES :: ksort = Y or ksort = N, operates on diff files and gives diff output files
## ----------------------------------------------------------------------------------------------------------

"""
    Two modes:

    A) ksort == 'Y' (default): stacked file mode (reaction every kstack lines)
       - reads from `kpath`
       - each block has `kstack` lines
       - the reaction line is the FIRST line of each block
       - keep the WHOLE block only if the LAST column of the reaction line == '0000'
       - output: netrxn_addition_<kpath.name>.dat
       - if no blocks pass, prints:
         "no 0000 entry in last column. Please check the files"
       - if clean_metadata == 'Y', ALSO write a "clean" file that contains ONLY the
         kept reaction lines (first line of each kept block), using the same meaning
         of clean_metadata as in mode B(just one change about keeping the whole line):
           * skip first line only if it starts with 'Reaction'
           * write exactly the whole line 
         Output name: netrxn_addition_<kpath.name with leading 'metadata_' removed>

    B) ksort == 'N': original one-line-reaction table mode
       1) Write a filtered file that removes rows where agm2007 != '0000'.
          Output: netrxn_addition_<inputfilename>

       2) If clean_metadata == 'Y' (default), also write a "clean" file built from
          the ORIGINAL INPUT, but ONLY for the rows that pass the agm2007 filter.
          - skips the first line ONLY if it starts with 'Reaction'
          - writes exactly 50 characters per kept line (truncate/pad)
          Output name: netrxn_addition_<inputfilename with leading 'metadata_' removed>
"""

icb.rm_nonzero_agm(
    input_path = rxn_netpath,
    clean_metadata = 'Y',
    ksort = 'Y',
    kpath = krxn_netpath,
    kstack = 9)


## ============================================ SORTING REACTIONS AS PER AGM CHEMISTRY.DAT ============================
## INPUT :: The output files produced from the rm_nonzero_agm; netrxn_addition_...dat
## OUTPUT :: agmfinal_{inputfilename}.dat

## NO DIFFERENCE :: between netrxn_addition....dat and agmfinal_netrxn_addition....dat other than the sorting of reactions




"""
Two modes:

    A) kmode == 'N' (default): one reaction per line
       - Reaction is the first 50 characters of each line
       - Tail is everything after character 50 and is preserved and written back
       - Output filename:
           agmfinal_<inputfilename>
           OR agmfinal_<metadata_inputfilename> (keeps name exactly)
         (i.e., output is always agmfinal_<in_path.name>)

    B) kmode == 'Y': stacked blocks (one reaction per kstack lines), separated by lines starting with '##########'
       - Input taken from `kpath`
       - Each reaction block has:
           * reaction header line (first line of block; reaction is chars 0:50)
           * the next (kstack-1) metadata lines
           * typically followed by a separator line '##########' (kept with the block)
       - Sorting is computed from the reaction text in the header line (first 50 chars)
       - When writing, the entire original block (header+metadata+separator if present) is written back as-is.

    Sorting order (same logic as before, multiplicity counts):
    1) Unimolecular: exactly 1 reactant token on LHS (HV not present)
    2) Photolysis: 1 reactant token + HV on LHS (exactly two tokens total, one is HV)
    3) Bimolecular: exactly 2 reactant tokens on LHS (counting multiplicity),
       excluding M and T as reactants, and excluding HV (already handled)
    4) Termolecular: exactly 3 reactant tokens on LHS (counting multiplicity) and includes M or T
    5) Any reaction containing M as reactant token (not already captured)
    6) Any reaction containing T as reactant token (not already captured)
    7) Anything else

"""
icb.sort_reactions_order(
    path = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/netrxn_addition_uniq_reactions_rateconst_H_HE_O.dat",
    kmode = "N",
    kpath = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/netrxn_addition_uniq_reactions_rateconst_metadata_H_HE_O.dat",
    kstack = 9
) 


