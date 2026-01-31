import pubnetworks as pb


"""
READ THIS BEFORE GOING FURTHER

NETWORK                                             FUNCTION
1. XYZ Network             --->            pb.xyz("input chem_net file path", {"element1", "element2"},metadata = N)

        xyz(path, allowed_elements, metadata='N')

        Filters reactions from a xyz-format chemical network file based on allowed elements.

        Arguments:
            path (str): Path to input chemical network file.
            allowed_elements (set of str): Elements to keep (case-insensitive). Secondary species
                like H2O, OH are included if composed of these elements. Example: {"H","O","He"}.
            metadata (str, optional): 'Y' to include metadata, 'N' to exclude (default).

        Output:
            Writes a .dat file named xyz_output_{ELEMENTS}.dat in the current directory.
            Preserves the full original reaction line.
        Example:
            import pubnetworks as pb

            allowed = {"H", "O", "He"}
            pb.xyz(
                path="Input_SNCHO_VULCAN.txt",
                allowed_elements=allowed,
                metadata='Y'
            )


"""

"""

    pb.xyz() 
    
    Here xyz can be:
    ================================================
     xyz                   Published Network
    ================================================
    1. vulcan              (Tsai+ 2024b).
    2. hu                  (Hu R., Seager S.+ 2012)
    3. moses               (Moses+ 2011)
    4. stand               (Rimmer P.B+ 2020)
    5. velliet_venot       (Veillet+  2024)
    6. agundez             (Agudez 2025)
    =================================================


"""



elements= {"O", "H", "He"}

elements_v= { "Cl", "e"}

"""
# For Moses data acquired via Shami paper on SO2
pb.moses_vII(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/MOSES_2011_data",
    allowed_elements=elements, metadata = 'N'
)
"""

pb.vulcan(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Input_VULCAN_fullnet.txt",
    allowed_elements=elements,
      metadata = 'N'
)


pb.hu(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Input_HU2012.txt",
    allowed_elements=elements, metadata = 'N'
)

pb.moses_vII(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/MOSES_2011_data",
    allowed_elements=elements, metadata = 'N'
)


pb.stand(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Input_STAND2020.dat",
    allowed_elements=elements_v, metadata = 'N'
)

pb.velliet_venot(
    path="/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/VELLIET_VENOT_2025",
    allowed_elements=elements, metadata = 'N'
)

pb.agundez(
    path = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Input_Agundez2025.dat",
    allowed_elements = elements, metadata = 'N')

