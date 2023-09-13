# Projet_court_M2BI
A Python script which implements Shrake-Rupley algorithm to calculate the solvent exposed surface of a protein.
Here it will do that for the human insulin for which the atomic coordinates are contained in 1h59.pdb file.

Steps are:
1) Extraction of each atom coordinates.
2) From each contained atom, creation of a cloud of points uniformly
on the surface of a sphere centered on the atom. 
3) Only keep points without contact any point of any sphere.
4) Deduce the solvent exposed surface.

Return:
- A absolute (Å) and relative (%) exposed surface for all the protein.
- A first file containing absolute (Å) and relative (%) exposed surface for each non hydrogen atom.
- A second file containing (Å) and relative (%) exposed surface for each Amino acid.

You can modify the script and change the points for each sphere number.
You should nee to install a Conda environment (see Projet_court_env.yml).
You can run it like a classic Python scrip: 
                                        $: python Exposed_surface.py
                        
