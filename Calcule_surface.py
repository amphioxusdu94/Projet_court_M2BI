"""A Python code that calculates the solvent exposed surface area of a protein.

Return:
- A value for all the protein
- A first file containing exposed surface for each non hydrogen atom
- A second file containing exposed surface for each Amino acid
"""

__author__ = "Aurélien Desvilles"
__contact__ = "desvillesaurelien@gmail.com"
__institution__ = "Université Paris Cité"
__date__ = "13/09/2023"
__version__ = "Python 3.9.13"


import numpy as np
import pandas as pd
import math
from Bio import PDB


"""Van der Waals radius of an oxygen atom (from a water molecule)."""
probe_radius = 1.52  
  
"""Van der Waals radius of each main atoms found in biological molecules (exept hydrogen)."""
VdW_radius = {"C" : 1.70, "N" : 1.55, "O" : 1.52, "F" : 1.47, "P" : 1.80, "S" : 1.80,
 "Cl" : 1.75, "Cu" : 1.40} 

"""Number of points for each sphere."""
nb_points = 20



"""
Functions
"""

def parse_PDB(file_pdb):
    """Parse pdb file and return for each no hydrogen :
    -atome name
    -atome type 
    -Amino acid type
    -amino acid number
    -center coordinates
        """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", file_pdb)

    atom_name = []
    atom_type = []
    residue_type = []
    residue_number = []
    coords_center = []
    protein_name = None

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if residue.resname == 'HOH':
                        continue 
                    if atom.element == 'H':
                        continue
                    if atom.element not in VdW_radius.keys():
                        continue
                    atom_name.append(atom.get_name())
                    atom_type.append(atom.element)
                    residue_type.append(residue.resname)
                    residue_number.append(residue.id[1])
                    coords_center.append(atom.coord)
                    
                    
    with open(file_pdb, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('HEADER'):
                protein_name = line[10:50].strip()                
    # Create DataFrame 
    data = {
        'Atom Name': atom_name,
        'Atom Type': atom_type,
        'Residue': residue_type,
        'Residue Number': residue_number,
        'Coordinates': coords_center
    }
    df = pd.DataFrame(data)
    df.name = protein_name
    return df



def compute_sphere_area(radius):
    """ Compute area for a spheric tool (for example an atom) """
    area = 4 * math.pi * (radius ** 2)
    return area



def generate_fibonacci_sphere_points(num_points, radius, center):
    """ Uniformly generate a number of points distributed around a sphere of predefined center and radius.
    Refereces: https://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
    """
    coords = []
    i = np.arange(0, num_points) 
    phi = (np.sqrt(5.0) + 1.0) / 2.0  # Golden ratio
    longitudes = 2.0 * np.pi * i / phi
    latitudes = np.arccos(1 - 2 * (i + 0.5) / (num_points))
    # Convert spherical coordinates to Cartesian
    for i in range(num_points):
        x = radius * np.sin(latitudes[i]) * np.cos(longitudes[i]) + center[0]
        y = radius * np.sin(latitudes[i]) * np.sin(longitudes[i]) + center[1]
        z = radius * np.cos(latitudes[i]) + center[2]
        coords.append([x, y, z])
    
    return coords



def distance_between_points(point1, point2):
    """ Compute euclidean distance between two 3D points. """ 
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2 + (point1[2] - point2[2])**2)



def points_without_proximity(points1, points2, distance_threshold):
    """Compute the points of a first list for which the
    distance with all the points of a second list is less 
    than a threshold distance.
    """
    points_without_proximity = []

    for point1 in points1:
        proximity_found = False
        for point2 in points2:
            distance = distance_between_points(point1, point2)
            if distance < distance_threshold:
                proximity_found = True
                break
        
        if not proximity_found:
            points_without_proximity.append(point1)
    
    return points_without_proximity



"""
Main Program
"""

if __name__ == "__main__":
    
    print("\n Hello")
    #File containing atom coordinates of human isulin.For example.
    file = "1h59.pdb"
    #Parsing file
    Atoms=parse_PDB(file) #A data frame containig eache no hydrogen atom features
    print(f"\n You are studying a protein named {Atoms.name}.")
    
    #Add empty columns intended to be filled in later
    Atoms['Atomic surface'] = None
    Atoms['Points on sphere'] = None
    Atoms['Points on sphere exposed'] = None
    Atoms['Exposed atomic surface (%)'] = None
    Atoms['Exposed atomic surface (Å2)'] = None

    #Creating sphere for each atom.
    print('\n Creating sphere for each atom.')
    for index, row in Atoms.iterrows():
        center = row['Coordinates']
        for key in VdW_radius:
            if row['Atom Type'] == key:
                radius = VdW_radius[key]
                sphere_area = compute_sphere_area(radius)
                Atoms.at[index,'Atomic surface'] = sphere_area
                points_on_sphere = generate_fibonacci_sphere_points(nb_points, radius + probe_radius, center)
                Atoms.at[index,'Points on sphere'] = points_on_sphere
    
    #Create list of each point set, it make the next step easier.
    list_points_on_sphere = []
    for i in Atoms['Points on sphere']:
        list_points_on_sphere.append(i)
    
    #Calculate for each sphere, points for which the distance to any point on the other spheres is less than the probe radus.
    print('\n Compare distances between points on spheres. It could take time...')
    for i in range(0,len(list_points_on_sphere)):
        for j in range(0,len(list_points_on_sphere)):
            if i != j:
                list_points_on_sphere[i] = points_without_proximity(list_points_on_sphere[i], list_points_on_sphere[j], probe_radius)
        Atoms.at[i,"Points on sphere exposed"] = list_points_on_sphere[i]

    #Calculate relative exposed surface for each atom
    for i in range(0,len(Atoms["Points on sphere"])):
        purcentage = 100 * len(Atoms.loc[i,"Points on sphere exposed"]) / len(Atoms.loc[i,"Points on sphere"])
        Atoms.at[i,'Exposed atomic surface (%)'] = purcentage
 
    #Calculate absolute exposed surface for each atom   
    for i in range(0,len(Atoms)):
        Atoms.at[i,'Exposed atomic surface (Å2)'] = Atoms.loc[i,"Atomic surface"] * Atoms.loc[i,"Exposed atomic surface (%)"] / 100
    
    #Calculate relative exposed surface for all protein.
    total_relative_exposed_surface = Atoms["Exposed atomic surface (%)"].sum() / len(Atoms)
    print(f"\n{total_relative_exposed_surface:.2f}% of the protein surface is exposed to solvent.")

    #Calculate absolute exposed surface for all protein.
    total_surface = Atoms["Exposed atomic surface (Å2)"].sum()
    print(f"\nThis is equivalent to an area of {total_surface:.2f} Å2.")
    
    
    #Create a first file containing exposed surface for each non hydrogen atom
    print("\n Generating file containing exposed surface for each non hydrogen atom... ")
    exposed_surface_per_atom = Atoms[['Atom Name','Atom Type','Residue','Residue Number','Exposed atomic surface (%)','Exposed atomic surface (Å2)']]
    file_name = file.split('.')[0]
    results_file_name_atoms = file_name + "_exposed_surface_atoms_" + str(nb_points) + "points.txt"
    exposed_surface_per_atom.to_csv(results_file_name_atoms, sep='\t', index=False, header=True)
    
    
    #Create a second file containing exposed surface for each Amino acid
    print("\n Generating file containing exposed surface for each residue... ")
    #create second data frame
    absolute_surface_per_AA = Atoms.groupby('Residue Number')['Exposed atomic surface (Å2)'].sum()
    type_residu = Atoms.groupby('Residue Number')['Residue'].first()
    sum_relative = Atoms.groupby('Residue Number')['Exposed atomic surface (%)'].sum()
    number_elements_col = Atoms.groupby('Residue Number')['Exposed atomic surface (%)'].count()
    relative_surface_per_AA = sum_relative / number_elements_col

    exposed_surface_per_resudual = pd.DataFrame({'Type': type_residu,
                            'Relative surface exposed (%)': relative_surface_per_AA,
                            'Absolute surface exposed (Å2)': absolute_surface_per_AA})
    
    results_file_name_residues = file_name + "_exposed_surface_AA_" + str(nb_points) + "points.txt"
    exposed_surface_per_resudual.to_csv(results_file_name_residues, sep='\t', index=True, header=True)

    #End word
    print("n\ Work completed. See yoo soon :)")
