"""Convert the SPR-KKR exchange interactions to vampire.
    Needed are the seedname.sys, seedname_SCF.out, seedname_JXC_XCPLTEN_Jij.dat and seedname_JXC_XCPLTEN_Dij.dat files."""

import pandas as pd
import numpy as np


import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '../ucf_crop_small_values_vampire')
from ucf_crop_small_values_vampire import ucf_crop



# ========================== USER INPUT =============================

# SPR-KKR calculation directory
path = "H:/2D/Cr2Te3/bulk/sprkkr_imitating_experimental/Gr" #"H:/2D/Cr2Te3/bulk/transport/PBE/NKTAB_1000"

# the seed name in SPR-KKR (the name of the system)
system_name = "POSCAR"

include_dmi = False
f_electrons = False

crop_threshold = 0.1    # meV; set value > 0 to perform the crop post-processing (via the ucf_crop_small_values_vampire.py script; creates a new file)

# ===================================================================



def write_mat_file(path, system_name, mag_moments, elements):
    """Write the .mat file for vampire, given the magnetic moments and types of elements."""
    if len(mag_moments) != len(elements):
        raise ValueError("elements and mag_moments arrays need to be of the same length!")

    path_in = f"{path}/vampire.mat"

    with open(path_in, 'w') as fw:
        fw.write(f"material:num-materials = {len(elements)}\n")
        for i in range(len(elements)):
            magmom = float(mag_moments[i])
            fw.write(
                      (  "#---------------------------------------------------\n"
                         f"# Material {i+1} \n"
                         "#---------------------------------------------------\n"
                         f"material[{i+1}]:material-name={elements[i]}\n"
                         f"material[{i+1}]:damping-constant=1.0\n"
                         f"material[{i+1}]:atomic-spin-moment={abs(magmom):.8f} !muB\n"
                         f"material[{i+1}]:uniaxial-anisotropy-constant=0.0\n"
                         f"material[{i+1}]:material-element={elements[i]}\n"                        
                         f"material[{i+1}]:initial-spin-direction = 0.0,0.0,{abs(magmom)/magmom:.1f}\n"
                         f"material[{i+1}]:uniaxial-anisotropy-direction = 0.0 , 0.0, 1.0\n"
                         "#---------------------------------------------------\n"
                       )
            )
    print(".mat file written")


def get_mag_moments(path, system_name, n_atoms, f_electrons=False):
    """Parse magnetic moments from _SCF.out file."""
    mag_moments = []
    lines_to_skip = 4 if f_electrons==False else 5

    path_in = f"{path}/{system_name}_SCF.out"
    IT_flag, read_flag = False, False
    skipped_lines = 0

    with open(path_in, 'r') as fr:
        for line in fr:
            if IT_flag == True:
                skipped_lines += 1

            if "IT=" in line:
                IT_flag = True

            if skipped_lines > lines_to_skip:
                mag_moments.append(line.split()[4])
                skipped_lines = 0
                IT_flag = False

    print("mag_moments obtained")
    return mag_moments[-n_atoms:]



def get_structure_from_sys_sprkkr(path, system_name):
    """From the SPR-KKR's .sys file get the unit cell, number of atoms and their basis vectors."""
    path_in = f"{path}/{system_name}.sys"

    with open(path_in, 'r') as fr:
        latt_param_flag, primitive_vectors_flag, number_flag, basis_vectors_flag, elements_flag = False, False, False, False, False
        i_basis   = 0
        i_primit  = 0
        i_elements = 0
        primit_arr = np.zeros((3,3))

        for line in fr:
            
            if latt_param_flag == True:
                latt_param = float(line.split()[0]) * 0.529177    # bohr radius to Angstrom
                latt_param_flag = False
                
            if number_flag == True:
                n_atoms = int(line.split()[0])
                number_flag = False

            if primitive_vectors_flag == True:
                primit_arr[i_primit,:] = np.array(line.split())
                i_primit += 1
                if i_primit > 2:
                    primitive_vectors_flag = False
                    primit_arr *= latt_param

            if basis_vectors_flag == True:
                basis_arr[i_basis, 1:4] = np.array(line.split())[2:5]
                i_basis += 1
                if i_basis >= n_atoms:
                    basis_vectors_flag = False

            if elements_flag == True:
                elements_arr.append(line.split()[2])
                i_elements += 1
                if i_elements >= n_atoms:
                    elements_flag = False

            if "lattice parameter A" in line:
                latt_param_flag = True

            if "primitive vectors" in line:
                primitive_vectors_flag = True

            if "number of sites NQ" in line:
                number_flag = True

            if "basis vectors" in line:
                basis_vectors_flag = True
                basis_arr = np.zeros((n_atoms, 5))

            if "IQAT (sites occupied)" in line:
                elements_flag = True
                elements_arr = []

    basis_arr[:,0] = np.array(range(n_atoms))
    basis_arr[:,4] = np.array(range(n_atoms))

    print(".sys file parsed")
    return n_atoms, primit_arr, basis_arr, elements_arr

            

def sprkkr_to_vampire_ucf(path, system_name, include_dmi=True):
    """Convert exchange and DMI from SPR-KKR into VAMPIRE UCF file."""

    meV = 1.602e-22  #J

    # SPR-KKR J and DMI .dat files in
    f_exchange_in = f"{path}/{system_name}_JXC_XCPLTEN_Jij.dat"
    f_dmi_in      = f"{path}/{system_name}_JXC_XCPLTEN_Dij.dat"

    # UCF file out
    fout = f"{path}/vampire.UCF"

    # get the structure data from .sys file
    n_atoms, primit_arr, basis_arr, elements_arr = get_structure_from_sys_sprkkr(path, system_name)

    # write the .mat file
    mag_moments = get_mag_moments(path, system_name, n_atoms, f_electrons=f_electrons)
    write_mat_file(path, system_name, mag_moments, elements_arr)

    print(elements_arr)
    #with open(f_exchange_in, 'r') as fr:
    df = pd.read_csv(f_exchange_in, skiprows=9+n_atoms, delim_whitespace=True, names= ['IT', 'IQ', 'JT', 'JQ', 'N1', 'N2', 'N3', 'DRX', 'DRY', 'DRZ', 'DR', 'J_xx', 'J_yy', 'J_xy', 'J_yx'] )
    df.drop(['IQ','JQ', 'J_yy', 'J_xy', 'J_yx', 'DRX', 'DRY', 'DRZ', 'DR'], axis=1, inplace=True)

    # indexing from 0 instead of 1 for VAMPIRE
    df['IT'] = df['IT']-1
    df['JT'] = df['JT']-1
    # make a tensor with zeros
    df.rename(columns={'J_xx':'J_11'}, inplace=True)
    df[['J_12', 'J_13', 'J_21', 'J_22', 'J_23', 'J_31', 'J_32', 'J_33']] = 0

    # ----- SET THE HEISENBERG MATRIX VALUES -----
    df['J_11'] = df['J_11'] * meV
    df['J_22'] = df['J_11']
    df['J_33'] = df['J_11']
    print("new matrix columns created")

    if include_dmi == True:
            df_dmi = pd.read_csv(f_dmi_in, skiprows=12+n_atoms, delim_whitespace=True, names= ['IT', 'IQ', 'JT', 'JQ', 'N1', 'N2', 'N3', 'DRX', 'DRY', 'DRZ', 'DR', 'DX', 'DY', 'DZ'] )
            print('check4')
            # see e.g. Coey for definition of DMI -> DMI tensor components ( (0, -Dz, Dy), (Dz, 0, -Dx), (-Dy, Dx, 0) )
            df['J_12'] += -df_dmi['DZ'] * meV
            df['J_13'] += df_dmi['DY'] * meV
            df['J_21'] += df_dmi['DZ'] * meV
            df['J_23'] += -df_dmi['DX'] * meV
            df['J_31'] += -df_dmi['DY'] * meV
            df['J_32'] += df_dmi['DX'] * meV
            print("DMI included")

    # --------------------------------------------

    df.sort_values(by=['IT', 'JT', 'N1', 'N2', 'N3'], inplace=True)
    df.reset_index(inplace=True)
    df.drop(labels='index', axis=1, inplace=True)
    print("datafield cleaned and sorted")

    # save to UCF file
    with open(fout, 'w') as fwrite:
        fwrite.write(f"# Unit cell size (Angstrom):\n1 1 1\n# Unit cell lattice vectors:\n")
        print('check')
        np.savetxt(fwrite, primit_arr, fmt='%.6f', delimiter=' ')
        fwrite.write(f"# Atoms\n{n_atoms} {n_atoms}\n")
        np.savetxt(fwrite, basis_arr, fmt='%d %.6f %.6f %.6f %d', delimiter=' ')
        fwrite.write(f"# Interactions\n{df.shape[0]} tensorial\n")
        df.to_csv(fwrite, mode='w', header=False, sep=" ", line_terminator='\n')
    print(".UCF file written")


    if crop_threshold > 0:
        ucf_crop(path+'/vampire.UCF', n_atoms+10, crop_threshold, save_file=True)


def main():
    sprkkr_to_vampire_ucf(path, system_name, include_dmi)


if __name__ == "__main__":
    main()
