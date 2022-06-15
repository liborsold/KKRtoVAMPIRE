"""Convert the SPR-KKR exchange interactions to vampire.
    Needed are the files: seedname.sys, seedname_SCF.out, seedname_JXC_XCPLTEN_Jij.dat and seedname_JXC_XCPLTEN_Dij.dat.
    Specify if DMI should be included, if f-electrons are present (for parsing) and if the interactions should be croppped."""

import pandas as pd
import numpy as np
from os.path import exists
import subprocess

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '../ucf_crop_small_values_vampire')
from ucf_crop_small_values_vampire import ucf_crop


# ========================== USER INPUT =============================

# SPR-KKR calculation directory
path = "H:/test_exchange_Fe_Co_Ni/Ni/a_lit" #"H:/2D/Cr2Te3/bulk/sprkkr_imitating_experimental/test_structure"   #  "H:/SPR-KKR/bcc_Fe" # "H:/2D/Cr2Te3/bulk/sprkkr/PBE/NKTAB_1000" #  
# the seed name in SPR-KKR (the name of the system)
system_name =  "Ni" #"POSCAR" # "Fe"   #

include_dmi = False
include_anisotropy = False

crop_threshold = 0 #0.1    # meV; set value > 0 to perform the crop post-processing (via the ucf_crop_small_values_vampire.py script; creates a new file)

# ===================================================================

def get_torques(path, system_name, n_atoms):
    """Parse torques from _TORQUE.out file."""
    torques = []
    lines_to_skip = 0

    path_in = f"{path}/{system_name}_TORQUE.out"

    torque_flag, read_flag = False, False
    skipped_lines = 0

    with open(path_in, 'r') as fr:
        for line in fr:
            if torque_flag == True:
                skipped_lines += 1

            if "torque (Ryd)" in line:
                torque_flag = True

            if n_atoms >= skipped_lines > 0:
                torques.append(line.split()[3])

    print("torques obtained")
    return torques



def get_latt_params_for_input_file(path, system_name):
    """Write part of the 'input' file for VAMPIRE containing the lattice parametres."""
    fread = f"{path}/{system_name}.sys"
    fwrite = f"{path}/input_latt-params"
    with open(fread, 'r') as fr:
        latt_params_flag = False
        for line in fr:
            if latt_params_flag == True:
                latt_params = np.array(line.split(), dtype=float) * 0.529177  # a.u. to Angstrom
                break

            if "lattice parameters  a b c" in line:
                latt_params_flag = True
        
    with open(fwrite, 'w') as fw:
        fw.write(f"dimensions:unit-cell-size-x = {latt_params[0]:.12f}\n")
        fw.write(f"dimensions:unit-cell-size-y = {latt_params[1]:.12f}\n")
        fw.write(f"dimensions:unit-cell-size-z = {latt_params[2]:.12f}\n")
    
    return latt_params


def write_mat_file(path, system_name, mag_moments, elements, torques):
    """Write the .mat file for vampire, given the magnetic moments and types of elements."""
    if len(mag_moments) != len(elements):
        raise ValueError("elements and mag_moments arrays need to be of the same length!")

    path_in = f"{path}/vampire.mat"

    with open(path_in, 'w') as fw:
        fw.write(f"material:num-materials = {len(elements)}\n")
        for i in range(len(elements)):
            magmom = float(mag_moments[i])
            torque = float(torques[i]) * 13.6 * 1.602e-19
            fw.write(
                      (  "#---------------------------------------------------\n"
                         f"# Material {i+1} \n"
                         "#---------------------------------------------------\n"
                         f"material[{i+1}]:material-name={elements[i]}\n"
                         f"material[{i+1}]:damping-constant=1.0\n"
                         f"material[{i+1}]:atomic-spin-moment={abs(magmom):.8f} !muB\n"
                         f"material[{i+1}]:uniaxial-anisotropy-constant={torque:.8e}\n"
                         f"material[{i+1}]:material-element={elements[i]}\n"                        
                         f"material[{i+1}]:initial-spin-direction = 0.0,0.0,{np.sign(magmom):.1f}\n"
                         f"material[{i+1}]:uniaxial-anisotropy-direction = 0.0 , 0.0, 1.0\n"
                         "#---------------------------------------------------\n"
                       )
            )
    print(".mat file written")


def get_mag_moments(path, system_name, n_atoms):
    """Parse magnetic moments from _SCF.out file."""
    mag_moments = []

    path_in = f"{path}/{system_name}_SCF.out"

    # # use bash to reduce the _SCF.out file
    # # ** NOT WORKING properly **
    # path_in_reduced = path_in + "_magmoms"
    # if not exists(path_in_reduced):
    #     bashCommand = f"grep -A5 m_spin {path_in} | tail -{5*n_atoms}"
    #     process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    #     output, error = process.communicate()
    #     with open(path_in_reduced, 'w') as fw:
    #         fw.write(str(output))
    #     print(error)

    IT_flag = False

    with open(path_in, 'r') as fr:
        for line in fr:
            if IT_flag == True:
                if line.split()[0] == "sum":
                    mag_moments.append(line.split()[4])
                    IT_flag = False

            if "IT=" in line:
                IT_flag = True


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
    return n_atoms, primit_arr, basis_arr, elements_arr, latt_param

            

def sprkkr_to_vampire_ucf(path, system_name, include_dmi=True, include_anisotropy=True):
    """Convert exchange and DMI from SPR-KKR into VAMPIRE UCF file."""

    meV = 1.602e-22  #J

    # SPR-KKR J and DMI .dat files in
    f_exchange_in = f"{path}/{system_name}_JXC_XCPLTEN_Jij.dat"
    f_dmi_in      = f"{path}/{system_name}_JXC_XCPLTEN_Dij.dat"

    # UCF file out
    fout = f"{path}/vampire.UCF"

    # get the structure data from .sys file
    n_atoms, primit_arr, basis_arr, elements_arr, latt_param = get_structure_from_sys_sprkkr(path, system_name)

    # write the .mat file
    mag_moments = get_mag_moments(path, system_name, n_atoms)
    torques = get_torques(path, system_name, n_atoms) if include_anisotropy else np.zeros((n_atoms,))
    print(mag_moments)
    write_mat_file(path, system_name, mag_moments, elements_arr, torques)

    # write the input file with the correct lattice parameters
    latt_params = get_latt_params_for_input_file(path, system_name)

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
    print("all matrix columns created")

    if include_dmi == True:
            df_dmi = pd.read_csv(f_dmi_in, skiprows=12+n_atoms, delim_whitespace=True, names= ['IT', 'IQ', 'JT', 'JQ', 'N1', 'N2', 'N3', 'DRX', 'DRY', 'DRZ', 'DR', 'DX', 'DY', 'DZ'] )
            # see e.g. Coey for definition of DMI -> DMI tensor components ( (0, -Dz, Dy), (Dz, 0, -Dx), (-Dy, Dx, 0) )
            df['J_12'] += -df_dmi['DZ'] * meV
            df['J_13'] += df_dmi['DY'] * meV
            df['J_21'] += df_dmi['DZ'] * meV
            df['J_23'] += -df_dmi['DX'] * meV
            df['J_31'] += -df_dmi['DY'] * meV
            df['J_32'] += df_dmi['DX'] * meV
            print("DMI included")
    else:
        print("DMI will *NOT* be included")

    # --------------------------------------------

    df.sort_values(by=['IT', 'JT', 'N1', 'N2', 'N3'], inplace=True)
    df.reset_index(inplace=True)
    df.drop(labels='index', axis=1, inplace=True)
    print("datafield cleaned and sorted")

    # save to UCF file
    primit_arr_reduced = np.array( [primit_arr[0,:]/latt_params[0], primit_arr[1,:]/latt_params[1], primit_arr[2,:]/latt_params[2]] )
    basis_arr_reduced = np.array( [basis_arr[:,0], basis_arr[:,1]*latt_param/latt_params[0], basis_arr[:,2]*latt_param/latt_params[1], basis_arr[:,3]*latt_param/latt_params[2], basis_arr[:,4]] ).T
    
    # coerce values  between  0 and 1 
    n_basis, m_basis = basis_arr_reduced.shape
    basis_arr_reduced = np.array( [[ basis_arr_reduced[i,j] if basis_arr_reduced[i,j] >= 0 else 1+basis_arr_reduced[i,j] for j in range(m_basis)] for i in range(n_basis) ]  )

    with open(fout, 'w') as fwrite:
        fwrite.write(f"# Unit cell size (Angstrom):\n")
        np.savetxt(fwrite, [latt_params], fmt='%.8f', delimiter=' ')
        fwrite.write("# Unit cell lattice vectors:\n")
        np.savetxt(fwrite, primit_arr_reduced, fmt='%.6f', delimiter=' ')
        fwrite.write(f"# Atoms\n{n_atoms} {n_atoms}\n")
        np.savetxt(fwrite, basis_arr_reduced, fmt='%d %.6f %.6f %.6f %d', delimiter=' ')
        fwrite.write(f"# Interactions\n{df.shape[0]} tensorial\n")
        df.to_csv(fwrite, mode='w', header=False, sep=" ", line_terminator='\n')
    print(".UCF file written")


    if crop_threshold > 0:
        ucf_crop(path, 'vampire.UCF', crop_threshold, save_file=True)


def main():
    sprkkr_to_vampire_ucf(path, system_name, include_dmi, include_anisotropy)


if __name__ == "__main__":
    main()
