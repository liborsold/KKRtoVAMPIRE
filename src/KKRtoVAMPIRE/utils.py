"""Convert the SPR-KKR exchange interactions to vampire.
    Needed are the files: seedname.pot_new, seedname_SCF.out, seedname_JXC_XCPLTEN_Jij.dat and seedname_JXC_XCPLTEN_Dij.dat.
    Specify if DMI and anisotropy should be included and if the interactions should be croppped.
    
    !!!! tested on .pot file with  >> FORMAT     9 (18.01.2019) <<
    
    """

from multiprocessing.sharedctypes import Value
import pandas as pd
import numpy as np
from os.path import exists
import subprocess
from copy import copy

class Converter():
    def __init__(self, path, system_name, include_dmi=False, include_anisotropy=True, crop_thresholds=[0]):
        self.path = path
        self.system_name = system_name
        self.include_dmi = include_dmi
        self.include_anisotropy = include_anisotropy
        self.crop_thresholds = crop_thresholds


def ucf_crop(path='.', file_name='vampire.UCF', threshold=0.00, save_file=True, \
             convert_tensorial_to_isotropic=True, thresholds_in_fraction_of_total=True):
    """Drop values from fread which

            -- if thresholds_in_fraction_of_total == True --
        (1) constitute the smallest sum of the <threshold> fraction of the total (in absolute value), or

            -- if thresholds_in_fraction_of_total == False --
        (2) are lower than <threshold>, given in meV."""

    def get_rows_to_skip_in_UCF(path, file_name):
        """Find how many lines wo skip in UCF file; needed for ucf_crop()."""
        fread = f"{path}/{file_name}"
        rows_to_skip = 0
        with open(fread, 'r') as fr:
            for line in fr:
                rows_to_skip += 1
                if "Interactions" in line:
                    return rows_to_skip + 1

    thresh_unit = '' if thresholds_in_fraction_of_total == True else 'meV'
    threshold_J = threshold * 1.602e-22    # convert meV to J, since J is used in .ucf
    print(f"----------- threshold {threshold:.3f} {thresh_unit} -----------")

    rows_to_skip = get_rows_to_skip_in_UCF(path, file_name)

    fread = f"{path}/{file_name}"
    with open(Path(fread), 'r') as fr:
        df = pd.read_csv(fread, skiprows=rows_to_skip, delim_whitespace=True, names= ['IID', 'i', 'j', 'dx', 'dy', 'dz', 'Jxx', 'Jxy', 'Jxz', 'Jyx', 'Jyy', 'Jyz', 'Jzx', 'Jzy', 'Jzz'] )
        n_orig = df.shape[0]

        # sort data by the decreasing value of sum of abs() of all columns in each row
          # abs sum of all columns
        df['abs_sum_all'] = df[['Jxx', 'Jxy', 'Jxz', 'Jyx', 'Jyy', 'Jyz', 'Jzx', 'Jzy', 'Jzz']].abs().sum(axis=1)
        df.sort_values('abs_sum_all', inplace=True, ascending=False)

        if thresholds_in_fraction_of_total == True:
            # crop all interactions that contribute in the lowest <threshold> fraction of the sum of abs() interactions
               # cumulative abs sum of all columns
            df['cum_abs_sum_all'] = df['abs_sum_all'].cumsum()
            sum_total = df['cum_abs_sum_all'].iloc[-1]
            df = df[ df['cum_abs_sum_all'] <= (1-threshold)*sum_total ]
            df.drop('cum_abs_sum_all', 1, inplace=True)
            df.drop('abs_sum_all', 1, inplace=True)

        else: 
            # crop all interactions which are smaller than <threshold> (in meV)
            df.drop('abs_sum_all', 1, inplace=True)
            df = df[ df[['Jxx', 'Jxy', 'Jxz', 'Jyx', 'Jyy', 'Jyz', 'Jzx', 'Jzy', 'Jzz']].abs().max(axis=1) >= threshold_J ]

        n_after_crop = df.shape[0]

        fwrite = fread + f"_cropped_{n_after_crop}_{threshold:.3f}{thresh_unit}"

        with open(fwrite, 'w') as fw:
            i = 0
            for line in fr:
                if i >= rows_to_skip-1:
                    break 
                fw.write(line)
                i += 1
            mode = 'isotropic' if convert_tensorial_to_isotropic == True else 'tensorial'
            fw.write(f"{n_after_crop} {mode}\n")
        df.reset_index(inplace=True)
        df.drop(labels=['IID', 'index'], axis=1, inplace=True)
        if convert_tensorial_to_isotropic == True:
            df.drop(labels=['Jxy', 'Jxz', 'Jyx', 'Jyy', 'Jyz', 'Jzx', 'Jzy', 'Jzz'], axis=1, inplace=True)
        df.to_csv(fwrite, mode='a', header=False, sep=" ")
        print("cropped .UCF file written")
        return [n_after_crop, n_orig]

def get_torques(path, system_name, n_atoms, types_arr):
    """Parse torques from _TORQUE.out file."""
    torques_atoms = []
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
                torques_atoms.append(line.split()[3])

    print("torques obtained")

    n_atoms = len(types_arr)
    n_types = max(types_arr)
    torques_types = np.zeros((n_types,))
    for i in range(n_atoms):
        torques_types[types_arr[i]-1] = torques_atoms[i]
    return torques_types

def write_latt_params_for_input_file(path, latt_params):
    """Write part of the 'input' file for VAMPIRE containing the lattice parametres."""
    fwrite = f"{path}/input_latt-params"
    with open(fwrite, 'w') as fw:
        fw.write(f"dimensions:unit-cell-size-x = {latt_params[0]:.12f}\n")
        fw.write(f"dimensions:unit-cell-size-y = {latt_params[1]:.12f}\n")
        fw.write(f"dimensions:unit-cell-size-z = {latt_params[2]:.12f}\n")
    return 0

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
                         f"material[{i+1}]:material-element={elements[i].split('_')[0]}\n"                        
                         f"material[{i+1}]:initial-spin-direction = 0.0,0.0,{np.sign(magmom):.1f}\n"
                         f"material[{i+1}]:uniaxial-anisotropy-direction = 0.0 , 0.0, 1.0\n"
                         "#---------------------------------------------------\n"
                       )
            )
    print(".mat file written")

def get_mag_moments(path, system_name, n_atoms, n_types, types_arr):
    """Parse magnetic moments from _SCF.out file."""
    mag_moments_types = []

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
        for i, line in enumerate(fr):
            if IT_flag == True:
                try:
                    if line.split()[0] == "sum":
                        mag_moments_types.append(line.split()[4])
                        IT_flag = False
                except IndexError:
                    print(f"IndexError on line {i} of file!")

            if "IT=" in line:
                IT_flag = True

    mag_moments_types = mag_moments_types[-n_types:]
    mag_moments_elements = [mag_moments_types[types_arr[i]-1] for i in range(n_atoms)]
    print(f"mag_moments obtained: {mag_moments_elements}")
    return mag_moments_elements, mag_moments_types

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

def get_structure_from_pot_sprkkr(path, system_name):
    """From the SPR-KKR's .pot file get the unit cell, number of atoms and their basis vectors."""
    path_in = f"{path}/{system_name}_SCF.pot_new"

    with open(path_in, 'r') as fr:
        basis_vectors_flag, elements_flag, types_flag = False, False, False
        delay_bas = True
        i_basis   = 0
        i_types = 0
        i_elements = 0
        primit_arr = np.zeros((3,3))

        for line in fr:

            if basis_vectors_flag == True:
                if delay_bas == False:
                    basis_arr[i_basis, 1:4] = np.array(line.split())[1:4]
                    i_basis += 1
                    if i_basis >= n_atoms:
                        basis_vectors_flag = False
                else:
                    delay_bas = False

            if types_flag == True:
                types_arr.append(int(line.split()[4]))
                i_types += 1
                if i_types >= n_atoms:
                    types_flag = False
                    
            if elements_flag == True:
                elements_arr.append(line.split()[1])
                i_elements += 1
                if i_elements >= n_types:
                    elements_flag = False

            if "ALAT" in line:
                latt_param = float(line.split()[1]) * 0.529177    # bohr radius to Angstrom

            if "A(1)" in line:
                lsplit = line.split()
                primit_arr[0,:] = np.array([float(lsplit[i]) for i in range(1,4)]) * latt_param  

            if "A(2)" in line:
                lsplit = line.split()
                primit_arr[1,:] = np.array([float(lsplit[i]) for i in range(1,4)]) * latt_param  

            if "A(3)" in line:
                lsplit = line.split()
                primit_arr[2,:] = np.array([float(lsplit[i]) for i in range(1,4)]) * latt_param  

            if "NQ" in line:
                n_atoms = int(line.split()[1])

            if line.startswith("NT"):
                n_types = int(line.split()[1])

            if "BASSCALE" in line:
                basis_vectors_flag = True
                delay_bas = True
                basis_arr = np.zeros((n_atoms, 5))

            if "ITOQ" in line:
                types_flag = True
                types_arr = []

            if "NSEMCORSHLT" in line:
                elements_flag = True
                elements_arr = []
                

    basis_arr[:,0] = np.array(range(n_atoms))
    basis_arr[:,4] = np.array(types_arr)-1

    types_names = copy(elements_arr)
    atoms_names = [elements_arr[types_arr[i]-1] for i in range(n_atoms)]

    print(".pot file parsed")
    return n_atoms, primit_arr, basis_arr, atoms_names, types_names, latt_param, types_arr, n_types


def structure_data_UCF(fname):
    """Return 'uc_vectors' and 'atom_coordinates' for futher distance calculations."""
    uc_vectors = np.zeros((3,3))
    atom_coordinates = []
    with open(fname, 'r') as fr:
        for i, line in enumerate(fr):
            if "Interactions" in line:
                break
            if i == 1:
                prefactors = [float(number) for number in line.split()]
            if i == 3:
                uc_vectors[0,:] = np.array([float(number) for number in line.split()]) * prefactors[0]
            if i == 4:
                uc_vectors[1,:] = np.array([float(number) for number in line.split()]) * prefactors[1]
            if i == 5:
                uc_vectors[2,:] = np.array([float(number) for number in line.split()]) * prefactors[2]
            if i >= 8:
                l_split = line.split()
                atom_coordinates.append([float(l_split[i]) for i in range(1,4)])
    
    # transform atom_coordinates from fractional to cartesian
    atom_coordinates = np.array(atom_coordinates) @ uc_vectors
    print(atom_coordinates)
    return uc_vectors, atom_coordinates


def distance_column(df, uc_vectors, atom_coordinates):
    """Take pandas df with UCF file data and calculate distance for each interaction.
    - df is the pandas datafield array
    - uc_vectors is a vector of the unit cell vectors, i.e., uc_vector[0] = ucx, uc_vector[1] = ucy, etc.
    - atom_coordinates is the fractional coordinates of all atoms (in terms of the uc_vectors)"""
    r1 = atom_coordinates[df['i']]
    r2 = atom_coordinates[df['j']] + uc_vectors[0,:]*df[['dx']].to_numpy() + uc_vectors[1,:]*df[['dy']].to_numpy() + uc_vectors[2,:]*df[['dz']].to_numpy()
    return np.sqrt(np.sum(np.power((r1-r2),2), axis=1))


def coerce_overflowing_atoms_to_unit_cell(basis_arr_reduced, df, primit_arr, latt_params):
    """Coerce fractional coordinates of atoms between 0.0 and 1.0.
         When fractional coordinate 
          - more than 1.0: subtract 1.0, but also change the positions of the unit cells in the UCF data (e.g. from 0 0 0 to 0 0 1)
          - less than 0.0: add 1.0, but also change the positions of the unit cells in the UCF data (e.g. from 0 0 0 to 0 0 -1)"""

    def shift_unit_cell_positions(df, atom_number, coordinate, shift):
        """In df, in all interactions involving 'atom_number', 
                            ADD      'shift' to    'coordinate' if atom is first in the interaction,
                            SUBTRACT 'shift' from  'coordinate' if atoms is second in the interaction. 
           
           ...                 
           (Given the parsed 'Jij.dat' in a pandas datafield 'df').
           (Column names: ['IT', 'IQ', 'JT', 'JQ', 'N1', 'N2', 'N3', 'DRX', 'DRY', 'DRZ', 'DR', 'J_xx', 'J_yy', 'J_xy', 'J_yx'] )"""
        coordinate_name = {0:'N1', 1:'N2', 2:'N3'}
        column = coordinate_name[coordinate]

        df.loc[df['IQ']==atom_number, [column]] += shift
        df.loc[df['JQ']==atom_number, [column]] -= shift
        return df

    def what_combination_of_latt_vectors_coerces_to_basic_unit_cell(basis_arr_reduced, latt_params, primit_arr):
        atom_coord_cart = basis_arr_reduced[1:4] * latt_params

        # if already coerced, don't do anything
        if 0 <= atom_coord_cart[0] < latt_params[0] and 0 <= atom_coord_cart[1] < latt_params[1] and 0 <= atom_coord_cart[2] < latt_params[2]:
            return basis_arr_reduced, 0, 0, 0

        # else look what transformation will coerce it
        for du1 in range(-1,2):
            for du2 in range(-1,2):
                for du3 in range(-1,2):
                    new_atom_coord_cart = atom_coord_cart + du1 * primit_arr[0,:] + du2 * primit_arr[1,:] + du3 * primit_arr[2,:]
                    if 0 <= new_atom_coord_cart[0] < latt_params[0] and 0 <= new_atom_coord_cart[1] < latt_params[1] and 0 <= new_atom_coord_cart[2] < latt_params[2]:
                        basis_arr_reduced[1:4] = (new_atom_coord_cart[0]/latt_params[0], new_atom_coord_cart[1]/latt_params[1], new_atom_coord_cart[2]/latt_params[2])
                        return (basis_arr_reduced, du1, du2, du3)

    n_basis = basis_arr_reduced.shape[0]

    for i in range(n_basis):
        basis_arr_reduced[i,:], du1, du2, du3 = what_combination_of_latt_vectors_coerces_to_basic_unit_cell(basis_arr_reduced[i,:], latt_params, primit_arr)
        df = shift_unit_cell_positions(df, atom_number=i, coordinate=0, shift=du1)
        df = shift_unit_cell_positions(df, atom_number=i, coordinate=1, shift=du2)
        df = shift_unit_cell_positions(df, atom_number=i, coordinate=2, shift=du3)

        # for j in range(1,4):
        #     # j: 1 - x, 2 - y, 3 - z
        #     if basis_arr_reduced[i,j] < 0.0:
        #         basis_arr_reduced[i,j] += 1.0
        #         df = shift_unit_cell_positions(df, atom_number=i, coordinate=j-1, shift=1)
            
        #     if basis_arr_reduced[i,j] >= 1.0 or f"{basis_arr_reduced[i,j]:.6f}" == "1.000000":
        #         basis_arr_reduced[i,j] -= 1.0
        #         df = shift_unit_cell_positions(df, atom_number=i, coordinate=j-1, shift=-1)

    basis_arr_reduced = np.array( basis_arr_reduced )

    return basis_arr_reduced, df

def sprkkr_to_vampire_ucf(path, system_name, crop_thresholds, include_dmi=True, include_anisotropy=True):
    """Convert exchange and DMI from SPR-KKR into VAMPIRE UCF file."""

    meV = 1.602e-22  #J

    # SPR-KKR J and DMI .dat files in
    f_exchange_in = f"{path}/{system_name}_JXC_XCPLTEN_Jij.dat"
    f_dmi_in      = f"{path}/{system_name}_JXC_XCPLTEN_Dij.dat"

    # UCF file out
    fout = f"{path}/vampire.UCF"

    # get the structure data from .sys file
    n_atoms, primit_arr, basis_arr, atoms_names, types_names, latt_param, types_arr, n_types = get_structure_from_pot_sprkkr(path, system_name) # get_structure_from_sys_sprkkr(path, system_name) <- frmo sys file outdated
    print(n_atoms)
    print(primit_arr)
    print(basis_arr)
    print(atoms_names)
    print(types_names)
    print(latt_param)

    # write the .mat file
    mag_moments_atoms, mag_moments_types = get_mag_moments(path, system_name, n_atoms, n_types, types_arr)

    # include anisotropy
    if include_anisotropy == True:
        torques = get_torques(path, system_name, n_atoms, types_arr)
    else:
        torques = np.zeros((n_atoms,))
        print("anisotropy will *NOT* be included")

    print(mag_moments_types)
    write_mat_file(path, system_name, mag_moments_types, types_names, torques)

    # calculate the length of the effective lattice parameters
    latt_params = np.array( [np.sqrt(np.sum(np.power(primit_arr[i,:], 2))) for i in range(3)] )
    print(f"latt_params: {latt_params}")

    # write the input file with the correct lattice parameters
    write_latt_params_for_input_file(path, latt_params)

    print(atoms_names)
    #with open(f_exchange_in, 'r') as fr:
    df = pd.read_csv(f_exchange_in, skiprows=9+n_atoms, delim_whitespace=True, names= ['IT', 'IQ', 'JT', 'JQ', 'N1', 'N2', 'N3', 'DRX', 'DRY', 'DRZ', 'DR', 'J_xx', 'J_yy', 'J_xy', 'J_yx'] )
    df.drop(['IT','JT', 'J_yy', 'J_xy', 'J_yx', 'DRX', 'DRY', 'DRZ', 'DR'], axis=1, inplace=True)

    # indexing from 0 instead of 1 for VAMPIRE
    df['IQ'] = df['IQ']-1
    df['JQ'] = df['JQ']-1
    # make a tensor with zeros
    df.rename(columns={'J_xx':'J_11'}, inplace=True)
    df[['J_12', 'J_13', 'J_21', 'J_22', 'J_23', 'J_31', 'J_32', 'J_33']] = 0

    # ----- SET THE HEISENBERG MATRIX VALUES -----
    df['J_11'] = 2 * df['J_11'] * meV
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

    # df.sort_values(by=['IT', 'JT', 'N1', 'N2', 'N3'], inplace=True)
    # df.reset_index(inplace=True)
    # df.drop(labels='index', axis=1, inplace=True)
    # print("datafield cleaned and sorted")

    # COMPLETENESS CHECK
    # If there is an interaction between two atoms of certain types 'a' and 'b', check if the same interaction is listed for all atoms of types 'a' and 'b'
    

    # convert to fractional
    primit_arr_reduced = np.array( [primit_arr[0,:]/latt_params[0], primit_arr[1,:]/latt_params[1], primit_arr[2,:]/latt_params[2]] )
    basis_arr_reduced = np.array( [basis_arr[:,0], basis_arr[:,1]*latt_param/latt_params[0], basis_arr[:,2]*latt_param/latt_params[1], basis_arr[:,3]*latt_param/latt_params[2], basis_arr[:,4]] ).T

    #basis_arr_reduced, df = coerce_overflowing_atoms_to_unit_cell(basis_arr_reduced, df, primit_arr, latt_params)

    # save to .UCF file
    with open(fout, 'w') as fwrite:
        fwrite.write(f"# Unit cell size (Angstrom):\n")
        np.savetxt(fwrite, np.array([[1], [1], [1]]).T, fmt='%.0f', delimiter=' ')
        fwrite.write("# Unit cell lattice vectors:\n")
        np.savetxt(fwrite, primit_arr, fmt='%.12f', delimiter=' ')
        fwrite.write(f"# Atoms\n{n_atoms} {n_atoms}\n")
        np.savetxt(fwrite, basis_arr_reduced, fmt='%d %.12f %.12f %.12f %d', delimiter=' ')
        fwrite.write(f"# Interactions; J values from SPR-KKR multiplied by 2; DMI not yet \n{df.shape[0]} tensorial\n")
        df.to_csv(fwrite, mode='w', header=False, sep=" ", line_terminator='\n')
    print(".UCF file written")

    for crop_threshold in crop_thresholds:
        if crop_threshold >= 0:
            ucf_crop(path, 'vampire.UCF', crop_threshold, save_file=True)


def main():
    sprkkr_to_vampire_ucf(path='.', system_name='POSCAR', crop_thresholds=[0], \
                          include_dmi=False, include_anisotropy=True)


if __name__ == "__main__":
    main()
