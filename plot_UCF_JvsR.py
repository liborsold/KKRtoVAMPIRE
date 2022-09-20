"""Plot UCF J vs R in the same way as TB2J does it, i.e., distance-based."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

plot_cumulative_J = True

labels = ['IID', 'i', 'j', 'dx', 'dy', 'dz', 'Jxx', 'Jxy', 'Jxz', 'Jyx', 'Jyy', 'Jyz', 'Jzx', 'Jzy', 'Jzz']

def get_n_atoms_UCF(fname):
    """Get number of atoms from UCF file."""
    with open(fname, 'r') as fr:
        for i, line in enumerate(fr):
            if i == 7:
                n_atoms = int(line.split()[0])
                return n_atoms


def get_rows_to_skip_UCF(fread):
    """Find how many lines to skip in UCF file."""
    n_atoms = get_n_atoms_UCF(fread)
    return n_atoms + 10


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
    atom_coordinates = np.array(atom_coordinates)
    # transform atom_coordinates from fractional to cartesian
    # wrong! : atom_coordinates_cartes = np.array(atom_coordinates) @ uc_vectors
    for i in range(3):
        atom_coordinates[:,i] *= prefactors[i]
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


def plot_JvsR_UCF(fread, fname_out):
    rows_to_skip = get_rows_to_skip_UCF(fread)

    uc_vectors, atom_coordinates = structure_data_UCF(fread)
    with open(fread, 'r') as fname:
        df = pd.read_csv(fname, skiprows=rows_to_skip, delim_whitespace=True, names= labels)
    df['dr'] = distance_column(df, uc_vectors, atom_coordinates)

    # order interactions by distance
    df.sort_values(by=['dr'], inplace=True)

    df['dr'] /= uc_vectors[0,0]
    df['Jxx_meV'] = df['Jxx']/1.602e-22
    df[['dr', 'Jxx_meV']].to_csv("Jxx_vs_dr_from_UCF.txt", sep='\t', float_format='%.10f')

    x = list(np.array(df['dr'], dtype=np.float32))
    y = list(np.array(df['Jxx'], dtype=np.float32) * 1000 / 1.602e-19)

    plt.plot(x, y, 'o')
    plt.title(r"$J_\mathrm{xx}$ vs. distance from vampire.UCF file" + f"\n{os.getcwd()}\n{fname_out}", fontsize=8)
    plt.ylabel(r"$J_\mathrm{xx}$ (meV)")
    plt.xlabel(r'$r$/$a$')
    # plt.xlim([-2, 50])
    # plt.ylim([-17, 17])
    plt.tight_layout()
    # plt.show()
    plt.savefig(f"{fname_out}_JvsR.png", dpi=400)
    plt.close()

    if plot_cumulative_J == True:
        y_summed = list(np.array(df['Jxx'].cumsum(), dtype=np.float32) * 1000 / 1.602e-19)
        plt.plot(x, y_summed, 'ok-', linewidth=2)
        plt.title(r"Cumulative $J_\mathrm{xx}$ vs. distance from vampire.UCF file" + f"\n{os.getcwd()}\n{fname_out}", fontsize=8)
        plt.ylabel(r"Cumulative $J_\mathrm{xx}$ (meV)")
        plt.xlabel(r'$r$/$a$')
        #plt.xlim([-2, 50])
        #plt.ylim([-50, 1200])
        plt.tight_layout()
        # plt.show()
        plt.savefig(f"{fname_out}_cumulative_JvsR.png", dpi=400)
        plt.close()


def main():
    files = glob.glob('vampire.UCF*', recursive = False)
    for file in files:
        if not "png" in file:
            plot_JvsR_UCF(file, fname_out=file)


if __name__ == "__main__":
    main()
    