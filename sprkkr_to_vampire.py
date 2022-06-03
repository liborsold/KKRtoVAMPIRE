"""Convert the SPR-KKR exchange interactions to vampire."""

import pandas as pd
import numpy as np




# ========================== USER INPUT =============================

# SPR-KKR calculation directory
path = "H:/2D/Cr2Te3/bulk/transport/PBE/NKTAB_1000"

# the seed name in SPR-KKR (the name of the system)
system_name = "Cr2Te3"
n_atoms     = 20

# ===================================================================




def sprkkr_to_vampire_ucf(path, system_name, n_atoms):
    """Convert exchange and DMI from SPR-KKR into VAMPIRE UCF file."""

    eV = 1.602e-19  #J

    # SPR-KKR J and DMI .dat files in
    f_exchange_in = f"{path}/{system_name}_JXC_XCPLTEN_Jij.dat"
    f_dmi_in      = f"{path}/{system_name}_JXC_XCPLTEN_Dij.dat"

    # UCF file out
    fout = f"vampire.UCF_from_{system_name}"

    with open(f_exchange_in, 'r') as fr:
        df = pd.read_csv(f_exchange_in, skiprows=9+n_atoms, delim_whitespace=True, names= ['IT', 'IQ', 'JT', 'JQ', 'N1', 'N2', 'N3', 'DRX', 'DRY', 'DRZ', 'DR', 'J_xx', 'J_yy', 'J_xy', 'J_yx'] )
        df.drop(['IQ','JQ', 'J_yy', 'J_xy', 'J_yx', 'DRX', 'DRY', 'DRZ', 'DR'], axis=1, inplace=True)

        # indexing from 0 instead of 1 for VAMPIRE
        df['IT'] = df['IT']-1
        df['JT'] = df['JT']-1

        # make a tensor with zeros
        df[['J_12', 'J_13', 'J_21', 'J_22', 'J_23', 'J_31', 'J_32', 'J_33']] = 0

        # ----- SET THE HEISENBERG MATRIX VALUES -----

        df['J_xx'] = df['J_xx'] * eV
        df['J_22'] = df['J_xx']
        df['J_33'] = df['J_xx']

        # --------------------------------------------

        df.sort_values(by=['IT', 'JT', 'N1', 'N2', 'N3'], inplace=True)
        df.reset_index(inplace=True)
        df.drop(labels='index', axis=1, inplace=True)

        # save to UCF file
        with open(fout, 'w') as fwrite:
            df.to_csv(fwrite, mode='w', header=False, sep=" ", line_terminator='\n')


def main():
    sprkkr_to_vampire_ucf(path, system_name, n_atoms)

if __name__ == "__main__":
    main()

