"""Convert the SPR-KKR exchange interactions to vampire."""

import pandas as pd
import numpy as np



# ========================== USER INPUT =============================

# SPR-KKR calculation directory
path = "H:/2D/Cr2Te3/bulk/sprkkr_imitating_experimental/Bi2Te3" #"H:/2D/Cr2Te3/bulk/transport/PBE/NKTAB_1000"

# the seed name in SPR-KKR (the name of the system)
system_name = "POSCAR"
n_atoms     = 20

include_dmi = True

# ===================================================================




def sprkkr_to_vampire_ucf(path, system_name, n_atoms, include_dmi=True):
    """Convert exchange and DMI from SPR-KKR into VAMPIRE UCF file."""

    eV = 1.602e-19  #J

    # SPR-KKR J and DMI .dat files in
    f_exchange_in = f"{path}/{system_name}_JXC_XCPLTEN_Jij.dat"
    f_dmi_in      = f"{path}/{system_name}_JXC_XCPLTEN_Dij.dat"

    # UCF file out
    fout = f"{path}/vampire.UCF_from_{system_name}"

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

    df['J_11'] = df['J_11'] * eV
    df['J_22'] = df['J_11']
    df['J_33'] = df['J_11']

    if include_dmi == True:
            df_dmi = pd.read_csv(f_dmi_in, skiprows=12+n_atoms, delim_whitespace=True, names= ['IT', 'IQ', 'JT', 'JQ', 'N1', 'N2', 'N3', 'DRX', 'DRY', 'DRZ', 'DR', 'DX', 'DY', 'DZ'] )

            # see e.g. Coey for definition of DMI -> DMI tensor components ( (0, -Dz, Dy), (Dz, 0, -Dx), (-Dy, Dx, 0) )
            df['J_12'] += -df_dmi['DZ'] * eV
            df['J_13'] += df_dmi['DY'] * eV
            df['J_21'] += df_dmi['DZ'] * eV
            df['J_23'] += -df_dmi['DX'] * eV
            df['J_31'] += -df_dmi['DY'] * eV
            df['J_32'] += df_dmi['DX'] * eV

    # --------------------------------------------

    df.sort_values(by=['IT', 'JT', 'N1', 'N2', 'N3'], inplace=True)
    df.reset_index(inplace=True)
    df.drop(labels='index', axis=1, inplace=True)

    # save to UCF file
    with open(fout, 'w') as fwrite:
        fwrite.write(f"# Interactions\n{df.shape[0]} tensorial\n")
        df.to_csv(fwrite, mode='w', header=False, sep=" ", line_terminator='\n')



def main():
    sprkkr_to_vampire_ucf(path, system_name, n_atoms, include_dmi)


if __name__ == "__main__":
    main()
