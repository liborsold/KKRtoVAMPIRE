# KKRtoVAMPIRE
Convert the Heisenberg Hamiltonian calculated with SPR-KKR to the VAMPIRE format.

- see `./src/KKRtoVAMPIRE/sprkkr_to_vampire.py`

Needed are the files: seedname.pot_new, seedname_SCF.out, seedname_JXC_XCPLTEN_Jij.dat and seedname_JXC_XCPLTEN_Dij.dat.
Specify if DMI and anisotropy should be included and if the interactions should be croppped.
Tested on .pot file with  >> FORMAT     9 (18.01.2019) <<
