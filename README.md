# KKRtoVAMPIRE

**<i>Convert the Heisenberg Hamiltonian calculated with SPR-KKR to the VAMPIRE format.</i>**

_If you find this package useful, please cite_ [Q. Guillet*, L. Vojáček* _et al._, Phys. Rev. Materials **7**, 054005 (2023)](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.7.054005).

## Usage

See `./examples/KKR_to_VAMPIRE_example.ipynb` for the example of use.

Needed files are: `seedname.pot_new`, `seedname_SCF.out`, `seedname_JXC_XCPLTEN_Jij.dat` and `seedname_JXC_XCPLTEN_Dij.dat`.

Specify if DMI and anisotropy should be included and if the interactions should be croppped. Tested on .pot file with  _FORMAT     9 (18.01.2019)_.
