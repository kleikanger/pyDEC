pyDEC
=====

The aim of this code is to perform DEC calculations on systems with periodic boundary conditions.
We plan to interface CP2K to generate MO's and LSDALTON to calculate localized Wannier orbitals and to run DEC calculations.

Python Library.
===============

Brief description
-----------------

Works as a link between CP2K and LSDALTON, to generate MO's + Fock
matrices with CP2K that can be ported to a DALTON simulation.

Functionality
-------------

-  Runs CP2K SCF routines to generate MO's.
-  Converts MO's to DALTON format.
-  Generates Basis input + system input (.mol files) to be sure that
   LSDALTON and CP2K runs with the same MO's.
-  Basis input formats

   -  DALTON
   -  CP2K

-  System Input Format

   -  DALTON
   -  CP2K (not implemented)

```python
>> cd ~/Dropbox/phd/pyCode/pyDEC/
```
```
/Users/leik/Dropbox/phd/pyCode/pyDEC
```

Porting basis files:
--------------------

Example: Basis STO-6G from DALTON to CP2K format.

```python
>> from basis import Basis
>> basis = Basis()
>> basis.read_basis_set('STO-3G', 'STO-3G', 'DALTON')  # inp: inputfile, basisname, inputformat
>> basis.read_basis_set('STO-6G', 'STO-6G', 'DALTON')  # inp: inputfile, basisname, inputformat
```
```
Reading basis from: STO-3G
Reading basis from: STO-6G
```

The basis sets 'STO-3G' and 'STO-6G' is now stored in the Basis object
in a list of BasisSet objects. In the next line we write to file.

```python
>> basis.print_basis_sets('STO_3G+STO_6G', 'CP2K')  # inp: output filename, format
```
```
Writing basis sets to file <STO_3G+STO_6G> in <CP2K> format.
Basis set <STO-3G> written to file.
Basis set <STO-6G> written to file.
```

Results:

'''python
>> !head -n26 STO-3G
'''
```
$Basis = STO-3G
$
$ Elements supported
$ H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu 
$ Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Rh Rh Pd Ag Cd 
$
$  REFERENCE
$ Elements      Contraction                       References                     
$ H - He: (3s)         -> [1s]       W.J. Hehre, R.F. Stewart and J.A. Pople,    
$ Li - Ne: (6s,3p)     -> [2s,1p]    J. Chem. Phys. 2657 (1969).                 
$ Na - Ar: (9s,6p)     -> [3s,2p]    W.J. Hehre, R. Ditchfield, R.F. Stewart,    
$                                    J.A. Pople, J. Chem. Phys.  2769 (1970).    
$ K,Ca - : (12s,9p)    -> [4s,3p]    W.J. Pietro, B.A. Levy, W.J. Hehre and R.F. 
$ Ga - Kr                            Stewart, Inorg. Chem. 19, 2225 (1980). 
$ Sc - Zn: (12s,9p,3d) -> [4s,3p,1d] W.J. Pietro and W.J. Hehre, J. Comp. Chem.  
$ Y - Cd: (15s,12p,6d) -> [5s,4p,2d] 4, 241 (1983).                              
$***********************************************************************
A 1
$ HYDROGEN     (3S) -> [1S]                                      
$ S-TYPE FUNCTIONS
     3    1    0
    3.4252509  0.15432897
    0.6239137  0.53532814
    0.1688554  0.44463454
A 2
```

```python
>> !head -n17 STO_3G+STO_6G
```
```
#
 # CP2K input file.
# This file is automatically generated by PROGNAME.
#
# ----------------- 
#
H STO-3G 
1
1 0 0 3 1 
3.4252509      0.15432897     
0.6239137      0.53532814     
0.1688554      0.44463454     
#
# ----------------- 
#
He STO-3G 
1
```

Load system information
-----------------------

Loading system from LSDALTON .mol file.

```python
>> !cat H2.mol
```
```
BASIS
STO-3G Aux=Ahlrichs-Coulomb-Fit
LDA molecular hessian without symmetry
Ethane LDA molecular hessian without symmetry
Atomtypes=1 Charge=0 angstrom 
charge=1.    atoms=2
H     0.50000000000      0.0000000000     0.0000000000
H     -0.50000000000      0.0000000000     0.0000000000
a1 = 2.0 0.0 0.0 active
a2 = 0.0 2.0 0.0 inactive
a3 = 0.0 0.0 2.0 inactive
'''

'''python
from System import System
s = System()
s.read_system_info('H2.mol', 'LSDALTON')
```

```
Reading DALTON .mol file: H2.mol
```

Result:

```python
>> for x in s.atoms:
     print('basis: ' + x.basisname, end=', ')
     print('a: %s' % x.a, end=' ')
     print('pos:', end=' '); print(x.position, end=', ')
     print('ghost: %s' % x.ghost, end=', ')
     print('charge: %s' % x.charge, end='\n') 
>> print('latvec:')
>> for x in s.latvec:
     print(x)
```

```
 basis: STO-3G, a: 1 pos: [0.5, 0.0, 0.0], ghost: False, charge: 1.0
 basis: STO-3G, a: 1 pos: [-0.5, 0.0, 0.0], ghost: False, charge: 1.0
 latvec:
 [2.0, 0.0, 0.0]
 [0.0, 2.0, 0.0]
 [0.0, 0.0, 2.0]
```

Manipulating the MO's
---------------------

class MOs: - Read MO's from CP2K file. - Permute to DALTON format. -
Print to DALTON input file.

```python
 from MOs import MOs
 # o = MOinfo('He_bulk2-RESTART.wfn')
 o = MOs('Ne_bulk2-RESTART.wfn', s, basis)  # inp: input file (CP2K restart file)
 o.transform_cp2k_to_dalton()
```

```
 Reading binary file Ne_bulk2-RESTART.wfn.
 [10, 10, 11, 20]
 Writing binary file orbitals_in.u.
```

The MO coefficients are read from file permuted to dalton order and
written to orbitals\_in.u. The permutations of the coefficients are
found by calling:

python
```
 perm = basis.get_mo_transformation('CP2K', 'DALTON', s.atoms)  # inp: format from, format to, __Atom objects
 print(perm)
```

```
    [0, 1]
```

##Example with several basis sets:

We use different basis sets for the atoms.

```python
 !cat N2.mol
```

```
 ATOMBASIS
 LDA molecular hessian without symmetry
 Ethane LDA molecular hessian without symmetry
 Atomtypes=1 Charge=0 angstrom 
 charge=10.    atoms=1    basis=aug-cc-pVTZ
 Ne     0.50000000000      0.0000000000     0.0000000000
 charge=10.    atoms=1    basis=STO-6G
 Ne     -0.50000000000      0.0000000000     0.0000000000
 a1 = 2.0 0.0 0.0 active
 a2 = 0.0 2.0 0.0 inactive
 a3 = 0.0 0.0 2.0 inactive
```


Load the .mol file and the basis sets.

```python
 systemN2 = System()
 systemN2.read_system_info('N2.mol', 'LSDALTON')
```
```
 Reading DALTON .mol file: N2.mol
```

```python
 basisN2 = Basis()
 basisN2.read_basis_set('aug-cc-pVTZ', 'aug-cc-pVTZ', 'DALTON')
 basisN2.read_basis_set('STO-6G', 'STO-6G', 'DALTON')
```
```
 Reading basis from: aug-cc-pVTZ
 Reading basis from: STO-6G
```

The permutations from CP2K to DALTON MOs:

```python
 print(basisN2.get_mo_transformation('CP2K', 'DALTON', systemN2.atoms))
```

```
 [44, 45, 35, 7, 3, 43, 41, 42, 27, 25, 26, 6, 4, 5, 2, 0, 1, 36, 37, 38, 39, 40, 20, 21, 22, 23, 24, 8, 9, 10, 11, 12, 28, 29, 30, 31, 32, 33, 34, 13, 14, 15, 16, 17, 18, 19, 50, 46, 49, 47, 48]
```

PS: Have printed the basis sets and manually checked that this is
correct.

Port MO's from CP2K to DALTON:
------------------------------

We have the files:

```
 CP2K_aug-cc-pvZT      MOLECULE.INP          t_c_g.dat
 He2.mol               aug-cc-pVTZ           transformMOS.py
 He_bulk2-RESTART.wfn  localized_orbitals.u
 LSDALTON.INP          orbitals_in.u
```

We look at the system (MOLECULE.INP):

```
 BASIS
 aug-cc-pVTZ
 Atomtypes=1
 Charge=2. Atoms=2
 He     0.50000000000      0.0000000000     0.0000000000
 He     -0.50000000000      0.0000000000     0.0000000000
```

-  We first run our code to generate the basis set input. NB: It is a
   good idea to generate both CP2K and DALTON files to make sure that
   the order of the orbitals is correct. The ao's does not need to be
   set up in any certain order, and this could result in wrong
   MO'permutation arrays because of the way the code is set up.
-  Molecular calculation with CP2K:
-  Note that we need to add the keyword ...:DFT:ADDED\_MOS 44 to the
   CP2K input. This keyword is needed to calculate virtual MO's (2 occ +
   44 virt) and write them to restart.
-  We run CP2K to get the restart file He\_bulk2-RESTART.wfn from the
   CP2K simulation.
-  To transform MO's to dalton format we run the transfomMOs.py script.

transformMOS.py:
```python
 import sys
 sys.path.append('/Users/leik/Dropbox/phd/pyCode/pyDEC/')
 from MOs import MOs
 from System import System
 from basis import Basis
 
 basis = Basis()
 basis.read_basis_set('aug-cc-pVTZ', 'aug-cc-pVTZ', 'DALTON')
 
 system = System()
 system.read_system_info('He2.mol', 'LSDALTON')
 
 o = MOs('He_bulk2-RESTART.wfn', system, basis)
 o.transform_cp2k_to_dalton()
```

We then run LSDALTON with the following input to calculate the molecular
orbitals (LSDALTON.INP)

```
 **WAVE FUNCTION
 .HF
 **LOCALIZE ORBITALS
 .PSM
 1 1
 .Only Loc
 *END OF INPUT
'''
