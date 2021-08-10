
# List of Bravais Lattices
1. Simple cubic: 
    * alpha-N2, one more 
2. BCC:
    * Many transition metals, i.e. K, however these all correspond to 
    elemental crystals 
    * SiO2 zeolite structure 
3. FCC:
    * silicon, germanium, diamond
    * cubic Boron nitride, copper, MgO, NaCl, PbS, LiH  
4. Simple tetragonal:
    * TiO2 rutile, add another
5. BC tetragonal:
    * TiO2 anatase
6. Hexagonal:
    * CdSe, unbuckled graphene
7. Rhombohedral
    * MoS2, other
8. Simple orthorhombic
    * Pick 2
9. Base-centred orthorhombic
    * NEED
10. Body-centred orthorhombic
    * NEED
11. Face-centred orthorhombic
    * NEED
12. Simple monoclinic
    * Try LiSn and PdCl2
13. Base-centred monoclinic
    * Ga2O2, other?
14. Triclinic 
    * NEED

###### Sources
* Materials Project: https://materialsproject.org
* AFLOW: http://aflow.org/search/advanced.php  **Find the link from Rui**

#Tests
* Tests should cover each Bravais lattice, listed above. 
* Also consider splitting in other types; covalent, ionic, sheets, metals
* Different temperatures: At 0k, 300k, 1000k
* SCF and SCC 

## Test Settings
* Use high precision
* Bohr is better than angstrom (internal unit conversion factor will change in 2020)
* Better to use fractions rather than decimals i.e. 1/3 instead of 0.333


## Trinclinc
Should put one crystal system in as a trinclinic lattice and confirm it gives the
correct results 


### Translational Invariance 
* Translate all atoms in the cell by some arbitrary amount 
* Expect the energy to be invariant w.r.t. translation. 
* TODO Finish adding orthorhombic, monoclinic and triclinic systems 


### Elemental Crystals
* Benchmark band structure and DOS for every elemental crystal that's parameterised
    * See README in the corresponding subfolder 


### Converged Total Energy of Primitive Cell 
* Using CP2K's method of cutting off the real-space potential.

###### Things to Check
* Does the energy converge?
* See how closely our energies agree with those of CP2K?
* How well do the band gaps a) agree with CP2K and b) experiment? 
* Do the band structures a) look sensible and b) in broad agreement with CP2K? 

######Issues
* What value should one cut off at? CP2K uses rcuta + rcutb, which comes from the basis functions, 
but I can't see exactly where the values are set, therefore I don't know what they are.
    * The way to address this is to print out for every pair and tabulate   
    => Compile CP2K from source and run. 
    * Not so viable if it's shell-resolved 
* Because 1/r**3 is is slowly decaying, one might expect truncation of the real-space potential
at any point is going to introduce a non-negligible error. 



### Energy vs volume 
* For crystals with a single lattice parameter, compute energy vs volume curves.
* This means any cubic (and potentially rhombohedral) structure. 
* Material candidates: 
    * Zincblende: Silicon, α-AgI, β-BN, diamond, CuBr, LiH(?)
    * Rock-salt: NaCl, LiF, PbS, MgO

######Issues
* Cubic (?) boron nitride showed good energy vs volume curve.
* See if these issues have been alleviated:
    * MgO energy vs volume curve just increased in energy w.r.t. lattice constant
    * Rutile energy vs volume curve also increased in energy w.r.t. lattice constant BUT to a ridiculous 
    energy. 



### Energy per atom
* Energy per atom is conserved as one moves from the primitive unit cell to conventional unit cell, and
subsequently super cells of increasing size.

* Can take any system and increase the super cell integers. k-sampling must be converged in each case 
although one can reduce the sampling density as the cell gets larger. 

* Compare the energy per atom of a primitive cell with converged k-sampling, with a very large cell (100-200 atoms) 
calculated at the Gamma point. Should also be equal to the molecular energy. 



### Molecules in a Box
* The energy of a molecule in a sufficiently large periodic box computed at Gamma should agree with 
the energy of a molecule in the molecular implementation. 

* Should already be a test for Co2 
    * I wrote some python to generate that box, but might be easier to use pymatgen
* Should be invariant w.r.t. position in the box.



 