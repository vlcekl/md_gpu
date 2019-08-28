# MD simulation code for a single GPU coded in CUDA
  
Allows simulations in NVE and NVT ensembles

Interaction potential terms
* Lennard-Jones non-bonded interactions
* Point charges
* Harmonic bonds and angles

Written in collaboration between Aaron Thompson (summer student at Vanderbilt University) and Lukas Vlcek

______

To do/done

* Neighbor list (Aaron)
* Add particle types to the neighbor list (Lukas)
* Bonded interactions input/output data organization (Lukas)
* Bonded particle list based on a simplified neighbor list (Aaron)
* Measurement - energy, pair distribution functions, dynamic properties (Lukas)
* Sorting algorithm SFCPACK based on: Anderson et al. J. Comp. Phys. 227 (2008) 5342-5359. (Aaron or who has time)
