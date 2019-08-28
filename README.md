# MD simulation code for a single GPU coded in CUDA
  
Allows simulations in NVE and NVT ensembles

Interaction potential terms
* Lennard-Jones non-bonded interactions
* Point charges
* Harmonic bonds and angles

Written in collaboration by Aaron Thompson (summer student at Vanderbilt University) and Lukas Vlcek

______
To do/done

1a) Neighbor list (Aaron)
1b) Add particle types to the neighbor list (Lukas)
2a) Bonded interactions input/output data organization (Lukas)
2b) Bonded particle list based on a simplified neighbor list (Aaron)
3) Measurement - energy, pair distribution functions, dynamic properties (Lukas)
3) Sorting algorithm SFCPACK based on: Anderson et al. J. Comp. Phys. 227 (2008) 5342-5359. (Aaron or who has time)
