# non collinear spin structure

The purpose of this work is to easily manipulate crystal structure with non collinear spin and to compute VASP workflows. It is done within the pymatgen/atomate framework. Non collinear spin are made by fusioning a Structure object with Magmom object for each site/species. In this work, the site localisation of spin is prefered because species localisation, with enable disordered strucutre, is harder not required in my analysis. 
Magmom takes 2 arguments : intensity and orientation. these 2 follow the SAXIS standard used by vasp. It follows the crystal vector axis.

# NON COLLINEAR STRUCTURE OBJECT
There are some function to easily make non collinear spin (from user array or predefined spin ice case). In the future, I can consider to implement other mean to establish non collinear situations.

# PLOT
It is possible to plot non collinear spins with vtk. It is particularly usefull for non collinear cases because vizualisation is tricky, especially in the case of spin ice. In the future, it would be possible to show angle from a certain plane, ...

# IMAGING BETWEEN 2 SPIN STATES
There is the possibility to establish images between 2 spin states, in a similar fashion than NEB. The energy between these 2 spin states can be computed with 2 VASP workflows (1) to compute magnetic anisotropy and (2) to compute energy profile and saddle point.

(-> COMPUTE U WITH COCOCCIONI. A workflow that enables to compute U with VASP via the Cococcioni method. Not for non collinear spin)
