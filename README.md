# Single-nucleosome imaging unveils condensin and nucleosome-nucleosome interactions constrain chromatin to organize mitotic chromosomes.

Kayo Hibino, Sachiko Tamura, Yuji Sakai, Masatoshi Takagi, Masa A. Shimazoe, Toyoaki Natsume, Masato Kanemaki, Naoko Imamoto, 
Kazuhiro Maeshima


## Code to calculate molecular dynamics of coase-grained chromatins composing mitotic chromosomes

This Code is used in the molecular dynamics (MD) simulation of coase-grained chromatins composing mitotic chromosomes. \
Figures 5 and Movie S3-S5 in the paper were calculated using this code.

Extensible Simulation Package for Research on Soft Matter (ESPResSo) [1] is an MD package, which features a broad range of interaction potentials. 
ESPResSo is used as an MD simulator in this study as in our previous works [2-4].

MD time evolution programs of ESPResSo are written in C. The scripting language, Tcl, provides the interface between the user and the simulation engine. Therefore, the user may interact with the parallelized package core, as well as modify simulation parameters during runtime via Tcl commands. 
Detailed guide to ESPResSo is found in the papers [1,4].


The tcl file "chromatin_motion.tcl" is the executable for the MD simulation. \
The file "polymer.init" is the initial coordination file for chromatins and condensins.

Execute at the command prompt with the following command;\
./Espresso chromatin_motion.tcl


One calculation usually takes several hours.\
Four output files are generated when the calculation is run.\
One is the "polymer.vtf" file, which is a time series of Cartesian coordinates for each of the 5000 chromatins.\
The others are the "MSD0", "MSD1", and "MSD2" file, which is the mean squared displacement of nucleosomes in the overall, axial, and peripheral regions, respectively.


An example output is shown in the "polymer.vtf" file, which is the Cartesian coordinates of each of the 5000 chromatins.


1. ESPResSo on Github; https://github.com/espressomd/espresso.
2. Y. Sakai, M. Tachikawa, A. Mochizuki, A simple model for eukaryotic chromosome segregation, Phys. Rev. E 2016, 94(4-1):042403.
3. Y. Sakai, A. Mochizuki, K. Kinoshita, T. Hirano, M. Tachikawa, Modeling the functions of condensin in chromosome shaping and　segregation, PLoS. Comp. Biol. 2018, 18;14(6):e1006152.
4. Y. Sakai, Molecular Dynamics Simulations of Condensin-Mediated Mitotic Chromosome Assembly, Meth. Mol Biol. 2019, 2004:319-334.
