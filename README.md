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


The tcl file "chromatin_motion.tcl" is the executable for the MD simulation.
The file "polymer.init" is the initial coordination file for chromatins and condensins.

One calculation usually takes several hours.\
An example output is shown in the "test.vtf" file, where the first, second, and third columns are x, y, and z-coordinates of each chromatin, respectively.

Since the code is written in Fortran, it should be executed on a machine with a Fortran compiler installed.



[1] ESPResSo on Github; https://github.com/espressomd/espresso. \
[2] Y. Sakai, M. Tachikawa, A. Mochizuki, A simple model for eukaryotic chromosome segregation, Phys. Rev. E 94, 042403 (2018).\
[3] Y. Sakai, A. Mochizuki, K. konoshita, T. hirano, M. Tachikawa, Modeling the functions of condensin in chromosome shaping and \
　　　　segregation, PLoS. Comp. Biol. \
[4] Y. Sakai, 
