# Single-nucleosome imaging unveils condensin and nucleosome-nucleosome interactions constrain chromatin to organize mitotic chromosomes.

Kayo Hibino, Sachiko Tamura, Yuji Sakai, Masatoshi Takagi, Masa A. Shimazoe, Toyoaki Natsume, Masato Kanemaki, Naoko Imamoto, 
Kazuhiro Maeshima


## Code to calculate molecular dynamics of coase-grained chromatins composing mitotic chromosomes

This Code is used in the molecular dynamics (MD) simulation of coase-grained chromatins composing mitotic chromosomes. \
Figures 5 and Movie S3-S5 in the paper were calculated using this code.

Extensible Simulation Package for Research on Soft Matter (ESPResSo) [1] is an MD package, which features a broad range of interaction potentials. 
ESPResSo is used as an MD simulator in this study as in our previous works [2,3].

MD time evolution programs of ESPResSo are written in C. The scripting language, Tcl, provides the interface between the user and the simulation engine. Therefore, the user may interact with the parallelized package core, as well as modify simulation parameters during runtime via Tcl commands. 

The tcl file "chromatin_motion.tcl" is the executable for the MD simulation.
The file "polymer.init" is the initial coordination file for chromatins and condensins.

One calculation usually takes several hours.\
An example output is shown in the output.dat file, \
where the first, second, and third columns are s, x-coordinates, and y-coordinates, respectively.

Since the code is written in Fortran, it should be executed on a machine with a Fortran compiler installed.


# Summary

Autophagy is an intracellular degradation process mediated by autophagosomes. The formation of autophagosomes involves dynamic morphological changes of the precursor phagophore, in which a disk-shaped membrane cisterna grows, bends into a cup-shaped intermediate, and finally becomes a spherical autophagosome. However, the physical mechanism behind these morphological changes remains elusive. This is in part because the shapes of phagophores (i.e., forming autophagosomes) are highly variable, and their standard shapes have not been characterized precisely. Here, we first determined the average shapes of phagophores by statistically investigating three-dimensional electron micrographs of more than 100 phagophores. The results showed that the cup-shaped structures adopted a characteristic morphology; they were longitudinally elongated, and the rim was catenoidal with an outwardly recurved shape. From the characteristic shape of the rim, we estimated the Gaussian modulus that determines the elastic properties of the phagophore membrane. Using the Gaussian modulus, we established a theoretical model of the shape of entire phagophores. The model quantitatively reproduced the average morphology observed by electron microscopy and revealed that the characteristic shape of phagophores (i.e., an elongated shape with a catenoidal rim) was primarily determined by the relative size of the open rim to the total surface area. These results suggest that autophagosomal membranes are highly flexible and that the morphological changes during autophagosome formation follow a stable path determined by elastic bending energy minimization.
