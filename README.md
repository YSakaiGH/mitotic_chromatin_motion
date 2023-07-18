Single-nucleosome imaging unveils condensin and nucleosome-nucleosome interactions constrain chromatin to organize mitotic chromosomes.

Yuji Sakai, Satoru Takahashi, Ikuko Koyama-Honda, Chieko Saito, Noboru Mizushima

## Code to calculate stable membrane shapes

This Code is used in the model calculations to determine the membrane shapes in the paper,\
https://www.biorxiv.org/content/10.1101/2022.07.20.500884v1. \
Figures 6, 7 and S3 in the paper were calculated using this code.

The membrane shape is determined from the elastic bending energy using the Euler-Lagrange equation (S3).\
Eq. (S3) is solved under the boundary conditions (S4)-(S5) to obtain the membrane coordinates (X, Z) at each s.\
A calculation with one fixed boundary conditions takes roughly a few minutes.\
An example output is shown in the output.dat file, \
where the first, second, and third columns are s, x-coordinates, and y-coordinates, respectively.

Since the code is written in Fortran, it should be executed on a machine with a Fortran compiler installed.


# Summary

Autophagy is an intracellular degradation process mediated by autophagosomes. The formation of autophagosomes involves dynamic morphological changes of the precursor phagophore, in which a disk-shaped membrane cisterna grows, bends into a cup-shaped intermediate, and finally becomes a spherical autophagosome. However, the physical mechanism behind these morphological changes remains elusive. This is in part because the shapes of phagophores (i.e., forming autophagosomes) are highly variable, and their standard shapes have not been characterized precisely. Here, we first determined the average shapes of phagophores by statistically investigating three-dimensional electron micrographs of more than 100 phagophores. The results showed that the cup-shaped structures adopted a characteristic morphology; they were longitudinally elongated, and the rim was catenoidal with an outwardly recurved shape. From the characteristic shape of the rim, we estimated the Gaussian modulus that determines the elastic properties of the phagophore membrane. Using the Gaussian modulus, we established a theoretical model of the shape of entire phagophores. The model quantitatively reproduced the average morphology observed by electron microscopy and revealed that the characteristic shape of phagophores (i.e., an elongated shape with a catenoidal rim) was primarily determined by the relative size of the open rim to the total surface area. These results suggest that autophagosomal membranes are highly flexible and that the morphological changes during autophagosome formation follow a stable path determined by elastic bending energy minimization.
