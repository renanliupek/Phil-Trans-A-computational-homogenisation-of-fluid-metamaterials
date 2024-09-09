# cohoflu

GitLab repository: https://gitlab.tue.nl/20204732/philosophical-transactions-a-computational-homogenisation-of-metamaterials

dataset DOI: https://data.4tu.nl/datasets/0c31cd57-7ea1-4587-84e7-b9b75ff5fa2e

### Renan Liupekevicius Carnielli [r.liupekevicius.carnielli@tue.nl](mailto::r.liupekevicius.carnielli@tue.nl)

Data Availability of the article "An efficient multiscale method for subwavelength transient analysis of acoustic metamaterials"

ASSOCIATED PEER-REVIEWED PUBLICATION

[https://doi.org/https://royalsocietypublishing.org/doi/10.1098/rsta.2023.0368
](https://royalsocietypublishing.org/doi/10.1098/rsta.2023.0368)
Author: Renan Liupekevicius Carnielli, TU Eindhoven, April 2024

## comsol.zip folder description (DOWNLOAD FROM https://data.4tu.nl)

This zip folder contains the comsol models for the direct numerical simulations (DNS) and the homogenised macro-scale model (HOM) shown in the article. 

-*dispersion DNS* is the comsol model that computes the Bloach Analysis (figure 5(a)) and the local resonance modes (figure 4);

-*transmission DNS* is the comsol model that computes the DNS of the transmission problem, figure 6(c) dashed curve.

-*transmission homogenized* is the comsol model that computes the homogenised transmission problem, figure 6(c) solid curve.

## Note

1- run `start.m` to include the the path to the folders *fun*, *tensorlab* and *meshes*.

2- export comsol 5.4 mesh file `.mphtxt` of your favorite geometry to the folder *meshes*, or use the mesh file provided.

3- FATAL: remember to rename the exported comsol mesh `.mphtxt` to `.txt`.

4- the files `computedworkspace_Liang.mat` and `computedworkspace_Liang_willis.mat` are the output workspaces of `main.m` for: the symmetric unit cell, figure 3(a); and the asymmetric unit cell, figure 7, respectively.

## `main.m` file description

1- `main.m` computes the homogenised material coefficients of equations (3.30) and equation (3.31).

2-  The input for `main.m` are: the mesh file in the folder *meshes* and the material properties of the saturating fluid which can be defined within the script `main.m`.

3-  The appropriate mesh is created by exporting a text file `.mphtxt` from comsol 5.4, see 'Note' 2nd and 3rd item.

## `effective_material_properties.m` file description
This script computes the frequency-domain material properties of equation (4.2) assuming you ran `main.m`, or loaded the provided workspaces.

## *fun* folder
Developed functions to assist the computations in `main.m`.

## *tensorlab* folder
Collection of functions and classes that define a tensor object.

## *meshes* folder
It contains the mesh files of the two unit cells shown in the article or any other unit cell you may want to homogenise. The symmetric unit cell is taken from Liang, Z., & Li, J. (2012). Extreme acoustic metamaterial by coiling up space. Physical review letters, 108(11), 114301.
