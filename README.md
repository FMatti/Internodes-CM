![](https://img.shields.io/badge/status-finished-green?style=flat-square)
![](https://img.shields.io/badge/licence-MIT-green?style=flat-square)
![](https://img.shields.io/badge/language-Python-blue?style=flat-square)
![](https://img.shields.io/badge/requirement-akantu-blue?style=flat-square)

# Internodes-CM

Implementation of the internodes method for two- and three-dimensional problems in contact mechanics.

## Introduction

This code served as a prototype for the inclusion of the method in [`akantu`](https://gitlab.com/akantu/akantu/-/tree/features/40-contact-using-the-internodes-method) in collaboration with [@Technici4n](https://github.com/Technici4n).

## Implementations
The INTERNODES method from [1] is implemented for static problems in contact mechanics [2]. The package [`akantu`](https://gitlab.com/akantu/akantu/) is used as a finite element engine to prime the algorithm, wherease the solution is obtained independently with `numpy` and `scipy`.

## Quick start

Clone the repository
```[bash]
git clone https://github.com/FMatti/Internodes-CM.git
cd Internodes-CM
```

Optionally create and activate a virtual environment
```[bash]
python -m venv .venv
source .venv/bin/activate (on Linux, macOS)
.venv\Scripts\activate.bat (on Windows)
```

Install the required packages
```[bash]
pip install --upgrade pip
pip install -r requirements.txt
```

Generate the meshes (replace the placeholders with appropriate values, i.e. for [DIMENSION] you should use 2 or 3 and for the filenames a valid string)
```[bash]
GMSH -[DIMENSION] -o [MESHNAME].msh [GEONAME].geo
```

## File structure
The main implementation is found in the file `contact_mechanics_internodes.py`. Some examples can be found in the jupyter notebooks.

```
Internodes-CM
│   README.md
|   contact_mechanics_internodes.py         (source file of all implementations)
|   example2d.ipynb                         (two-dimensional example)
|   example3d.ipynb                         (three-dimensional example)
|   helper.py                               (helper file containing utility functions) 
|   plots.ipynb                             (notebook used to generate plots)
|   requirements.txt                        (list of required packages)
|   test_contact_mechanics_internodes.py    (tests for source file)
|
└───material                                (akantu material file)
|   |   material.dat                       
|
└───mesh                                    (geometry files to generate meshes)
|   |   contact2d_circle.geo                 
|   |   contact2d_circle_circle.geo
|   |   contact2d_plane_circle.geo
|   |   contact3d_plane_sphere.geo
|   |   contact3d_sphere.geo
|   |   contact3d_sphere_sphere.geo 
|
└───presentation
|   |   ...
|
└───report
|   |   ...
```

## References

[1] Simone Deparis, Davide Forti, and Alfio Quarteroni. A rescaled localized radial basis function interpolation on non-cartesian and nonconforming grids. SIAM Journal on Scientific Computing, 36(6):A2745–A2762, 2014. doi: 10.1137/130947179

[2] Yannis Voet, Guillaume Anciaux, Simone Deparis, and Paola Gervasio. The internodes method for applications in contact mechanics and dedicated preconditioning techniques. Computers & Mathematics with Applications, 127:48–64, 2022. doi: 10.1016/j.camwa.2022.09.019
