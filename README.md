# Introduction to opacfractal

This `opacfractal` package allows to compute optical properties of fractal dust aggregates 
by means of a statistical distribution model of monomers.
This code has been tested by comparing the results obtained by the T-Matrix Method in [Tazaki et al. 2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...823...70T) for scattering matrix elements and [Tazaki and Tanala 2018](https://ui.adsabs.harvard.edu/abs/2016ApJ...823...70T) for the opacity and asymmetry parameter.

`opacfractal` has been also implemented in a commad-line opacity tool [optool](https://github.com/cdominik/optool.git).

# Terms of use

`opacfractal` is distributed under the [MITlicense](https://opensource.org/licenses/MIT) and can be used, changed
and redistributed freely. If you use this package to publish papers, please cite the relevant papers.
Because this code includes some formulations by previous authors and can be called by the code options (`iqsca`, `iqcor`, and `iqgeo`). Depending on each code option, the relevant paper is as follows.
 - `iqsca=1` : [Tazaki et al. 2016, ApJ, 823, 70](https://ui.adsabs.harvard.edu/abs/2016ApJ...823...70T)
 - `iqsca=2` : [Botet et al. 1997, ApOpt, 36, 8791](https://ui.adsabs.harvard.edu/abs/1997ApOpt..36.8791B)
 - `iqsca=3` : [Tazaki & Tanaka 2018, ApJ, 860, 79](https://ui.adsabs.harvard.edu/abs/2018ApJ...860...79T)
 - `iqcor=1` : [Tazaki et al. 2016, ApJ, 823, 70](https://ui.adsabs.harvard.edu/abs/2016ApJ...823...70T)
 - `iqcor=2` : [Berry & Percival 1986, AcOpt, 33, 577](https://ui.adsabs.harvard.edu/abs/1986AcOpt..33..577B)
 - `iqcor=3` : [Botet et al. 1995, JPhA, 28, 297](https://ui.adsabs.harvard.edu/abs/1995JPhA...28..297B)
 - `iqgeo=2` : [Okuzumi et al. 2009, ApJ, 707, 1247](https://ui.adsabs.harvard.edu/abs/2009ApJ...707.1247O) 
 - `iqgeo=3` : [Tazaki (submitted to MNRAS)](https://github.com/rtazaki1205/geofractal)

# How to use it? 

There are three examples of calling routines: `call.f90`, `call2.f90`, and `call3.f90`. The user needs to specify one of them in `Makefile`.
`call.f90` produces an output file `out_smat.dat`, which contains scattering elements at a single wavelength. `call2.f90` produces `out_opc.dat`, which contains wavelength dependent opacities. `call3.f90` produces `dustkapscatmat.inp`, which contains produces scattering matrix elements and opacities averaged over size distribution, and this file can be used as a input file of [RADMC-3D](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/). In `call2.f90` and `call3.f90`, refractive index of astronomical silicate is used (Draine 2003).

In the calling routine, the user needs specify following input parameters:

- `df` : Fractal dimension (1 ≦ df ≦ 3)
- `k0` : Fractal prefactor
- `PN` : Number of monomers (1 ≦ PN)
- `R0` : Monomer radius (micron)
- `lmd` : Wavelength (micron)
- `ref` : Complex refractive index
- `nang` : Number of angle mesh (=91) 

In addition, the user also needs to specify following four options:

- Light scattering solver   
  `iqsca=1` : Rayleigh-Gans-Debye theory  
  `iqsca=2` : Mean field theory  
  `iqsca=3` : Modified mean field theory  
- The two-point correlation function of fractal aggregates  
  `iqcor=1` : The Gaussian cut-off model  
  `iqcor=2` : The exponential cut-off model  
  `iqcor=3` : The fractal dimension cut-off model  
- Geometric cross section of fractal aggregates (needed only when `iqsca=3`)  
  `iqgeo=1` : The characteristic cross sections  
  `iqgeo=2` : Empirical formula by [Okuzumi et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...707.1247O)  
  `iqgeo=3` : Analytical formula by Tazaki (in prep.) from [geofractal code](https://github.com/rtazaki1205/geofractal)
- standard output  
  `iquiet=0` : show standard output  
  `iquiet=1` : suppress standard output (including warnings)  
	
I recommend following set of options: `iqsca=3`,`iqcor=1`,`iqgeo=3` (default).  
For the two-point correlation function, now I suggest to use `iqcor=1` (Gaussian type) because it is numerically stable and is also possible to reproduce optical properties of fractal aggregates. 

To run the code, first of all, compile the codes by
```
make
```
This will create an executable file `fracsca.x`. Then, perform
```
./fracsca.x
```
As a result, the output file `out_smat.dat` or `out_opc.dat` is created. 

# Limitation of the code 

`opacfractal` is a code based on an approximate technique to solve light scattering, and therefore, it is not applicable for all possible parameters. Beyond the limitation of the code, I suggest the users to use more rigorous numerical approach, such as [MSTM](https://www.eng.auburn.edu/~dmckwski/scatcodes/) and [DDSCAT](http://ddscat.wikidot.com/).

I summarize some limitations of `opacfractal` for each light scattering solver (`iqsca`). The limitation can be simply judged by the phase shift induced by an aggregate (Equation 9 in Tazaki & Tanaka 2018; see also Section 3.2 in this paper).

- `iqsca=1`   
  All outputs would be physically reasonable for the phase shift <~ 1.   
 
- `iqsca=2`  
 The extinction cross section would be calculated without limitation.  
 However, the other outputs would be reliable for the phase shift <~ 1.  

- `iqsca=3`  
  The extinction cross section would be calculated without limitation. Scattering and absorption cross sections could be calculated for the phase shift > 1, however too large phase shift may cause some problem. The asymmetry parameter and the sattering matrix elements would be reliable for the phase shift <~ 1.  

The code also neglects the enhancement of absorption opacity due to the presence of adjacent monomers in an aggregate (see e.g., Section 5 in [Tazaki & Tanaka 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...860...79T)). Therefore, the code underestimates the absorptin opacity in the Rayleigh domain particularly when monomers have high refractive index, such as amorphous carbon.

For safety, the code always returns a value of the phase shift, and therefore, the user can check whether a set of inputted parameters is within the limitation or not. The code also produces a warning when the phase shift is above unity, although this warning is suppressed when `iquiet = 1`.

# Acknowledgement 

`opacfractal` contains following subroutines provided by other authors.
I thank the authors for the availability of the subroutines.  

A list of subroutines by other authors:  

`lorenz_mie` : Modified version of BHMIE code written by Bruce T. Draine  
`renorm_mie` : Modified version of BHMIE code written by Bruce T. Draine  
`gamma`      : Zhang, S. and Jin, J. (1996) "Computation of Special Functions.  
`chgm`       : Zhang, S. and Jin, J. (1996) "Computation of Special Functions.  
`lpmns`      : Zhang, S. and Jin, J. (1996) "Computation of Special Functions.  
`lpn`        : Zhang, S. and Jin, J. (1996) "Computation of Special Functions.  
`ludcmp`     : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"    
`lubksb`     : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"  
`gauleg`     : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"  


# History

*Version 3.0*
- Public release 
