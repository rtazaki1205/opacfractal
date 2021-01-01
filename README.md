# Introduction to opacfractal

This package allows to compute optical properties of fractal dust aggregates 
by means of a statistical distribution model of monomers.

# Terms of use

`opacfractal` is distributed under the [MITlicense](https://opensource.org/licenses/MIT) and can be used, changed
and redistributed freely. If you use this package to publish papers, please cite the relevant papers depending on code options (`iqsca` and `iqcor`):

 - `iqsca = 1`: [Tazaki et al. 2016, ApJ, 823, 70](https://ui.adsabs.harvard.edu/abs/2016ApJ...823...70T)
 - `iqsca = 2`: [Botet et al. 1997, ApOpt, 36, 8791](https://ui.adsabs.harvard.edu/abs/1997ApOpt..36.8791B)
 - `iqsca = 3`: [Tazaki & Tanaka 2018, ApJ, 860, 79](https://ui.adsabs.harvard.edu/abs/2018ApJ...860...79T/)
 - `iqcor = 1` : [Tazaki et al. 2016, ApJ, 823, 70](https://ui.adsabs.harvard.edu/abs/2016ApJ...823...70T)
 - `iqcor = 2` : [Berry & Percival 1986, AcOpt, 33, 577](https://ui.adsabs.harvard.edu/abs/1986AcOpt..33..577B)
 - `iqcor = 3` : [Botet et al. 1995, JPhA, 28, 297](https://ui.adsabs.harvard.edu/abs/1995JPhA...28..297B)


# How to use it? 

The user must specify following input parameters in either `call.f90` or `call2.f90`:

- `df` : Fractal dimension (1 ≦ df ≦ 3)
- `k0` : Fractal prefactor
- `PN` : Number of monomers (1 ≦ PN)
- `R0` : Monomer radius (micron)
- `lmd`: wavelength (micron)
- `ref`: complex refractive index
- `nang`: Number of angle mesh (=91) 

In addition, the user also needs to specify following four options:

- Light scattering solver   
  `iqsca = 1` : Rayleigh-Gans-Debye theory  
  `iqsca = 2` : Mean field theory  
  `iqsca = 3` : Modified mean field theory  
- The two-point correlation function of fractal aggregates  
  `iqcor = 1` : The Gaussian cut-off model  
  `iqcor = 2` : The exponential cut-off model  
  `iqcor = 3` : The fractal dimension cut-off model  
- Geometric cross section of fractal aggregates (needed only when `iqsca=3`)  
  `iqgeo = 1` : The characteristic area   
  `iqgeo = 2` : Empirical formula by [Okuzumi et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...707.1247O)  
- standard output  
  `iquiet = 0` : show standard output  
  `iquiet = 1` : suppress standard output (including warnings)  
	
I recommend following set of options: `iqsca=3`,`iqcor=1`,`iqgeo=2` (default).  

To run the code, first, you make the codes by
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

- `iqsca = 1`   
  All outputs would be physically reasonable for the phase shift <~ 1.   
 
- `iqsca = 2`  
 The extinction cross section would be calculated without limitation.  
 However, the other outputs would be reliable for the phase shift <~ 1.  

- `iqsca = 3`  
  The extinction cross section would be calculated without limitation. Scattering and absorption cross sections could be calculated for the phase shift > 1, however too large phase shift may cause some problem. The asymmetry parameter and the sattering matrix elements would be reliable for the phase shift <~ 1.  

 For the two-point correlation function, I suggest to use `iqcor=1` (Gaussian type) because it is numerically stable and is likely to reproduce optical properties of fractal aggregates. For safety, the code always outputs a value of the phase shift, and therefore, the user can check the value of the phase shift. The code also produces a warning when the phase shift is above unity. 

# Acknowledgement 

`opacfractal` contains following subroutines provided by other authors.
I thank the authors for the availability of the subroutines.  

A list of subroutines by other authors:  

`lorenz_mie`  : Modified version of BHMIE code written by Bruce T. Draine  
`renorm_mie`  : Modified version of BHMIE code written by Bruce T. Draine  
`gamma`       : Zhang, S. and Jin, J. (1996) "Computation of Special Functions.  
`chgm`        : Zhang, S. and Jin, J. (1996) "Computation of Special Functions.  
`lpmns`       : Zhang, S. and Jin, J. (1996) "Computation of Special Functions.  
`lpn`         : Zhang, S. and Jin, J. (1996) "Computation of Special Functions.  
`ludcmp`      : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"    
`lubksb`      : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"  
`gauleg`      : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"  


# History

*Version 3.0*
- Initial release 
