# Introduction to opacfractal

This package allows to compute optical properties of fractal dust aggregates 
by means of a statistical distribution model of monomers.

# Terms of use

`opacfractal` is distributed under the [MITlicense](https://opensource.org/licenses/MIT) and can be used, changed
and redistributed freely. If you use this package to publish papers, please cite the following paper

> Ryo Tazaki and Hidekazu Tanaka  
> *Light Scattering by Fractal Dust Aggregates. II. Opacity and Asymmetry parameter*  
> The Astrophysical Journal, 860:79 (17pp), 2018

# Examples 

The user must specify following input parameters in `call.f90` or `call2.f90`

- `df` : Fractal dimension (1 ≦ df ≦ 3)
- `k0` : Fractal prefactor
- `PN` : Number of monomers (1 ≦ PN)
- `R0` : Monomer radius (micron)
- `lmd`: wavelength (micron)
- `ref`: complex refractive index
- `nang`: Number of angle mesh (=91) 

In addition, the user also needs to specify following three options

- Light scattering solver   
  `iqsca=1` : Rayleigh-Gans-Debye theory  
	`iqsca=2` : Mean field theory  
	`iqsca=3` : Modified mean field theory  
- The two-point correlation function of fractal aggregates  
  `iqcor=1` : The Gaussian cut-off model  
	`iqcor=2` : The exponential cut-off model  
	`iqcor=3` : The fractal dimension cut-off model  
- Geometric cross section of fractal aggregates  
	`iqgeo=1` : The characteristic area   
	`iqgeo=2` : Empitical formula by Okuzumi et al. (2009)  
- standard output  
	`iquiet=0` : show standard output  
	`iquiet=1` : suppress standard output (inclusing warnings)  
	
I recommend following set of options: `iqsca=3`,`iqcor=1`,`iqgeo=2` (default).  

To run the code, first, you make the codes by
```
Make
```
This will create an executable file `results.x`. Then, perform
```
./fracsca.x
```
As a result, the output file `out_smat.dat` or `out_opc.dat` is created. 

# History

*Version 3.0*
- Initial release
