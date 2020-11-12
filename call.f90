!--------------------------------------------------------------------------------
!
! An example to call opacfractal.f90
!
! The default parameter set is the same as the benchmark test 
! used in Appendix A in Tazaki & Tanaka (2018).
!
!--------------------------------------------------------------------------------
program call_opacfractal
implicit none
integer, parameter                      :: dp = selected_real_kind(P=15)
real(kind=dp),parameter                 :: pi = 3.1415926535897932384_dp
integer                                 :: iqsca,iqcor,iqgeo,iquiet,i,nang
complex(kind=dp)                        :: ref
real(kind=dp)                           :: lmd,R0,PN,DF,k0
real(kind=dp)                           :: ang,dang
real(kind=dp)                           :: Cext,Csca,Cabsp,Gsca,dphi
real(kind=dp),allocatable,dimension(:,:):: Smat
!--------------------------------------------------------------------------------
! SET INPUT PARAMETERS
!--------------------------------------------------------------------------------
iqsca  = 3                               ! a method to solve light scattering
iqcor  = 1                               ! correlation function
iqgeo  = 1                               ! geometric cross section
iquiet = 0                               ! stdout
R0     = 0.5_dp                          ! Monomer radius (micron)
PN     = 64.0_dp                         ! Number of monomers
nang   = 91                              ! Angle mesh
lmd    = 0.8_dp                          ! wavelength (micron)
k0     = 0.825_dp                        ! fractal prefactor
df     = 3.0_dp                          ! fractal dimension
ref    = cmplx(1.4_dp,0.0001_dp,kind=dp) ! complex refractive index
!--------------------------------------------------------------------------------
! RUN OPACFRACTAL
!--------------------------------------------------------------------------------
allocate(Smat(1:4,1:2*nang-1))
call meanscatt(lmd,R0,PN,Df,k0,ref,iqsca,iqcor,iqgeo,iquiet,nang,&
                Cext,Csca,Cabsp,Gsca,Smat,dphi)
!--------------------------------------------------------------------------------
! OUTPUT 
!--------------------------------------------------------------------------------
DANG=0.5d0*PI/DBLE(NANG-1)
open(10,file="out_smat.dat",status="unknown")
write(10,1100) iqsca,   " = iqsca"
write(10,1100) iqcor,   " = iqcor"
write(10,1100) iqgeo,   " = iqgeo"
write(10,1000) lmd,     " = wavelength (um)"
write(10,1000) R0,      " = Monomer radius (um)"
write(10,1000) PN,      " = Number of monomers"
write(10,1000) DF,      " = Fractal dimension"
write(10,1000) k0,      " = Fractal prefactor"
write(10,1000) real(ref)," = Real part of refractive index"
write(10,1000) aimag(ref)," = Imag part of refractive index"
write(10,*) "---------------------------------"
write(10,1000) Cext," = Extinction cross-section (um^2)"
write(10,1000) Csca," = Scattering cross-section (um^2)"
write(10,1000) Cabsp," = Absorption cross-section (um^2)"
write(10,1000) gsca," = asymmetry parameter"
write(10,2000) "degree","S11","S12","S33","S34"
do i=1,2*nang-1
        ang = dang*dble(i-1)*(180.0/pi)
        write(10,2100) ang,Smat(1,i),Smat(2,i),Smat(3,i),Smat(4,i)
enddo
deallocate(Smat)
!--------------------------------------------------------------------------------
print *, "finished!"
1000 format('#',1PE15.5,A)
1100 format('#',I15,A)
2000 format('#',6A15)
2100 format(' ',1P5E15.5)
stop
end program call_opacfractal
