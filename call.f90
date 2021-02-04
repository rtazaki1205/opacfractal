!--------------------------------------------------------------------------------
!
! An example to call opacfractal.f90
!
! The default parameter set is the same as the benchmark test used in Appendix A 
! in Tazaki & Tanaka (2018), which is based on RCCA model.
! 
!--------------------------------------------------------------------------------
program call_opacfractal
implicit none
integer, parameter                      :: dp = selected_real_kind(P=15)
real(kind=dp),parameter                 :: pi = 3.1415926535897932384_dp
integer                                 :: iqsca,iqcor,iqgeo,iquiet,i,nang,unit
complex(kind=dp)                        :: ref
real(kind=dp)                           :: lmd,R0,PN,DF,k0,ang,dang
real(kind=dp)                           :: Cext,Csca,Cabsp,Gsca,dphi
real(kind=dp),allocatable,dimension(:,:):: Smat
!--------------------------------------------------------------------------------
! SET INPUT PARAMETERS
!--------------------------------------------------------------------------------
iqsca  = 3                               ! a method to solve light scattering
iqcor  = 1                               ! correlation function
iqgeo  = 3                               ! geometric cross section
iquiet = 0                               ! stdout
R0     = 0.5_dp                          ! Monomer radius (micron)
PN     = 64.0_dp                         ! Number of monomers
nang   = 91                              ! Angle mesh
lmd    = 0.8_dp                          ! wavelength (micron)
k0     = 0.825_dp                        ! fractal prefactor
df     = 2.0_dp                          ! fractal dimension
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
unit = 10
DANG=0.5d0*PI/DBLE(NANG-1)
open(unit,file="out_smat.dat",status="unknown")
call out_headerinfo(unit,r0,lmd,PN,k0,df,ref,iqsca,iqcor,iqgeo)
write(unit,'("#",1PE18.8,A)') Cext," = Extinction cross-section (um^2)"
write(unit,'("#",1PE18.8,A)') Csca," = Scattering cross-section (um^2)"
write(unit,'("#",1PE18.8,A)') Cabsp," = Absorption cross-section (um^2)"
write(unit,'("#",1PE18.8,A)') gsca," = asymmetry parameter"
write(unit,'("#",A17,4A18)') "degree","S11","S12","S33","S34"
do i=1,2*nang-1
        ang = dang*dble(i-1)*(180.0/pi)
        write(unit,'(1P5E18.8)') ang,Smat(1,i),Smat(2,i),Smat(3,i),Smat(4,i)
enddo
deallocate(Smat)
!--------------------------------------------------------------------------------
print *, "finished!"

stop
end program call_opacfractal

!--------------------------------------------------------------------------------
!       output header information
!--------------------------------------------------------------------------------
subroutine out_headerinfo(unit,r0,lmd,PN,k0,df,ref,iqsca,iqcor,iqgeo)
use types
implicit none
integer::unit,iqsca,iqcor,iqgeo
real(kind=dp):: df,k0,r0,PN,rhod,lmd
complex(kind=dp)::ref

write(unit,'("# ================================================")')
write(unit,'("#              OPACFRACTAL v3.0                  ")') 
write(unit,'("# ================================================")')
write(unit,'("# iqsca = ",I1," iqcor = ",I1," iqgeo = ",I1)') iqsca,iqcor,iqgeo
write(unit,'("# df           = ",1PE12.5)') df
write(unit,'("# k0           = ",1PE12.5)') k0
write(unit,'("# r0   [um]    = ",1PE12.5)') r0
write(unit,'("# N            = ",1PE12.5)') PN
write(unit,'("# lambda [um]  = ",1PE12.5)') lmd
write(unit,'("# Re(m) = ",1PE12.5," Im(m) = ",1PE12.5)') real(ref),aimag(ref)

return
end subroutine out_headerinfo

