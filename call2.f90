!--------------------------------------------------------------------------------
!
! An example to call opacfractal.f90
!      
!--------------------------------------------------------------------------------
program call_opacfractal
use omp_lib
implicit none
integer, parameter                      :: dp = selected_real_kind(P=15)
real(kind=dp),parameter                 :: pi = 3.1415926535897932384_dp
real(kind=dp),parameter                 :: k0_bpca=0.3_dp
real(kind=dp),parameter                 :: k0_bcca=1.04_dp
real(kind=dp),parameter                 :: df_bpca=3.0_dp
real(kind=dp),parameter                 :: df_bcca=1.9_dp
integer                                 :: iqsca,iqcor,iqgeo,iquiet,i,j,nang,unit
complex(kind=dp)                        :: ref
real(kind=dp)                           :: lmd,R0,Rv,PN,df_min,df_max
real(kind=dp)                           :: Cext,Csca,Cabsp,Gsca,dphi,rhod,mass
real(kind=dp)                           :: wlmin,wlmax,dwl
real(kind=dp),allocatable,dimension(:,:):: Smat
integer                                 :: iwl
integer,parameter                       :: nwl=100
integer,parameter                       :: ndf=7
real(kind=dp),dimension(1:ndf)          :: DF,k0
real(kind=dp),dimension(1:nwl)          :: wl,re,im
real(kind=dp),dimension(1:ndf,0:5,1:nwl):: Q

!--------------------------------------------------------------------------------
!       SET INPUT PARAMETERS
!--------------------------------------------------------------------------------
iqsca  = 3                           ! A method to solve light scattering
iqcor  = 1                           ! correlation function
iqgeo  = 3                           ! geometric cross section
iquiet = 1                           ! stdout
R0     = 1.e-1_dp                    ! Monomer radius                    [um]
Rv     = 1.e1_dp                     ! Volume equivalent radius (micron) [um]
rhod   = 3.3_dp                      ! material density                  [g/cm^3]
PN     = (Rv/R0)**3.0_dp             ! Number of monomer in an aggregate
nang   = 90*100+1                    ! Angle mesh

!--------------------------------------------------------------------------------
!        AGGREGATE STRUCTURE MODEL 
!--------------------------------------------------------------------------------
df_min = 1.8_dp
df_max = 3.0_dp
do j=1,ndf
    DF(j) = df_min + real(j-1,kind=dp)*(df_max-df_min)/real(ndf-1,kind=dp)
    k0(j) = 0.716_dp * (1.0_dp - DF(j)) + sqrt(3.0_dp)
enddo
mass = 4.0_dp * pi * rv ** 3.0_dp * rhod / 3.0_dp

!--------------------------------------------------------------------------------
!        Wavelength mesh
!--------------------------------------------------------------------------------
open(16,file="astrosil.nk",status="old")
do iwl=1,nwl
    read(16,*) wl(iwl),re(iwl),im(iwl)
enddo
close(16)

!--------------------------------------------------------------------------------
!        RUN OPACFRACTAL
!--------------------------------------------------------------------------------
call out_headerinfo(6,r0,rv,iqsca,iqcor,iqgeo)
allocate(Smat(1:4,1:2*nang-1))
do j=1,ndf
    !$OMP parallel do schedule (static,1)                     &
    !$OMP default(none)                                       &
    !$OMP shared(wl,re,im,mass,Q)                             &
    !$OMP shared(j,R0,PN,Df,k0,iqsca,iqcor,iqgeo,iquiet,nang) &
    !$OMP private(lmd,ref,Cext,Csca,Cabsp,Gsca,Smat,dphi)   
    do iwl=1,nwl
        lmd=wl(iwl)
        ref=cmplx(re(iwl),im(iwl),kind=dp)
        write(*,'("df = ",f8.3," | k0 = ",f8.3," | wavelength (um) = ",1PE18.8)') df(j),k0(j),lmd
        call meanscatt(lmd,R0,PN,Df(j),k0(j),ref,iqsca,iqcor,iqgeo,iquiet,nang,&
                        Cext,Csca,Cabsp,Gsca,Smat,dphi)
        Q(j,1,iwl) = Cext  * 1.e4 / mass
        Q(j,2,iwl) = Csca  * 1.e4 / mass
        Q(j,3,iwl) = Cabsp * 1.e4 / mass
        Q(j,4,iwl) = gsca
        Q(j,5,iwl) = dphi
    enddo
    !$OMP end parallel do
enddo
deallocate(Smat)

!--------------------------------------------------------------------------------
!       OUTPUT 
!--------------------------------------------------------------------------------
unit=10
open(unit,file="out_opc.dat",status="unknown")
call out_headerinfo(unit,r0,rv,iqsca,iqcor,iqgeo)
do j=1,ndf
write(10,'("#",A17,7A18)') "df","k0","lambda (um)",&
    "kext (cm^2/g)","ksca (cm^2/g)","kabs (cm^2/g)","<cos>","phase shift"
do iwl=1,nwl 
    write(10,'(1P8E18.8)') df(j),k0(j),wl(iwl),(Q(j,i,iwl),i=1,5)
enddo
write(10,*)
enddo
close(unit)

print *, "finished!"
stop
end program call_opacfractal

!--------------------------------------------------------------------------------
!       output header information
!--------------------------------------------------------------------------------
subroutine out_headerinfo(unit,r0,rv,iqsca,iqcor,iqgeo)
use types
implicit none
integer::unit,iqsca,iqcor,iqgeo
real(kind=dp):: df,k0,r0,rv,p,rhod

write(unit,'("# ================================================")')
write(unit,'("#              OPACFRACTAL v3.0                  ")') 
write(unit,'("# ================================================")')
write(unit,'("# iqsca = ",I1," iqcor = ",I1," iqgeo = ",I1)') iqsca,iqcor,iqgeo
write(unit,'("# monomer radius (um)                   = ",1PE15.5)') r0
write(unit,'("# volume-equivalent radius (um)         = ",1PE15.5)') rv

return
end subroutine out_headerinfo

