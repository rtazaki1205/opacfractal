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
integer                                 :: iqsca,iqcor,iqgeo,iquiet,i,j,nang
complex(kind=dp)                        :: ref
real(kind=dp)                           :: lmd,R0,Rv,PN,df_min,df_max
real(kind=dp)                           :: Cext,Csca,Cabsp,Gsca,dphi
real(kind=dp)                           :: wlmin,wlmax,dwl
real(kind=dp),allocatable,dimension(:,:):: Smat
integer                                 :: iwl
integer,parameter                       :: nwl=100
integer,parameter                       :: ndf=16
real(kind=dp),dimension(1:ndf)          :: DF,k0
real(kind=dp),dimension(1:nwl)          :: wl,re,im
real(kind=dp),dimension(1:ndf,0:5,1:nwl):: CS
real(kind=dp)::Gratio,tmp
!--------------------------------------------------------------------------------
! SET INPUT PARAMETERS
!--------------------------------------------------------------------------------
iqsca  = 3                               ! A method to solve light scattering
iqcor  = 1                               ! correlation function
iqgeo  = 3                               ! geometric cross section
iquiet = 1                               ! stdout
R0     = 1.e-1_dp                        ! Monomer radius (micron)
Rv     = 3.e1_dp                         ! Volume equivalent radius (micron)
PN     = (Rv/R0)**3.0_dp                 ! Number of monomer in an aggregate
nang   = 90*100+1                        ! Angle mesh
!--------------------------------------------------------------------------------
! AGGREGATE STRUCTURE MODEL 
!--------------------------------------------------------------------------------
df_min = 1.5_dp
df_max = 3.0_dp
do j=1,ndf
        DF(j) = df_min + real(j-1,kind=dp)*(df_max-df_min)/real(ndf-1,kind=dp)
        k0(j) = (k0_bpca-k0_bcca)/(df_bpca-df_bcca)*(Df(j)-df_bpca)+k0_bpca
enddo
!--------------------------------------------------------------------------------
! Wavelength mesh
!--------------------------------------------------------------------------------
open(16,file="astrosil.nk",status="old")
do iwl=1,nwl
        read(16,*) wl(iwl),re(iwl),im(iwl)
enddo
close(16)
!--------------------------------------------------------------------------------
! RUN OPACFRACTAL
!--------------------------------------------------------------------------------
allocate(Smat(1:4,1:2*nang-1))
do j=1,ndf
        !$OMP parallel do schedule (static,1)                     &
        !$OMP default(none)                                       &
        !$OMP shared(wl,re,im,CS)                                 &
        !$OMP shared(j,R0,PN,Df,k0,iqsca,iqcor,iqgeo,iquiet,nang) &
        !$OMP private(lmd,ref,Cext,Csca,Cabsp,Gsca,Smat,dphi,tmp)   
        do iwl=1,nwl
                lmd=wl(iwl)
                ref=cmplx(re(iwl),im(iwl),kind=dp)
                write(*,fmt='(A5,f3.1,A15,1PE15.5,A15,I3)') "Df = ",df(j),&
                        "wavel (mic) = ",lmd,&
                        ": thread = ",omp_get_thread_num()
                call meanscatt(lmd,R0,PN,Df(j),k0(j),ref,iqsca,iqcor,iqgeo,iquiet,nang,&
                                Cext,Csca,Cabsp,Gsca,Smat,dphi)
                CS(j,1,iwl) = Cext
                CS(j,2,iwl) = Csca
                CS(j,3,iwl) = Cabsp
                CS(j,4,iwl) = gsca
                CS(j,5,iwl) = dphi
        enddo
        !$OMP end parallel do
enddo
deallocate(Smat)
!--------------------------------------------------------------------------------
! OUTPUT 
!--------------------------------------------------------------------------------
open(10,file="out_opc.dat",status="unknown")
do j=1,ndf
write(10,1100) iqsca,   " = iqsca"
write(10,1100) iqcor,   " = iqcor"
write(10,1100) iqgeo,   " = iqgeo"
write(10,1000) R0,      " = Monomer radius (um)"
write(10,1000) PN,      " = Number of monomers"
write(10,1000) real(ref), " = Re(m)"
write(10,1000) aimag(ref)," = Im(m)"
write(10,1000) DF(j),   " = Fractal dimension"
write(10,1000) k0(j),   " = Fractal prefactor"
write(10,*) "---------------------------------"
write(10,2000) "wl (um)","Cext (um^2)","Csca (um^2)","Cabs (um^2)","<cos>","dphi"
do iwl=1,nwl 
        write(10,2100) wl(iwl),(CS(j,i,iwl),i=1,5)
enddo
write(10,*)
enddo
!--------------------------------------------------------------------------------
print *, "finished!"
1000 format('#',1PE15.5,A)
1100 format('#',I15,A)
2000 format('#',6A15)
2100 format(' ',1P6E15.5)
stop
end program call_opacfractal
