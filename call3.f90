!--------------------------------------------------------------------------------
!
!   optical properties averaged over power-law aggregate size distribution
!
!--------------------------------------------------------------------------------
program opacfractal_distribution_average
use types; use const
use omp_lib
implicit none
real(kind=dp),parameter :: mic2cm = 1.0e-4_dp    
real(kind=dp),parameter :: cm2mic = 1.0e4_dp     
integer :: i,j,ia,na,nang,nwl,unit
integer :: iqsca,iqcor,iqgeo,iquiet
real(kind=dp) :: r0,rmax,rmin,p,pn,k0,df,lmd,Cext,Csca,Cabsp,gsca,dphi
real(kind=dp) :: lnrmin,lnrmax,dlna,rhod,wn,dang,tot,tot2,mu1,mu2,dmu
real(kind=dp),allocatable,dimension(:)   :: kabs,ksca,g,ang
real(kind=dp),allocatable,dimension(:)   :: wl,refre,refim,kabsave,kscaave,gsave,ri,nd,mass
real(kind=dp),allocatable,dimension(:,:) :: Z11,Z12,Z22,Z33,Z34,Z44,Smat
real(kind=dp),allocatable,dimension(:,:) :: Z11ave,Z12ave,Z22ave,Z33ave,Z34ave,Z44ave
complex(kind=dp)::refrel

!--------------------------------------------------------------------------------
!       Input parameters
!--------------------------------------------------------------------------------
k0     = 1.04_dp            ! fractal pre-factor
df     = 1.90_dp            ! fractal dimension
r0     = 0.01_dp * mic2cm   ! monomer radius                        [cm]
rmin   = 0.10_dp * mic2cm   ! minimum volume-equivalent radius      [cm]
rmax   = 3.00_dp * mic2cm   ! maximum volume-equivalent radius      [cm]
p      = 3.50_dp            ! power-law index of size distribution
rhod   = 3.30_dp            ! material density                      [g/cm^3]
nang   = 91                 ! number of angle grids 
nwl    = 100                ! number of wavelength in optical constant file
na     = 15                 ! number of grain radius bin
!--------------------------------------------------------------------------------
iqsca  = 3                  
iqcor  = 1                  
iqgeo  = 3                  
iquiet = 1

!--------------------------------------------------------------------------------
!       write message
!--------------------------------------------------------------------------------
call out_headerinfo(6,r0*cm2mic,rmin*cm2mic,rmax*cm2mic,p,k0,df,iqsca,iqcor,iqgeo,rhod)

!--------------------------------------------------------------------------------
!       allocate arrays
!--------------------------------------------------------------------------------

allocate(wl(1:nwl),refre(1:nwl),refim(1:nwl),kabsave(1:nwl),kscaave(1:nwl),gsave(1:nwl))
allocate(ri(1:na),nd(1:na),mass(1:na),kabs(1:na),ksca(1:na),g(1:na))
allocate(ang(1:2*nang-1),Smat(1:4,1:2*nang-1))
allocate(Z11(1:na,1:2*nang-1),Z12(1:na,1:2*nang-1),Z22(1:na,1:2*nang-1),&
         Z33(1:na,1:2*nang-1),Z34(1:na,1:2*nang-1),Z44(1:na,1:2*nang-1))
allocate(Z11ave(1:nwl,1:2*nang-1),Z12ave(1:nwl,1:2*nang-1),&
         Z22ave(1:nwl,1:2*nang-1),Z33ave(1:nwl,1:2*nang-1),&
         Z34ave(1:nwl,1:2*nang-1),Z44ave(1:nwl,1:2*nang-1))

!--------------------------------------------------------------------------------
!       Load optical constant
!--------------------------------------------------------------------------------
open(50,file="astrosil.nk",status="old")
do i=1,nwl
    read(50,*) wl(i),refre(i),refim(i)
end do
close(50)
wl = wl * mic2cm        

!--------------------------------------------------------------------------------
!       make size distribution
!--------------------------------------------------------------------------------
lnrmin = log(rmin)
lnrmax = log(rmax)
dlna   = (lnrmax-lnrmin) / real(na-1,kind=dp)
do ia=1,na
    ri(ia)   = exp(lnrmin+real(ia-1,kind=dp)*dlna)
    PN       = (ri(ia)/r0)**3.0
    mass(ia) = 4.0_dp * pi * r0 ** 3.0 * rhod * PN / 3.0_dp
    nd(ia)   = ri(ia)**(1.0-p)      ! logarithmic bin : n(a)dlna
end do
! Now, dlna = constant, so it is going to be omitted hereafter as it cancels out.
tot = 0.0_dp
do ia=1,na-1
    tot = tot + 0.5_dp * (mass(ia)*nd(ia)+mass(ia+1)*nd(ia+1))
enddo
nd = nd / tot

!--------------------------------------------------------------------------------
!       make angle mesh
!--------------------------------------------------------------------------------
dang  = halfpi/real(nang-1,kind=dp)
do j=1,2*nang-1
    ang(j)=dang*real(j-1,kind=dp)*r2d
end do

!--------------------------------------------------------------------------------
!       perform distribution average of optical properties
!--------------------------------------------------------------------------------
Z11ave  = 0.0_dp
Z12ave  = 0.0_dp 
Z22ave  = 0.0_dp 
Z33ave  = 0.0_dp 
Z34ave  = 0.0_dp 
Z44ave  = 0.0_dp 
kabsave = 0.0_dp
kscaave = 0.0_dp
gsave   = 0.0_dp

do i=1,nwl
    
    print *, '... wavelength',i,'/',nwl

    Z11 = 0.0_dp
    Z12 = 0.0_dp
    Z22 = 0.0_dp
    Z33 = 0.0_dp
    Z34 = 0.0_dp
    Z44 = 0.0_dp
    refrel = cmplx(refre(i),refim(i))
    wn     = twopi/wl(i)
    !
    ! storing opacities and scattering matrix at each aggregate radius
    !
    !$OMP parallel do schedule (static,1)                   &
    !$OMP default(none)                                     &
    !$OMP shared(refrel,wn,ri,wl,mass,i,nang,na)            &
    !$OMP shared(kabs,ksca,g,Z11,Z12,Z22,Z33,Z34,Z44)       &
    !$OMP shared(lmd,r0,Df,k0,iqsca,iqcor,iqgeo,iquiet)     &
    !$OMP private(PN,Cext,Csca,Cabsp,Gsca,Smat,dphi)
    do ia=1,na

        PN = (ri(ia)/r0)**3.0
        call meanscatt(wl(i)*cm2mic,r0*cm2mic,PN,Df,k0,refrel,iqsca,iqcor,iqgeo,iquiet,nang,&
                        Cext,Csca,Cabsp,Gsca,Smat,dphi)

        kabs(ia) = Cabsp * mic2cm * mic2cm / mass(ia)
        ksca(ia) = Csca  * mic2cm * mic2cm / mass(ia)
        g(ia)    = Gsca
        do j=1,2*nang-1
            Z11(ia,j) = Smat(1,j) / mass(ia) / wn / wn 
            Z12(ia,j) = Smat(2,j) / mass(ia) / wn / wn 
            Z22(ia,j) = Z11(ia,j)
            Z33(ia,j) = Smat(3,j) / mass(ia) / wn / wn 
            Z34(ia,j) = Smat(4,j) / mass(ia) / wn / wn 
            Z44(ia,j) = Z33(ia,j)
        end do

    end do
    !$OMP end parallel do
    !
    ! opacities and scattering matrix averaging over the size distribution
    !
    do ia=1,na-1
        kabsave(i) = kabsave(i) + 0.5_dp * (kabs(ia)*mass(ia)*nd(ia)+&
               kabs(ia+1)*mass(ia+1)*nd(ia+1))
        kscaave(i) = kscaave(i) + 0.5_dp * (ksca(ia)*mass(ia)*nd(ia)+&
               ksca(ia+1)*mass(ia+1)*nd(ia+1)) 
        do j=1,2*nang-1
            Z11ave(i,j) = Z11ave(i,j) + 0.5_dp * (Z11(ia,j)*mass(ia)*nd(ia)+&
                             Z11(ia+1,j)*mass(ia+1)*nd(ia+1))  
            Z12ave(i,j) = Z12ave(i,j) + 0.5_dp * (Z12(ia,j)*mass(ia)*nd(ia)+&
                             Z12(ia+1,j)*mass(ia+1)*nd(ia+1))   
            Z22ave(i,j) = Z22ave(i,j) + 0.5_dp * (Z22(ia,j)*mass(ia)*nd(ia)+&
                             Z22(ia+1,j)*mass(ia+1)*nd(ia+1))      
            Z33ave(i,j) = Z33ave(i,j) + 0.5_dp * (Z33(ia,j)*mass(ia)*nd(ia)+&
                             Z33(ia+1,j)*mass(ia+1)*nd(ia+1))      
            Z34ave(i,j) = Z34ave(i,j) + 0.5_dp * (Z34(ia,j)*mass(ia)*nd(ia)+&
                             Z34(ia+1,j)*mass(ia+1)*nd(ia+1))      
            Z44ave(i,j) = Z44ave(i,j) + 0.5_dp * (Z44(ia,j)*mass(ia)*nd(ia)+&
                             Z44(ia+1,j)*mass(ia+1)*nd(ia+1))        
        enddo
    enddo
    !
    ! compute asymmetry parameter by integrating Z11
    !   in mu-space (RADMC-3D-like weightning)
    !
    tot  = 0.0_dp
    tot2 = 0.0_dp
    do j=1,2*nang-2 
        mu1 = cos(ang(j)*d2r); mu2 = cos(ang(j+1)*d2r); dmu = twopi * (mu1-mu2)
        tot  = tot  + 0.5_dp * (Z11ave(i,j)+Z11ave(i,j+1)) * dmu
        tot2 = tot2 + 0.5_dp * (Z11ave(i,j)+Z11ave(i,j+1)) * dmu * 0.5_dp * (mu1 + mu2)
    enddo
    gsave(i) = tot2 / tot

enddo ! end wavelength loop

!--------------------------------------------------------------------------------
!       output results in format of RADMC-3D input file
!--------------------------------------------------------------------------------
print *, 'dustkapscatmat.inp'
unit = 10
open(unit,file='dustkapscatmat.inp',status='unknown')
call out_headerinfo(unit,r0*cm2mic,rmin*cm2mic,rmax*cm2mic,p,k0,df,iqsca,iqcor,iqgeo,rhod)
call out_radmc3d_format(unit,nwl,2*nang-1,wl,kabsave,kscaave,gsave,ang,&
            Z11ave,Z12ave,Z22ave,Z33ave,Z34ave,Z44ave)
close(unit)

deallocate(wl,refre,refim,ri,nd,mass)
deallocate(Smat,kabs,ksca,g,Z11,Z12,Z22,Z33,Z34,Z44)
deallocate(ang,kabsave,kscaave,gsave,Z11ave,Z12ave,Z22ave,Z33ave,Z34ave,Z44ave)

stop
end program opacfractal_distribution_average

!--------------------------------------------------------------------------------
!       output header information
!--------------------------------------------------------------------------------
subroutine out_headerinfo(unit,r0,rmin,rmax,p,k0,df,iqsca,iqcor,iqgeo,rhod)
use types
implicit none
integer::unit,iqsca,iqcor,iqgeo
real(kind=dp):: df,k0,r0,rmin,rmax,p,rhod

write(unit,'("# ================================================")')
write(unit,'("#              OPACFRACTAL v3.0                  ")') 
write(unit,'("# ================================================")')
write(unit,'("# iqsca = ",I1," iqcor = ",I1," iqgeo = ",I1)') iqsca,iqcor,iqgeo
write(unit,'("# df           = ",1PE12.5)') df
write(unit,'("# k0           = ",1PE12.5)') k0
write(unit,'("# r0   [um]    = ",1PE12.5)') r0
write(unit,'("# rmin [um]    = ",1PE12.5)') rmin
write(unit,'("# rmax [um]    = ",1PE12.5)') rmax
write(unit,'("# powind       = ",1PE12.5)') p
write(unit,'("# rho (g/cm^3) = ",1PE12.5)') rhod

return
end subroutine out_headerinfo

!--------------------------------------------------------------------------------
!
!       output RADMC-3D input file
!
!--------------------------------------------------------------------------------
subroutine out_radmc3d_format(unit,nwl,nang,wl,kabs,ksca,g,ang,Z11,Z12,Z22,Z33,Z34,Z44)
use types; use const
implicit none
integer::i,j,nwl,nang,unit
real(kind=dp),parameter :: cm2mic = 1.0e4_dp     
real(kind=dp),dimension(1:nwl)::wl,kabs,ksca,g
real(kind=dp),dimension(1:nwl,1:nang)::Z11,Z12,Z22,Z33,Z34,Z44
real(kind=dp)::ang(1:nang),tot,mu1,mu2,dmu,f

write(unit,*) 1       !ifotmat
write(unit,*) nwl     !number of wavelength
write(unit,*) nang    !number of anguler mesh
write(unit,*)
do i=1,nwl
    write(unit,'(1P4E18.8)') wl(i)*cm2mic,kabs(i),ksca(i),g(i)
end do
write(unit,*)
do j=1,nang
    write(unit,'(1PE18.8)') ang(j)
end do
write(unit,*)
do i=1,nwl
    tot = 0.0_dp
    ! mu-space integration of Z11 (as done in RADMC3D)
    do j=1,nang-1
        mu1 = cos(ang(j)*d2r); mu2 = cos(ang(j+1)*d2r); dmu=twopi*(mu1-mu2)
        tot = tot + 0.5_dp * ( Z11(i,j) + Z11(i,j+1) ) * dmu
    enddo
    ! Rescaling factor of Zij to match scattering opacity
    f = ksca(i) / tot
    do j=1,nang
        write(unit,'(1P6E18.8)') f*Z11(i,j),f*Z12(i,j),f*Z22(i,j),&
                                 f*Z33(i,j),f*Z34(i,j),f*Z44(i,j)
    end do        
end do

return
end subroutine out_radmc3d_format
