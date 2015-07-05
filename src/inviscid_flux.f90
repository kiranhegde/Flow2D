! =======================================================
subroutine fluxinv_int(nx,ny,idir,r1,r2,r3, & 
            r4,aml,amr,qleft,qright,qfleft,qfright,fintr)
! =======================================================
!calculates the interface flux (fintr) using ausm-dv flux
!function, given qleft,qright,qfleft,qright and normals
!aml and amr
!
!variables r1-r4 are dummy in the present implementation


!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!nx    = same as n for this code
!ny    = same as n for this code
!idir  = 1 for xi direction, 2 for eta direction
!r1,r2,r3,r4 = distances of cell centres on the left and right of
!        the cell interface under consideration
!aml,amr = left and right face normals, amr is for cell face
!          under consideration, the code always calculates fluxes
!          for the face on the right of the cell
!qleft,qright = left the right cell centre values of the face
!               under consideration
!fintr = output, flux at the interface
!===============================================================
use dims
use inf
use pri

implicit none
integer(kind=i2):: nx,ny,idir,i,j
real(kind=dp):: r1,r2,r3,r4,aml(ndim),amr(ndim),fintr(ndim+2)
real(kind=dp):: qleft(ndim+2),qright(ndim+2)
real(kind=dp):: qfleft(ndim+2),qfright(ndim+2)
real(kind=dp):: qneg(ndim+2),qpos(ndim+2)

! calculating the muscl outputs q+ and q- at interface
! m+1/2

call muscl(qleft,qfleft,qright,qfright,aml,amr,qpos,qneg,r1,r2,r3,r4)

! qpos(:) = qright(:)
! qneg(:) = qleft(:)

call flux_ausmPlus_Up(idir,qpos,qneg,amr,fintr)

end subroutine fluxinv_int 
!==============================================================================
!
!
!
!
!==============================================================================
real(kind=dp) function m1p(mach)
use dims
implicit none

real(kind=dp):: mach

m1p=0.50*(mach+dabs(mach))

end function m1p
!====================================
real(kind=dp) function m1m(mach)
use dims
implicit none

real(kind=dp):: mach

m1m=0.5d0*(mach-dabs(mach))

end function m1m
!====================================
real(kind=dp) function m2p(mach)
use dims
implicit none

real(kind=dp):: mach

m2p=0.25d0*(mach+1)*(mach+1)

end function m2p
!====================================
real(kind=dp) function m2m(mach)
use dims
implicit none

real(kind=dp):: mach

m2m=-0.25d0*(mach-1)*(mach-1)

end function m2m
!====================================
real(kind=dp) function m4p(mach,beta)
use dims
implicit none

real(kind=dp):: mach,beta,m2m,m2p

if(dabs(mach)>=1) then
   m4p=0.5d0*(mach+dabs(mach))
else
   m4p=m2p(mach)*(1-16*beta*m2m(mach))
endif

end function m4p
!====================================
real(kind=dp) function m4m(mach,beta)
use dims
implicit none

real(kind=dp):: mach,beta,m2m,m2p

if(dabs(mach)>=1) then
   m4m=0.5d0*(mach-dabs(mach))
else
   m4m=m2m(mach)*(1+16*beta*m2p(mach))
endif

end function m4m
!====================================
real(kind=dp) function p5p(mach,alpha)
use dims
implicit none

real(kind=dp):: mach,alpha,m2m,m2p

if(dabs(mach)>=1) then
   p5p=0.5d0*(mach+dabs(mach))/mach
else
   p5p=m2p(mach)*( (2-mach)-16*alpha*mach*m2m(mach))
endif

end function p5p
!====================================
real(kind=dp) function p5m(mach,alpha)
use dims
implicit none

real(kind=dp):: mach,alpha,m2m,m2p

if(dabs(mach)>=1) then
   p5m=0.5d0*(mach-dabs(mach))/mach
else
   p5m=m2m(mach)*( (-2-mach)+16*alpha*mach*m2p(mach))
endif

end function p5m
!==============================================================================
!
!
!
!
!==============================================================================
subroutine flux_ausmPlus_Up(idir,qpos,qneg,nxyz,flux)
!==============================================================================
!#     computes total convective flux across a face using van leer
!#  flux vector splitting method given left and right conserved states
!#
!#     nxyz  - components of face normal vector. 
!#     ql,qr - left & right vector of conserved variables.
!#     flux  - vector of flux across face.
!------------------------------------------------------------------------------
use dims
use inf
use pri
implicit none
!------------------------------------------------------------------------------
integer(kind=i2):: i,j,idir,bk 
real(kind=dp):: nxyz(ndim)
real(kind=dp):: qpos(ndim+2), qneg(ndim+2)
real(kind=dp):: flux(ndim+2)
real(kind=dp):: nxd,nyd
real(kind=dp):: nx, ny, area
real(kind=dp):: rol, ul, vl, pl, cl, unorml, ml, hl
real(kind=dp):: ror, ur, vr, pr, cr, unormr, mr, hr
real(kind=dp):: m12, mlp, mrm, aml, amr, mrp, mlm,p5p,p5m,m4p,m4m,alpha,beta
real(kind=dp):: p12, plp, prm,dmm, dpp, ro12, mp, m1m, m1p, m2p, m2m
real(kind=dp):: fluxL(ndim+2),fluxR(ndim+2),fluxP(ndim+2)
real(kind=dp):: mco, mref, mo, mbar, mbar2, mstar, mstar2, fa, c12, mass_p, mass_m, mass12, dissi
real(kind=dp):: sigma, kp, ku, pu
real(kind=dp):: aL_star,aR_star,aL_hat,aR_hat,atilR,atilL,astarL,astarR
real(kind=dp):: mass,fluxN,fluxT 
real(kind=dp):: term1,term2 
!------------------------------------------------------------------------------
bk=1
fa=0.d0
sigma=1.00
ku=0.75
kp=0.25

nx=nxyz(1)
ny=nxyz(2)
area = dsqrt(nx*nx + ny*ny)
nx = nx/area
ny = ny/area

!.......left state
call q_conv(qneg)
rol = rho   
ul  = u         
vl  = v         
pl  = p 
hl  = h 
cl  = a
unorml = ul*nx + vl*ny
!.......right state
call q_conv(qpos)
ror = rho   
ur  = u         
vr  = v         
pr  = p 
hr  = h 
cr  = a
unormr = ur*nx + vr*ny

!aL_star = cl
!aR_star = cr
aL_star=dsqrt(2.0*(gamma-1.0)*hl/(gamma+1.0))
aR_star=dsqrt(2.0*(gamma-1.0)*hr/(gamma+1.0))
aL_hat = aL_star*aL_star/dmax1(aL_star,unorml)
aR_hat = aR_star*aR_star/dmax1(aR_star,-1.0*unormr)

c12=dmin1(aL_hat,aR_hat)
!c12=0.5d0*(cl+cr)

ml = unorml/c12
mr = unormr/c12
mbar2  = 0.50*(ml*ml+mr*mr)

! Cut-off Mach number
!mref = dmax1(0.32,0.5*m_inf)
!mref=0.25*(dmin1(m_inf,1.0d0))
!mref=dmax1(m_inf,dsqrt(0.5*(ml+mr)))
mref=dmax1(m_inf,dsqrt(mbar2))
mref =dmin1(mref,1.0)
!mref=m_inf

if(mbar2 >=1.0) then
fa = 1.0
else
Mo = dsqrt(dmin1(1.0,dmax1(mbar2,mref*mref)))
!if(mstar2<0.0.or.mstar2>1.0) then 
!   print *,'mstar2:',mstar2
!   stop
!endif
fa = Mo*(2.0-Mo)
endif




beta=1.0/8.0
alpha=3.0/16.0*(-4.0+5.0*fa*fa)
plp=p5p(ml,alpha)
prm=p5m(mr,alpha)
mlp=m4p(ml,beta)
mrm=m4m(mr,beta)

ro12= 0.5d0*(rol+ror)

term1=(pr-pl)/(ro12*c12*c12)
mp = (kp/fa)*dmax1((1.0-sigma*mbar2),0.0)*term1
!if (iter2<=2.and.it==201.and.jt==98) then
if (isnan(mp)) then
   print*,iter2
   print '(2(i3,1x),4(f15.8,1x))',it,jt,pr-pl,fa,c12,ro12
   print '(1x,"Q-",6(f15.8,2x))',rol,ul,vl,hl,pl,cl 
   print '(1x,"Q+",6(f15.8,2x))',ror,ur,vr,hr,pr,cr 
   call get_cell(bk)
   stop
endif  


m12 = mlp + mrm - mp

term2=c12*(unormr-unorml)
pu = ku*plp*prm*2.0*ro12*fa*term2
p12 = plp*pl + prm*pr - pu

mass=c12*m12*ror
if(m12>0.0) mass=c12*m12*rol

!-----------------------------> Compute total flux <--------------------------!
flux(:)  = 0.0d0
fluxP(:) = 0.0d0
!------>  Pressure Flux
fluxP(3)=nx*P12
fluxP(4)=ny*P12

flux(1) = hr
flux(2) = 1.0d0
flux(3) = ur
flux(4) = vr

if(m12>0.0) then
   flux(1) = hl
   flux(2) = 1.0d0
   flux(3) = ul
   flux(4) = vl
endif   

!------>  total  flux
flux(:) =area*(mass*flux(:)+fluxP(:))
!
!
!!------->  Mass Flux
!mass = 0.5d0*(m12*c12*(rol+ror)-dabs(m12*c12)*(ror-rol))
!flux(1) = 0.5d0*(mass*(hr+hl)-dabs(mass)*(hr-hl)) 
!flux(2) = mass 
!
!!------>  total  flux
!if(idir==1)  then
!   fluxN = 0.5d0*(mass*(ur+ul)-dabs(mass)*(ur-ul)) 
!   fluxT = 0.5d0*(mass*(vr+vl)-dabs(mass)*(vr-vl)) 
!   flux(3) = fluxN + p12*nx
!   flux(4) = fluxT + p12*ny
!elseif(idir==2) then
!   fluxT = 0.5d0*(mass*(ur+ul)-dabs(mass)*(ur-ul)) 
!   fluxN = 0.5d0*(mass*(vr+vl)-dabs(mass)*(vr-vl)) 
!   flux(3) = fluxT + p12*nx
!   flux(4) = fluxN + p12*ny
!endif
!
!flux(:) = area*flux(:) 

end subroutine flux_ausmPlus_Up 
!==============================================================================


subroutine get_cell(bk)
use flux
use st_end
implicit none
integer(kind=i2):: i,j,k,bk

k=1
call start_end(bk)
!print *,is,ie
!print *,js,je
do j=jt,jt+1
do i=it,it+1
   print *,i,j
   write(22,100)blk(bk)%mesh(i,j,k)%x,blk(bk)%mesh(i,j,k)%y
enddo
   write(22,*)
enddo
do i=it,it+1
do j=jt,jt+1
   write(22,100)blk(bk)%mesh(i,j,k)%x,blk(bk)%mesh(i,j,k)%y
enddo
   write(22,*)
enddo
100 format(1x,2(f15.8,1x))   
end subroutine get_cell

