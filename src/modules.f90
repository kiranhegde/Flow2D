module data_kinds
implicit none
integer, parameter :: i1=selected_int_kind(2)
integer, parameter :: i2=selected_int_kind(4)
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(18)
integer, parameter :: sp=selected_real_kind(6,37)
integer, parameter :: dp=selected_real_kind(15,307)
integer, parameter :: qp=selected_real_kind(31,307)

end module data_kinds
!---------------------------------------------------------
module dims
use data_kinds
implicit none

integer (kind=i2),parameter :: ndim = 2
integer (kind=i2):: iter2,error=0,nblk=1,it,jt
integer (kind=i2):: bcc(0:8)
integer (kind=i2),parameter:: one =1 , two=2, three=3, four=4 
real(kind=dp):: pi = 4.0d0*datan(1.d0)

end module dims
!---------------------------------------------------------
module vis 
use data_kinds
implicit none

integer (kind=i2):: ieuler,iwall,turb
real(kind=dp):: t_w,cfl_fact,pr_l,re_l,pr_t

type turbulence
     real(kind=dp) :: yy,rh,uu,vo,vi,tauw,tv,bc
end type turbulence

type(turbulence),allocatable,dimension(:,:,:)::Blmx

end module vis 
!---------------------------------------------------------
module inf
use data_kinds
implicit none

integer (kind=i2):: initer 
real(kind=dp):: mu_inf,m_inf,p_inf,pinf,t_inf,rho_inf,gamma,gas_const,pout,aoa
real(kind=dp):: qinf(4),ent_inf,a_inf,Pdyna
real(kind=dp):: M_max,M_min 

end module inf
!---------------------------------------------------------
module pri 
use data_kinds
implicit none

real(kind=dp):: rho,u,v,e,q,p,h,a,t,znd

contains 

subroutine q_conv(qtemp)
use dims
use inf
!converts conserved variable to primitive variables

implicit none
real(kind=dp):: qtemp(ndim+2)

rho = qtemp(2)
u   = qtemp(3)/qtemp(2)
v   = qtemp(4)/qtemp(2)
e   = qtemp(1)

q   = 0.5d0*(u*u + v*v)
p   = (gamma-1.d0)*(e - rho*q)
h   = (e+p)/rho
a   = dsqrt(gamma*p/rho)
znd = 1.d0/(gamma*m_inf*m_inf)
t   = p/(rho*znd)
end subroutine q_conv

end module pri 
!---------------------------------------------------------
module  flux 
use data_kinds
use dims
implicit none

real    (kind=dp),allocatable :: delta_w(:,:,:)

type point
     real(kind=dp) ::x,y,z
end type point

type meshes
     real(kind=dp) ::x,y,z
     !type(point) :: p
end type meshes

type consrv
     real(kind=dp) :: rho,rhou,rhov,E
end type consrv

!type viscosity
!     real(kind=dp) :: mu,mut 
!end type viscosity

type finitevolume
     real(kind=dp) :: vol,dt
     real(kind=dp),dimension(ndim+2) ::  qvar,qold,dflux
     real(kind=dp) :: mu,mut
end type finitevolume


type  facs
     integer(kind=i2) :: is,js,ie,je,inc,jnc
end type facs


type fluxes
     real(kind=dp) :: nx,ny
     real(kind=dp) :: u_x,u_y,v_x,v_y,t_x,t_y
end type fluxes

type blocks
     integer (kind=i2):: imax,jmax,kmax
     integer (kind=i2):: cx,cy,cz
     integer (kind=i2):: is,ie,js,je,ks,ke
     integer (kind=i2):: bc(ndim*2),fc(ndim*2),mbc(ndim*2)
     type(facs),dimension(ndim*2,ndim*2)::fcs
     type(meshes),allocatable,dimension(:,:,:) :: mesh
     type(finitevolume),allocatable,dimension(:,:,:):: fv
     type(fluxes),allocatable,dimension(:,:,:):: zhi,eta
end type blocks

type walldata
     real(kind=dp) :: tauw(ndim),pwall(ndim),sf(ndim),cp(ndim) 
end type walldata 

type cells
     real(kind=dp) :: qvar,qold,dt,delq,delFlux,vol,mu
end type cells

type(blocks),allocatable,dimension(:):: blk
type(walldata),allocatable,dimension(:):: zhi_wall,eta_wall

!==============
end module flux
!==============

module st_end
use data_kinds
implicit none

integer (kind=i2):: is,ie,js,je,ks,ke
integer (kind=i2):: cx,cy,cz 
!integer (kind=i2):: ci1,cj1,ci2,cj2

contains
subroutine start_end(bk)
use flux
implicit none

integer (kind=i2):: bk

is=blk(bk)%is
ie=blk(bk)%ie
js=blk(bk)%js
je=blk(bk)%je
ks=blk(bk)%ks
ke=blk(bk)%ke

cx=blk(bk)%cx
cy=blk(bk)%cy
cz=blk(bk)%cz

end subroutine start_end


end module st_end
!---------------------------------------------------------

!=============
module Vect2D
!=============
implicit none

type points
   real(kind=8) :: x, y
end type points

type lines 
    type (points) :: ps, pe
end type lines 

type v2d      
   real(kind=8) :: x, y
end type v2d      

! -- INTERFACES -------------------------------------------------------------

interface v2dto1d
  module procedure v2d_2dto1d     
endinterface

interface vecxy
  module procedure v2d_xycomp     
endinterface

interface vmag
  module procedure v2d_magnitude
endinterface

interface operator(+)
  module procedure v2d_addition
endinterface

interface operator(-)
  module procedure v2d_substraction, v2d_opp
endinterface

interface operator(*)
  module procedure v2d_multiplydp
endinterface

interface operator(/)
  module procedure v2d_divisiondp
endinterface

interface operator(.scal.)
  module procedure v2d_scalar_product
endinterface

interface operator(.vect.)
  module procedure v2d_vectorial_product
endinterface

interface rot
  module procedure   v2d_rot_a
endinterface

interface norm
  module procedure v2d_nrm           
endinterface

interface intvect
  module procedure v2d_int           
endinterface



contains

!------------------------------------------------------------------------------!
! Initialise vector components 
!------------------------------------------------------------------------------!
type(v2d) function v2d_2dto1d(v)
implicit none
type(v2d), intent(in) :: v
!real(kind=8) :: x, y  

  v2d_2dto1d%x = v%x
  v2d_2dto1d%y = v%y

endfunction v2d_2dto1d

!------------------------------------------------------------------------------!
! Initialise vector components 
!------------------------------------------------------------------------------!
type(v2d) function v2d_xycomp(x, y)
implicit none
real(kind=8), intent(in) :: x, y  

  v2d_xycomp%x = x 
  v2d_xycomp%y = y

endfunction v2d_xycomp

!------------------------------------------------------------------------------!
! Add two vectors 
!------------------------------------------------------------------------------!
type(v2d) function v2d_addition(v1, v2)
implicit none
type(v2d), intent(in) :: v1, v2

  v2d_addition%x = v1%x + v2%x 
  v2d_addition%y = v1%y + v2%y 

endfunction v2d_addition

!------------------------------------------------------------------------------!
! Difference of two vectors
!------------------------------------------------------------------------------!
type(v2d) function v2d_substraction(v1, v2)
implicit none
type(v2d), intent(in) :: v1, v2

  v2d_substraction%x = v1%x - v2%x 
  v2d_substraction%y = v1%y - v2%y 

endfunction v2d_substraction

!------------------------------------------------------------------------------!
! Fonction : calcul de l'oppose d'un vecteur
!------------------------------------------------------------------------------!
type(v2d) function v2d_opp(v)
implicit none
type(v2d), intent(in) :: v

  v2d_opp%x = - v%x 
  v2d_opp%y = - v%y 

endfunction v2d_opp

!------------------------------------------------------------------------------!
! Multiply a vector by real no. 
!------------------------------------------------------------------------------!
type(v2d) function v2d_multiplydp(x, v)
implicit none
real(kind=8),   intent(in) :: x
type(v2d), intent(in) :: v

  v2d_multiplydp%x = x * v%x 
  v2d_multiplydp%y = x * v%y 

endfunction v2d_multiplydp

!------------------------------------------------------------------------------!
! Divide a vector by real 
!------------------------------------------------------------------------------!
type(v2d) function v2d_divisiondp(v,x)
implicit none
real(kind=8),   intent(in) :: x
type(v2d), intent(in) :: v

  v2d_divisiondp%x = v%x / x 
  v2d_divisiondp%y = v%y / x

endfunction v2d_divisiondp

!------------------------------------------------------------------------------!
! Magnitude of a vector        
!------------------------------------------------------------------------------!
real(kind=8) function v2d_magnitude(v)
implicit none
type(v2d), intent(in) :: v

  v2d_magnitude = dsqrt(v%x*v%x + v%y*v%y)

endfunction v2d_magnitude 

!------------------------------------------------------------------------------!
! Dot product of a vector                   
!------------------------------------------------------------------------------!
real(kind=8) function v2d_scalar_product(v1, v2)
implicit none
type(v2d), intent(in) :: v1, v2

  v2d_scalar_product = v1%x*v2%x + v1%y*v2%y

endfunction v2d_scalar_product

!------------------------------------------------------------------------------!
! Cross product of  a vector
!------------------------------------------------------------------------------!
real(kind=8) function v2d_vectorial_product(v1, v2)
implicit none
type(v2d), intent(in) :: v1, v2

  v2d_vectorial_product = v1%x*v2%y - v1%y*v2%x

endfunction v2d_vectorial_product

!------------------------------------------------------------------------------!
! Rotation of a vector               
!------------------------------------------------------------------------------!
type(v2d) function v2d_nrm(v)
implicit none
type(v2d), intent(in) :: v

  v2d_nrm%x =   v%y 
  v2d_nrm%y = - v%x 

endfunction v2d_nrm

!------------------------------------------------------------------------------!
! Fonction : rotation d'un vecteur
!------------------------------------------------------------------------------!
type(v2d) function v2d_rot_a(v, a)
implicit none
type(v2d), intent(in) :: v  ! vecteur a tourner
real(kind=8), intent(in) :: a  ! angle de rotation en radians
real(kind=8) :: ca, sa

  ca = dcos(a)
  sa = dsin(a)
  v2d_rot_a%x = v%x*ca - v%y*sa
  v2d_rot_a%y = v%x*sa + v%y*ca

endfunction v2d_rot_a

!------------------------------------------------------------------------------!
! Intersection of two vectors               
!------------------------------------------------------------------------------!
type(v2d) function v2d_int(la,lb)
implicit none
type(lines), intent(in) :: la,lb
type(v2d) :: a,b,c
real(kind=8)::ratio

  a=vecxy(la%pe%x-la%ps%x,la%pe%y-la%ps%y)
  b=vecxy(lb%pe%x-lb%ps%x,lb%pe%y-lb%ps%y)
  c=vecxy(lb%ps%x-la%ps%x,lb%ps%y-la%ps%y)

  !ratio=(c%x*a%y-c%y*a%x)/(a%y*b%x-a%x*b%y)
  ratio=(c%x*b%y-c%y*b%x)/(b%y*a%x-b%x*a%y)
  v2d_int%x = a%x*ratio
  v2d_int%y = a%y*ratio

endfunction v2d_int

!------------------------------------------------------------------------------!

end module  Vect2D

