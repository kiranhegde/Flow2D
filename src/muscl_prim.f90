! muscl_prim.f.old
!==========================================================
subroutine muscl(qleft,qfleft,qright,qfright,aml,amr,qpos,qneg,r1,r2,r3,r4)
!==========================================================

!muscl extrapolation using primitive variables

!===============================================================
!ndim  = space dimension (2)
!qleft   = q in the 1st left cell
!qfleft  = q in the 2nd left cell
!qright  = q in the 1st right cell
!qfright = q in the 2nd right cell
!aml,amr = normals, amr is the normal for face under
!           consideration
!qpos    = q+, extrapolated using muscl
!qneg    = q-, extrapolated using muscl
!r1,r2,r3,r4 = cell centre distances from the interface
!r1    = distance of the interface to the cell centre of 2nd cell
!        to the left
!r2    = distance of the interface to the cell centre of 1st cell
!        to the left
!r3    = distance of the interface to the cell centre of 1st cell
!        to the right
!r4    = distance of the interface to the cell centre of 2nd cell
!        to the right

!===============================================================
use dims
use inf
use pri

implicit none
integer(kind=i2):: ndim1,itvd,j

real(kind=dp):: qleft(ndim+2),qfleft(ndim+2)
real(kind=dp):: qright(ndim+2),qfright(ndim+2)
real(kind=dp):: amr(ndim),aml(ndim)
real(kind=dp):: r1,r2,r3,r4,iphi,h1,h2,h3,h4
real(kind=dp):: qpos(ndim+2),qneg(ndim+2)

integer(kind=i2):: jterm
real(kind=dp):: gterm,ametav(ndim),term1,term2,small,beta,s1,s2
real(kind=dp):: u_left,v_left,p_left,rho_left
real(kind=dp):: u_fleft,v_fleft,p_fleft,rho_fleft
real(kind=dp):: u_right,v_right,p_right,rho_right
real(kind=dp):: u_fright,v_fright,p_fright,rho_fright
real(kind=dp):: a1(ndim+2),a2(ndim+2),a3(ndim+2),a4(ndim+2)
real(kind=dp):: a11(ndim+2),a12(ndim+2),a21(ndim+2),a22(ndim+2)
real(kind=dp):: pneg,roneg,uneg,vneg,ppos,ropos,upos,vpos

itvd = 1
gterm = gamma - 1.0

!itvd = 0 from first-order
if(itvd .eq. 0) then
  qpos(:) = qright(:)
  qneg(:) = qleft(:)
  return
endif

!iphi = -1.0
!iphi = 0.0 
!iphi = 1.0
iphi = 1.0/3.0
small = 1e-12 

call q_conv(qleft)
rho_left = rho
u_left   = u
v_left   = v
p_left   = p

call q_conv(qright)
rho_right = rho
u_right   = u
v_right   = v
p_right   = p

call q_conv(qfleft)
rho_fleft = rho
u_fleft   = u
v_fleft   = v
p_fleft   = p

call q_conv(qfright)
rho_fright = rho
u_fright   = u
v_fright   = v
p_fright   = p

a1(1) = p_left   - p_fleft
a1(2) = rho_left - rho_fleft
a1(3) = u_left   - u_fleft
a1(4) = v_left   - v_fleft

a2(1) = p_right   - p_left
a2(2) = rho_right - rho_left
a2(3) = u_right   - u_left
a2(4) = v_right   - v_left

a3(1) = p_right   - p_left
a3(2) = rho_right - rho_left
a3(3) = u_right   - u_left
a3(4) = v_right   - v_left

a4(1) = p_fright   - p_right
a4(2) = rho_fright - rho_right
a4(3) = u_fright   - u_right
a4(4) = v_fright   - v_right

! For non uniform mesh
! Cell length method
! h - cell length
h2=2*r2
h3=2*r3
h1=2*(r1-h2)
h4=2*(r4-h3)
a1(:)=2.0*a1(:)/(h1+h2)
a2(:)=2.0*a2(:)/(h2+h3)
a3(:)=2.0*a3(:)/(h2+h3)
a4(:)=2.0*a4(:)/(h3+h4)


! van Albada limiter
do j = 1,ndim+2
   s1=(2*a2(j)*a1(j)+small)/(a2(j)*a2(j)+a1(j)*a1(j)+small)
   s2=(2*a4(j)*a3(j)+small)/(a4(j)*a4(j)+a3(j)*a3(j)+small)
   a11(j) =s1*(1.d0-s1*iphi)*a1(j)*h2
   a12(j) =s1*(1.d0+s1*iphi)*a2(j)*h2
   a21(j) =s2*(1.d0-s2*iphi)*a4(j)*h3
   a22(j) =s2*(1.d0+s2*iphi)*a3(j)*h3
enddo

pneg  = p_left   + 0.25d0*(a11(1) + a12(1))
roneg = rho_left + 0.25d0*(a11(2) + a12(2))
uneg  = u_left   + 0.25d0*(a11(3) + a12(3))
vneg  = v_left   + 0.25d0*(a11(4) + a12(4))

qneg(1) = pneg/gterm + 0.5d0*roneg*(uneg*uneg + vneg*vneg)
qneg(2) = roneg
qneg(3) = roneg*uneg
qneg(4) = roneg*vneg

ppos  = p_right   - 0.25d0*(a21(1) + a22(1))
ropos = rho_right - 0.25d0*(a21(2) + a22(2))
upos  = u_right   - 0.25d0*(a21(3) + a22(3))
vpos  = v_right   - 0.25d0*(a21(4) + a22(4))

qpos(1) = ppos/gterm + 0.5d0*ropos*(upos*upos + vpos*vpos)
qpos(2) = ropos
qpos(3) = ropos*upos
qpos(4) = ropos*vpos


if ((pneg .le. 0.d0) .or. (ppos .le. 0.d0)) then
  do jterm = 1,ndim+2
    qneg(jterm) = qleft(jterm)
    qpos(jterm) = qright(jterm)
  enddo
endif

!====================
end subroutine muscl
!====================
