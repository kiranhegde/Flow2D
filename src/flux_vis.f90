! flux_vis.f.old
!=====================================================================================
subroutine fluxvis_int_xi(bk,icell,jcell,idir,qleft,qright,amr,aml,fvintr)  
!=====================================================================================
use flux
use pri
use vis
use inf
use st_end
!viscous flux calculation for a given xi cell interfaces

implicit none
integer(kind=i2)::bk
integer(kind=i2):: icell,jcell,idir,jterm,k
real(kind=dp):: qleft(ndim+2),qright(ndim+2),fvintr(ndim+2)
real(kind=dp):: amr(ndim),aml(ndim)
real(kind=dp):: uleft,vleft,tleft,visl,voll
real(kind=dp):: uright,vright,tright,visr,volr
real(kind=dp):: zhi_x,zhi_y
real(kind=dp):: com_fact,ratio
real(kind=dp):: tau_xx,tau_yy,tau_xy,tau_yx
real(kind=dp):: term1,pr_eff,term,term2,u_int,v_int
real(kind=dp):: beta_x,beta_y
real(kind=dp):: vistl,vistr 
real(kind=dp) :: u_x,u_y,v_x,v_y,t_x,t_y

k=1

call q_conv(qleft)
uleft = u
vleft = v
tleft = t
!visl  = mu(icell,jcell)
visl  = blk(bk)%fv(icell,jcell,k)%mu

call q_conv(qright)
uright = u
vright = v
tright = t
!visr   = mu(icell+1,jcell)
visr  = blk(bk)%fv(icell+1,jcell,k)%mu
voll = blk(bk)%fv(icell  ,jcell,k)%vol
volr = blk(bk)%fv(icell+1,jcell,k)%vol

zhi_x = amr(1)
zhi_y = amr(2)

call deriv_int_xi(bk,icell,jcell,idir,qleft,qright,amr,aml) 

u_x=blk(bk)%zhi(icell,jcell,k)%u_x
u_y=blk(bk)%zhi(icell,jcell,k)%u_y
v_x=blk(bk)%zhi(icell,jcell,k)%v_x
v_y=blk(bk)%zhi(icell,jcell,k)%v_y
t_x=blk(bk)%zhi(icell,jcell,k)%t_x
t_y=blk(bk)%zhi(icell,jcell,k)%t_y

com_fact = 2.d0/3.d0*(u_x + v_y)

tau_xx = 2.d0*u_x - com_fact
tau_yy = 2.d0*v_y - com_fact
tau_xy = u_y + v_x
tau_yx = tau_xy

ratio = 0.5d0*(visl/voll + visr/volr)/pr_l

!kkr
!vistl = visturb(icell  ,jcell)
!vistr = visturb(icell+1,jcell)
vistl  = blk(bk)%fv(icell,jcell,k)%mut
vistr  = blk(bk)%fv(icell+1,jcell,k)%mut

ratio = ratio + 0.5d0*(vistl/voll + vistr/volr)/pr_t
visl = visl + vistl
visr = visr + vistr
!kkr

!calculating the viscous interface fluxes
term1  = 0.5d0*(visl/voll + visr/volr)
pr_eff = term1/ratio
term   = term1/re_l
term2  = 1.d0/((gamma-1.d0)*m_inf*m_inf*pr_eff)
u_int  = 0.5d0*(uright + uleft)
v_int  = 0.5d0*(vright + vleft)

beta_x = u_int*tau_xx + v_int*tau_xy + term2*t_x
beta_y = u_int*tau_yx + v_int*tau_yy + term2*t_y
!
fvintr(2) = 0.d0

fvintr(3) = term*(zhi_x*tau_xx + zhi_y*tau_xy)
fvintr(4) = term*(zhi_x*tau_yx + zhi_y*tau_yy)
fvintr(1) = term*(zhi_x*beta_x + zhi_y*beta_y)

end subroutine fluxvis_int_xi



!===============================================================
subroutine fluxvis_int_eta(bk,icell,jcell,idir,qleft,qright,amr,aml,fvintr)
!===============================================================

!viscous flux calculation for a given eta cell interfaces
use flux
use pri
use vis
use inf
use st_end

implicit none
integer(kind=i2):: bk
integer(kind=i2):: icell,jcell,idir,jterm,k
real(kind=dp):: qleft(ndim+2),qright(ndim+2),fvintr(ndim+2)
real(kind=dp):: amr(ndim),aml(ndim)
real(kind=dp):: uleft,vleft,tleft,visl,voll
real(kind=dp):: uright,vright,tright,visr,volr
real(kind=dp):: eta_x,eta_y
real(kind=dp):: com_fact,ratio,pr_eff
real(kind=dp):: tau_xx,tau_yy,tau_xy,tau_yx
real(kind=dp):: term,term1,term2,u_int,v_int
real(kind=dp):: beta_x,beta_y!,fvintr(ndim+2)
real(kind=dp):: vistl,vistr 
real(kind=dp) :: u_x,u_y,v_x,v_y,t_x,t_y

k=1
call q_conv(qleft)
uleft = u
vleft = v
tleft = t
!visl  = mu(icell,jcell)
visl  = blk(bk)%fv(icell,jcell,k)%mu

call q_conv(qright)
uright = u
vright = v
tright = t
!visr   = mu(icell,jcell)
visr  = blk(bk)%fv(icell,jcell+1,k)%mu
voll = blk(bk)%fv(icell,jcell  ,k)%vol
volr = blk(bk)%fv(icell,jcell+1,k)%vol

eta_x = amr(1)
eta_y = amr(2)


call deriv_int_eta(bk,icell,jcell,idir,qleft,qright,amr,aml) 

u_x=blk(bk)%eta(icell,jcell,k)%u_x
u_y=blk(bk)%eta(icell,jcell,k)%u_y
v_x=blk(bk)%eta(icell,jcell,k)%v_x
v_y=blk(bk)%eta(icell,jcell,k)%v_y
t_x=blk(bk)%eta(icell,jcell,k)%t_x
t_y=blk(bk)%eta(icell,jcell,k)%t_y

com_fact = 2.d0/3.d0*(u_x + v_y)

tau_xx = 2.d0*u_x - com_fact
tau_yy = 2.d0*v_y - com_fact
tau_xy = u_y + v_x
tau_yx = tau_xy

ratio = 0.5d0*(visl/voll + visr/volr)/pr_l

! kkr
!vistl = visturb(icell,jcell  )
!vistr = visturb(icell,jcell+1)
vistl  = blk(bk)%fv(icell,jcell,k)%mut
vistr  = blk(bk)%fv(icell,jcell+1,k)%mut

ratio = ratio + 0.5d0*(vistl/voll + vistr/volr)/pr_t
visl = visl + vistl
visr = visr + vistr

!ratio = ratio + 0.5d0*(vistl/voll + vistr/volr)/pr_t
!visl = visl + vistl
!visr = visr + vistr
! kkr

term1  = 0.5d0*(visl/voll + visr/volr)
pr_eff = term1/ratio
term   = term1/re_l
term2  = 1.d0/((gamma-1.d0)*m_inf*m_inf*pr_eff)
u_int  = 0.5d0*(uright + uleft)
v_int  = 0.5d0*(vright + vleft)

beta_x = u_int*tau_xx + v_int*tau_xy + term2*t_x
beta_y = u_int*tau_yx + v_int*tau_yy + term2*t_y

fvintr(2) = 0.d0

fvintr(3) = term*(eta_x*tau_xx + eta_y*tau_xy)
fvintr(4) = term*(eta_x*tau_yx + eta_y*tau_yy)
fvintr(1) = term*(eta_x*beta_x + eta_y*beta_y)
!=============================
end subroutine fluxvis_int_eta
!=============================


!=====================================================================================
subroutine deriv_int_xi(bk,icell,jcell,idir,qleft,qright,amr,aml)  
!=====================================================================================
use flux
use pri
use vis
use inf
use st_end
!viscous flux calculation for a given xi cell interfaces
implicit none
integer(kind=i2):: bk,icell,jcell,idir,jterm,k
real(kind=dp):: qleft(ndim+2),qright(ndim+2)
real(kind=dp):: amr(ndim),aml(ndim)
real(kind=dp):: uleft,vleft,tleft,visl,voll
real(kind=dp):: uright,vright,tright,visr,volr
real(kind=dp):: zhi_x,zhi_y,u_zhi,v_zhi,t_zhi
real(kind=dp):: eta_x1,eta_x2,eta_x3,eta_x4
real(kind=dp):: eta_y1,eta_y2,eta_y3,eta_y4
real(kind=dp):: ubotl,vbotl,tbotl,ubotr,vbotr,tbotr
real(kind=dp):: utopl,vtopl,ttopl,utopr,vtopr,ttopr
real(kind=dp):: u_eta1,u_eta2,u_eta3,u_eta4
real(kind=dp):: v_eta1,v_eta2,v_eta3,v_eta4
real(kind=dp):: t_eta1,t_eta2,t_eta3,t_eta4
real(kind=dp):: eta_xu_eta,eta_yu_eta
real(kind=dp):: eta_xv_eta,eta_yv_eta
real(kind=dp):: eta_xt_eta, eta_yt_eta

k=1
call q_conv(qleft)
uleft = u
vleft = v
tleft = t

call q_conv(qright)
uright = u
vright = v
tright = t

zhi_x = amr(1)
zhi_y = amr(2)

eta_x1 = blk(bk)%eta(icell  ,jcell  ,k)%nx
eta_x2 = blk(bk)%eta(icell  ,jcell+1,k)%nx
eta_x3 = blk(bk)%eta(icell+1,jcell+1,k)%nx
eta_x4 = blk(bk)%eta(icell+1,jcell  ,k)%nx
eta_y1 = blk(bk)%eta(icell  ,jcell  ,k)%ny
eta_y2 = blk(bk)%eta(icell  ,jcell+1,k)%ny
eta_y3 = blk(bk)%eta(icell+1,jcell+1,k)%ny
eta_y4 = blk(bk)%eta(icell+1,jcell  ,k)%ny

u_zhi   = uright - uleft
v_zhi   = vright - vleft
t_zhi   = tright - tleft

do jterm = 1,ndim+2
  qleft(jterm)  = blk(bk)%fv(icell,jcell-1,k)%qvar(jterm)
  qright(jterm) = blk(bk)%fv(icell+1,jcell-1,k)%qvar(jterm)
enddo

call q_conv(qleft)
ubotl = u
vbotl = v
tbotl = t

call q_conv(qright)
ubotr = u
vbotr = v
tbotr = t

do jterm = 1,ndim+2
  qleft(jterm)  = blk(bk)%fv(icell,jcell+1,k)%qvar(jterm)
  qright(jterm) = blk(bk)%fv(icell+1,jcell+1,k)%qvar(jterm)
enddo

call q_conv(qleft)
utopl = u
vtopl = v
ttopl = t

call q_conv(qright)
utopr = u
vtopr = v
ttopr = t

u_eta1 = uleft  - ubotl
u_eta2 = utopl  - uleft
u_eta3 = utopr  - uright
u_eta4 = uright - ubotr
v_eta1 = vleft  - vbotl
v_eta2 = vtopl  - vleft
v_eta3 = vtopr  - vright
v_eta4 = vright - vbotr
t_eta1 = tleft  - tbotl
t_eta2 = ttopl  - tleft
t_eta3 = ttopr  - tright
t_eta4 = tright - tbotr

eta_xu_eta = 0.25d0*( eta_x1*u_eta1 + eta_x2*u_eta2 &
                    + eta_x3*u_eta3 + eta_x4*u_eta4)
eta_yu_eta = 0.25d0*( eta_y1*u_eta1 + eta_y2*u_eta2 &
                    + eta_y3*u_eta3 + eta_y4*u_eta4)
eta_xv_eta = 0.25d0*( eta_x1*v_eta1 + eta_x2*v_eta2 &
                    + eta_x3*v_eta3 + eta_x4*v_eta4)
eta_yv_eta = 0.25d0*( eta_y1*v_eta1 + eta_y2*v_eta2  &
                    + eta_y3*v_eta3 + eta_y4*v_eta4)
eta_xt_eta = 0.25d0*( eta_x1*t_eta1 + eta_x2*t_eta2 &
                    + eta_x3*t_eta3 + eta_x4*t_eta4)
eta_yt_eta = 0.25d0*( eta_y1*t_eta1 + eta_y2*t_eta2 &
                    + eta_y3*t_eta3 + eta_y4*t_eta4)

blk(bk)%zhi(icell,jcell,k)%u_x= zhi_x*u_zhi + eta_xu_eta
blk(bk)%zhi(icell,jcell,k)%u_y= zhi_y*u_zhi + eta_yu_eta
blk(bk)%zhi(icell,jcell,k)%v_x= zhi_x*v_zhi + eta_xv_eta
blk(bk)%zhi(icell,jcell,k)%v_y= zhi_y*v_zhi + eta_yv_eta
blk(bk)%zhi(icell,jcell,k)%t_x= zhi_x*t_zhi + eta_xt_eta
blk(bk)%zhi(icell,jcell,k)%t_y= zhi_y*t_zhi + eta_yt_eta

!==========================
end subroutine deriv_int_xi
!==========================


!===============================================================
subroutine deriv_int_eta(bk,icell,jcell,idir,qleft,qright,amr,aml)
!===============================================================

!viscous flux calculation for a given eta cell interfaces
use flux
use pri
use vis
use inf
use st_end
implicit none
integer(kind=i2):: bk,k,icell,jcell,idir,jterm
real(kind=dp):: qleft(ndim+2),qright(ndim+2)
real(kind=dp):: amr(ndim),aml(ndim)
real(kind=dp):: uleft,vleft,tleft
real(kind=dp):: uright,vright,tright
real(kind=dp):: eta_x,eta_y,u_eta,v_eta,t_eta
real(kind=dp):: zhi_x1,zhi_x2,zhi_x3,zhi_x4
real(kind=dp):: zhi_y1,zhi_y2,zhi_y3,zhi_y4
real(kind=dp):: ubotl,vbotl,tbotl,ubotr,vbotr,tbotr
real(kind=dp):: utopl,vtopl,ttopl,utopr,vtopr,ttopr
real(kind=dp):: u_zhi1,u_zhi2,u_zhi3,u_zhi4
real(kind=dp):: v_zhi1,v_zhi2,v_zhi3,v_zhi4
real(kind=dp):: t_zhi1,t_zhi2,t_zhi3,t_zhi4
real(kind=dp):: zhi_xu_zhi,zhi_yu_zhi
real(kind=dp):: zhi_xv_zhi,zhi_yv_zhi
real(kind=dp):: zhi_xt_zhi,zhi_yt_zhi
k=1
call q_conv(qleft)
uleft = u
vleft = v
tleft = t

call q_conv(qright)
uright = u
vright = v
tright = t

eta_x = amr(1)
eta_y = amr(2)

zhi_x1 = blk(bk)%zhi(icell  ,jcell  ,k)%nx
zhi_x2 = blk(bk)%zhi(icell  ,jcell+1,k)%nx
zhi_x3 = blk(bk)%zhi(icell+1,jcell+1,k)%nx
zhi_x4 = blk(bk)%zhi(icell+1,jcell  ,k)%nx
zhi_y1 = blk(bk)%zhi(icell  ,jcell  ,k)%ny
zhi_y2 = blk(bk)%zhi(icell  ,jcell+1,k)%ny
zhi_y3 = blk(bk)%zhi(icell+1,jcell+1,k)%ny
zhi_y4 = blk(bk)%zhi(icell+1,jcell  ,k)%ny

u_eta   = uright - uleft
v_eta   = vright - vleft
t_eta   = tright - tleft

do jterm = 1,ndim+2
  qleft(jterm)  = blk(bk)%fv(icell-1,jcell,k)%qvar(jterm)
  qright(jterm) = blk(bk)%fv(icell+1,jcell,k)%qvar(jterm)
enddo

call q_conv(qleft)
ubotl = u
vbotl = v
tbotl = t

call q_conv(qright)
ubotr = u
vbotr = v
tbotr = t

do jterm = 1,ndim+2
  qleft(jterm)  = blk(bk)%fv(icell-1,jcell+1,k)%qvar(jterm)
  qright(jterm) = blk(bk)%fv(icell+1,jcell+1,k)%qvar(jterm)

enddo

call q_conv(qleft)
utopl = u
vtopl = v
ttopl = t

call q_conv(qright)
utopr = u
vtopr = v
ttopr = t

u_zhi1 = uleft  - ubotl
u_zhi2 = uright - utopl
u_zhi3 = utopr  - uright
u_zhi4 = ubotr  - uleft
v_zhi1 = vleft  - vbotl
v_zhi2 = vright - vtopl
v_zhi3 = vtopr  - vright
v_zhi4 = vbotr  - vleft
t_zhi1 = tleft  - tbotl
t_zhi2 = tright - ttopl
t_zhi3 = ttopr  - tright
t_zhi4 = tbotr  - tleft

zhi_xu_zhi = 0.25d0*( zhi_x1*u_zhi1 + zhi_x2*u_zhi2 &
                    + zhi_x3*u_zhi3 + zhi_x4*u_zhi4)
zhi_yu_zhi = 0.25d0*( zhi_y1*u_zhi1 + zhi_y2*u_zhi2 &
                    + zhi_y3*u_zhi3 + zhi_y4*u_zhi4)
zhi_xv_zhi = 0.25d0*( zhi_x1*v_zhi1 + zhi_x2*v_zhi2 &
                    + zhi_x3*v_zhi3 + zhi_x4*v_zhi4) 
zhi_yv_zhi = 0.25d0*( zhi_y1*v_zhi1 + zhi_y2*v_zhi2 &
                    + zhi_y3*v_zhi3 + zhi_y4*v_zhi4)  
zhi_xt_zhi = 0.25d0*( zhi_x1*t_zhi1 + zhi_x2*t_zhi2&
                    + zhi_x3*t_zhi3 + zhi_x4*t_zhi4)
zhi_yt_zhi = 0.25d0*( zhi_y1*t_zhi1 + zhi_y2*t_zhi2 &
                    + zhi_y3*t_zhi3 + zhi_y4*t_zhi4)  

blk(bk)%eta(icell,jcell,k)%u_x= zhi_xu_zhi + eta_x*u_eta
blk(bk)%eta(icell,jcell,k)%u_y= zhi_yu_zhi + eta_y*u_eta
blk(bk)%eta(icell,jcell,k)%v_x= zhi_xv_zhi + eta_x*v_eta
blk(bk)%eta(icell,jcell,k)%v_y= zhi_yv_zhi + eta_y*v_eta
blk(bk)%eta(icell,jcell,k)%t_x= zhi_xt_zhi + eta_x*t_eta
blk(bk)%eta(icell,jcell,k)%t_y= zhi_yt_zhi + eta_y*t_eta

!==========================
end subroutine deriv_int_eta
!==========================
