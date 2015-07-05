! bc_vis.f.old
!============================================================
subroutine boundary_vis(bk)
!============================================================

!defines the fictitous cell values of conserved variable qvar
!for wall boundary, computes mu (viscosity) for the entire field
!and calculates the boundary fluxes fin_zhiv,fin_etav

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!nx    = same as n for this code
!ny    = same as n for this code
!mbc   = boundary condition vector
!n_zhi = normals on the xi faces
!n_eta = normals on the eta faces
!qvar  = conserved variable q
!fin_zhiv = fluxes on the zhi faces (output)
!fin_etav = fluxes on the eta faces (output)
!x = x coordinates
!y = y coordinates
!area = cell areas
!mu   = viscosity
!===============================================================
use st_end
use flux
implicit none

!integer(kind=i2):: nx,ny,nz
integer(kind=i2):: bk,i,j,k,icell,jcell
integer(kind=i2):: cxp1,cxm1,cyp1,cym1,idir
!real(kind=dp):: q1(ndim+2),q2(ndim+2),q3(ndim+2)
real(kind=dp):: q2(ndim+2),q3(ndim+2)
real(kind=dp):: aml(ndim),amr(ndim),amfr(ndim)
real(kind=dp):: qleft(ndim+2),qright(ndim+2)
real(kind=dp):: fvintr(ndim+2)
k=1
!do bk = 1,nblk
   call start_end(bk)
cxp1 = ie+two !cx+1
cxm1 = ie   !cx-1
cyp1 = je+two !cy+1
cym1 = je   !cy-1

!xi-direction
idir = one
if(blk(bk)%mbc(1) .eq. 2) then
  icell=is-one
  do j = js,je
    q2(:) = blk(bk)%fv(icell+1,j,k)%qvar(:)
    call bound_type_vis(q2,q3)
    blk(bk)%fv(icell,j,k)%qvar(:) = q3(:)
  enddo
endif

if(blk(bk)%mbc(2) .eq. 2) then
  icell=ie
  do j = js,je
    q2(:) = blk(bk)%fv(icell,j,k)%qvar(:)
    call bound_type_vis(q2,q3)
    blk(bk)%fv(icell+1,j,k)%qvar(:) = q3(:)
  enddo
endif

!eta-direction
idir = two
if(blk(bk)%mbc(3) .eq. 2) then
  jcell=js-one
  do i = is,ie
    q2(:) = blk(bk)%fv(i,jcell+1,k)%qvar(:)
    call bound_type_vis(q2,q3)
    blk(bk)%fv(i,jcell,k)%qvar(:) = q3(:)
  enddo
endif

if(blk(bk)%mbc(4) .eq. 2) then
  jcell=je  
  do i = is,ie
    q2(:) = blk(bk)%fv(i,jcell,k)%qvar(:)
    call bound_type_vis(q2,q3)
    blk(bk)%fv(i,jcell+1,k)%qvar(:) = q3(:)
  enddo
endif

!computes viscosity for the entire field
call transpt(bk,cx,cy)

!xi-direction boundary fluxes
!left boundary
idir = one
!if(blk(bk)%mbc(1) .eq. 0) then
!  fin_zhiv(is,:,:) = 0.d0
!else
  do j = js,je
    aml(1)  = blk(bk)%zhi(is-1,j,k)%nx
    aml(2)  = blk(bk)%zhi(is-1,j,k)%ny
    amr(1)  = blk(bk)%zhi(is,j,k)%nx
    amr(2)  = blk(bk)%zhi(is,j,k)%ny
    amfr(1)  = blk(bk)%zhi(is+1,j,k)%nx
    amfr(2)  = blk(bk)%zhi(is+1,j,k)%ny

    qleft(:)  = blk(bk)%fv(is-1,j,k)%qvar(:)
    qright(:) = blk(bk)%fv(is,j,k)%qvar(:)

    icell = is-one
    jcell = j
    
    call fluxvis_int_xi(bk,icell,jcell,1,qleft,qright,amr,aml,fvintr)

   !fin_zhiv(is,j,:) = -fvintr(:)
   !fin_zhi(icell+1,j,:) = fintr(:)
   blk(bk)%fv(is,j,k)%dflux(1:ndim+2)=blk(bk)%fv(is,j,k)%dflux(1:ndim+2)+fvintr(1:ndim+2)
  enddo
!endif

!right boundary
!if(blk(bk)%mbc(2) .eq. 0) then
!  fin_zhiv(ie+1,:,:) = 0.d0
!else
  do j = js,je
    aml(1)  = blk(bk)%zhi(ie,j,k)%nx
    aml(2)  = blk(bk)%zhi(ie,j,k)%ny
    amr(1)  = blk(bk)%zhi(ie+1,j,k)%nx
    amr(2)  = blk(bk)%zhi(ie+1,j,k)%ny
    amfr(1)  = blk(bk)%zhi(ie+2,j,k)%nx
    amfr(2)  = blk(bk)%zhi(ie+2,j,k)%ny

    qleft(:)  = blk(bk)%fv(ie,j,k)%qvar(:)
    qright(:) = blk(bk)%fv(ie+1,j,k)%qvar(:)

    icell = ie 
    jcell = j

    call fluxvis_int_xi(bk,icell,jcell,1,qleft,qright,amr,aml,fvintr)

    !fin_zhiv(ie+1,j,:) = -fvintr(:)
    !fin_zhi(icell+1,j,:) = fintr(:)
    blk(bk)%fv(ie,j,k)%dflux(1:ndim+2)=blk(bk)%fv(ie,j,k)%dflux(1:ndim+2)-fvintr(1:ndim+2)

  enddo
!endif

!eta-direction boundary fluxes
!left boundary
idir = two
!if(blk(bk)%mbc(3) .eq. 0) then
!  fin_etav(:,js,:) = 0.d0
!else
  do i = is,ie
    aml(1)  = blk(bk)%eta(i,js-1,k)%nx
    aml(2)  = blk(bk)%eta(i,js-1,k)%ny
    amr(1)  = blk(bk)%eta(i,js,k)%nx
    amr(2)  = blk(bk)%eta(i,js,k)%ny
    amfr(1)  = blk(bk)%eta(i,js+1,k)%nx
    amfr(2)  = blk(bk)%eta(i,js+1,k)%ny

    qleft(:)  = blk(bk)%fv(i,js-1,k)%qvar(:)
    qright(:) = blk(bk)%fv(i,js,k)%qvar(:)

    icell = i
    jcell = js-one

    call fluxvis_int_eta(bk,icell,jcell,2,qleft,qright,amr,aml,fvintr)

    !fin_etav(i,js,:) = -fvintr(:)
    blk(bk)%fv(i,js,k)%dflux(1:ndim+2)=blk(bk)%fv(i,js,k)%dflux(1:ndim+2)+fvintr(1:ndim+2)
  enddo
!endif

!right boundary
!if(blk(bk)%mbc(4) .eq. 0) then
!  fin_etav(:,je+1,:) = 0.d0
!else
  do i = is,ie
    aml(1)  = blk(bk)%eta(i,je,k)%nx
    aml(2)  = blk(bk)%eta(i,je,k)%ny
    amr(1)  = blk(bk)%eta(i,je+1,k)%nx
    amr(2)  = blk(bk)%eta(i,je+1,k)%ny
    amfr(1)  = blk(bk)%eta(i,je+2,k)%nx
    amfr(2)  = blk(bk)%eta(i,je+2,k)%ny

    qleft(:)  = blk(bk)%fv(i,je,k)%qvar(:)
    qright(:) = blk(bk)%fv(i,je+1,k)%qvar(:)

    icell = i
    jcell = je

    call fluxvis_int_eta(bk,icell,jcell,2,qleft,qright,amr,aml,fvintr)

    !fin_etav(i,je+1,:) = -fvintr(:)
    blk(bk)%fv(i,je,k)%dflux(1:ndim+2)=blk(bk)%fv(i,je,k)%dflux(1:ndim+2)-fvintr(1:ndim+2)
  enddo
!endif

!enddo
!===========================
end subroutine boundary_vis
!===========================


!=====================================
!subroutine bound_type_vis(ndim1,q2,q3)
subroutine bound_type_vis(q2,q3)
!=====================================

!defines the fictitous cell values for a given wall boundary cell

!===============================================================
!ndim  = space dimension (2)
!q2    = q in the cell at the boundary
!q3    = q in the fictitous cell
!===============================================================
use dims
use inf
use vis
use pri
implicit none

!integer(kind=i2):: ndim1
real(kind=dp):: q2(ndim+2),q3(ndim+2)
real(kind=dp):: gterm,r
!real(kind=dp):: rho2,u2,v2,e2,p2,t2,c2
!real(kind=dp):: rho3,u3,v3,e3,p3,t3,c3
real(kind=dp):: rho2,u2,v2,p2,t2
real(kind=dp):: rho3,u3,v3,p3,t3

gterm = gamma - 1.d0
r     = znd

call q_conv(q2)

p2   = p
rho2 = rho
t2   = t
u2   = u
v2   = v

u3  = -u2
v3  = -v2
p3  =  p2

if(iwall .eq. 0) then
  rho3  = rho2
 else
  !t3a  = 2.d0*t_w - t2
  t3  = 2.d0*t_w - t2
  if(t3 .lt. 0.d0) t3 = 0.01d0*t_w
  rho3  = p3/(r*t3)
endif

q3(1) = p3/gterm + 0.5d0*rho3*(u3*u3 +  v3*v3)
q3(2) = rho3
q3(3) = rho3*u3
q3(4) = rho3*v3

!=============================
end subroutine bound_type_vis
!=============================


!==========================================
subroutine transpt(bk,cx,cy)
!==========================================

!calculates viscosity in the entire flow-field using sutherlands
!law

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!cx    = same as n for this code
!cy    = same as n for this code
!qvar  = q in the domain
!mu    = viscosity
!===============================================================
use pri
use flux
implicit none

integer(kind=i2):: cx,cy
integer(kind=i2):: bk,i,j,k
real(kind=dp):: qleft(ndim+2),vis
k=1
do j = 1,cy
  do i = 1,cx
    qleft(:) = blk(bk)%fv(i,j,k)%qvar(:) 
    call q_conv(qleft)
    call viscous(t,vis)
    blk(bk)%fv(i,j,k)%mu = vis
  enddo
enddo

!======================
end subroutine transpt
!======================


!============================
subroutine viscous(ts,mu)
!============================

!sutherlands law

!===============================================================
!temp = temperature in the cell
!vis  = viscosity in the cell (output)
!===============================================================
use inf
!use vis
implicit none

real(kind=dp):: ts,mu
real(kind=dp):: sut,term1,term2

!sut   = 110.4d0/t_inf
!term1 = dsqrt(ts*ts*ts)
!term2 = (1.d0 + sut)/(ts + sut)
!mu  = term1*term2

sut=110.4d0

term1=(ts/t_inf)**1.5
term2=(t_inf+sut)/(ts+sut)
mu  = term1*term2

!======================
end subroutine viscous
!======================
