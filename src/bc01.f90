! bc.f.old
!========================================================
subroutine boundary(bk)
!========================================================
!calculates the fictitous cell values of conserved variable
!qvar as well as the invisid boundary interface fluxes fin_zhi and
!fin_eta for the complete boundary
!
!delta_w is dummy in the present implementation
!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!cx    = same as n for this code
!cy    = same as n for this code
!mbc   = boundary condition vector
!n_zhi = normals on the xi faces
!n_eta = normals on the eta faces
!qvar  = conserved variable q
!fin_zhi = fluxes on the zhi faces (output)
!fin_eta = fluxes on the eta faces (output)
!x = x coordinates
!y = y coordinates
!area = cell areas
!fintr = flux at the given cell interface
!       delta_w = dummy, not used
!===============================================================
use flux
use st_end
implicit none

integer(kind=i2):: bk,ibound,i
!integer(kind=i2):: idir,iface,iside,ibound
integer(kind=i2):: k,icell,jcell
real(kind=dp):: q1(ndim+2),q2(ndim+2),q3(ndim+2),fintr(ndim+2)
real(kind=dp):: aml(ndim),amr(ndim),amfr(ndim)
k=1

!do bk = one,nblk
call start_end(bk)

do i=1,4
   ibound = blk(bk)%mbc(i)
   if (ibound /= 2) call BCond(i,ibound)
enddo       
   
! Applying  wall BC at the end
do i=1,4
   ibound = blk(bk)%mbc(i)
   if (ibound == 2) call BCond(i,ibound)
enddo       
    
contains


subroutine BCond(i,ibound)
implicit none
integer(kind=i2):: i,idir,iface,iside,ibound


if(i==1) then
!xi direction; left boundary
idir   = one
iface  = one
iside  = one
!ibound = blk(bk)%mbc(1)
icell  = is-one
do jcell = js,je

   aml(1)  = blk(bk)%zhi(icell  ,jcell,k)%nx
   aml(2)  = blk(bk)%zhi(icell  ,jcell,k)%ny
   amr(1)  = blk(bk)%zhi(icell+1,jcell,k)%nx
   amr(2)  = blk(bk)%zhi(icell+1,jcell,k)%ny
   amfr(1) = blk(bk)%zhi(icell+2,jcell,k)%nx
   amfr(2) = blk(bk)%zhi(icell+2,jcell,k)%ny
     q1(:) = blk(bk)%fv( icell+2,jcell,k)%qvar(:)
     q2(:) = blk(bk)%fv( icell+1,jcell,k)%qvar(:)
     q3(:) = blk(bk)%fv( icell  ,jcell,k)%qvar(:)
 
   call bound(bk,idir,iside,icell,jcell,ibound,aml,amr,amfr,q1,q2,q3,fintr) 
                  
   blk(bk)%fv(icell  ,jcell,k)%qvar(:)         = q3(:)
   blk(bk)%fv(icell+1,jcell,k)%dflux(1:ndim+2) = blk(bk)%fv(icell+1,jcell,k)%dflux(1:ndim+2)-fintr(1:ndim+2)

enddo

elseif (i==2) then
!xi direction; right boundary
idir   = one
iface  = two
iside  = two
!ibound = blk(bk)%mbc(2)
icell  = ie 
do jcell = js,je

   aml(1)  = blk(bk)%zhi(icell  ,jcell,k)%nx
   aml(2)  = blk(bk)%zhi(icell  ,jcell,k)%ny
   amr(1)  = blk(bk)%zhi(icell+1,jcell,k)%nx
   amr(2)  = blk(bk)%zhi(icell+1,jcell,k)%ny
   amfr(1) = blk(bk)%zhi(icell+2,jcell,k)%nx
   amfr(2) = blk(bk)%zhi(icell+2,jcell,k)%ny
     q1(:) = blk(bk)%fv( icell-1,jcell,k)%qvar(:)
     q2(:) = blk(bk)%fv( icell  ,jcell,k)%qvar(:)
     q3(:) = blk(bk)%fv( icell+1,jcell,k)%qvar(:)
 
   call bound(bk,idir,iside,icell,jcell,ibound,aml,amr,amfr,q1,q2,q3,fintr)
 
   blk(bk)%fv(icell+1,jcell,k)%qvar(:) = q3(:)
   blk(bk)%fv(icell  ,jcell,k)%dflux(1:ndim+2)=blk(bk)%fv(icell,jcell,k)%dflux(1:ndim+2)+fintr(1:ndim+2)

enddo

elseif (i==3) then
!eta direction; left boundary
idir   = two
iface  = three
iside  = one
!ibound = blk(bk)%mbc(3)
jcell  = js-one
do icell = is,ie

   aml(1)  = blk(bk)%eta(icell,jcell  ,k)%nx
   aml(2)  = blk(bk)%eta(icell,jcell  ,k)%ny
   amr(1)  = blk(bk)%eta(icell,jcell+1,k)%nx
   amr(2)  = blk(bk)%eta(icell,jcell+1,k)%ny
   amfr(1) = blk(bk)%eta(icell,jcell+2,k)%nx
   amfr(2) = blk(bk)%eta(icell,jcell+2,k)%ny
     q1(:) = blk(bk)%fv( icell,jcell+2,k)%qvar(:)
     q2(:) = blk(bk)%fv( icell,jcell+1,k)%qvar(:)
     q3(:) = blk(bk)%fv( icell,jcell  ,k)%qvar(:)
 
   call bound(bk,idir,iside,icell,jcell,ibound,aml,amr,amfr,q1,q2,q3,fintr)
 
   blk(bk)%fv(icell,jcell  ,k)%qvar(:) = q3(:)
   blk(bk)%fv(icell,jcell+1,k)%dflux(1:ndim+2)=blk(bk)%fv(icell,jcell+1,k)%dflux(1:ndim+2)-fintr(1:ndim+2)

enddo

elseif (i==4) then
!eta direction; right boundary
idir   = two
iface  = four
iside  = two
!ibound = blk(bk)%mbc(4)
jcell  = je  
do icell = is,ie

   aml(1)  = blk(bk)%eta(icell,jcell  ,k)%nx
   aml(2)  = blk(bk)%eta(icell,jcell  ,k)%ny
   amr(1)  = blk(bk)%eta(icell,jcell+1,k)%nx
   amr(2)  = blk(bk)%eta(icell,jcell+1,k)%ny
   amfr(1) = blk(bk)%eta(icell,jcell+2,k)%nx
   amfr(2) = blk(bk)%eta(icell,jcell+2,k)%ny
     q1(:) = blk(bk)%fv( icell,jcell-1,k)%qvar(:)
     q2(:) = blk(bk)%fv( icell,jcell  ,k)%qvar(:)
     q3(:) = blk(bk)%fv( icell,jcell+1,k)%qvar(:)
 
   call bound(bk,idir,iside,icell,jcell,ibound,aml,amr,amfr,q1,q2,q3,fintr)
 
   blk(bk)%fv(icell,jcell+1,k)%qvar(:) = q3(:)
   blk(bk)%fv(icell,jcell  ,k)%dflux(1:ndim+2)=blk(bk)%fv(icell,jcell,k)%dflux(1:ndim+2)+fintr(1:ndim+2)

enddo
endif
end subroutine bcond

!enddo
!=======================
end subroutine boundary
!=======================




!==============================================================
subroutine bound(bk,idir,iside,icell,jcell,ibound,aml,amr,amfr,q1,q2,q3,fintr)
!==============================================================
!calculates the fictitous cell values of conserved variable
!qvar as well as the boundary interface fluxes fin_zhi and
!fin_eta for a given cell
!boundary types:
!ibound = 0 : block interface
!ibound = 1 : characteristic
!ibound = 2 : wall
!ibound = 3 : extrapolation
!ibound = 4 : freestream
!ibound = 5 : symmetry
!ibound = 6 : periodic
!ibound = 7 : outflow pressure specified
!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!cx    = same as n for this code
!cy    = same as n for this code
!idir  = 1 for xi direction, 2 for eta direction
!iside = 1 for left interface, 2 for right interface
!icell = xi index
!jcell = eta index
!ibound = boundary conditin (taken from vector mbc)
!aml,amr,amfr = normal on the left,right(present),far-right
!                interface
!q1 = q for second computation cell on the boundary
!q2 = q for first computation cell on the boundary
!q3 = fictitous cell q
!fintr = flux at the interface under consideration
!n_zhi = normals on the xi faces
!n_eta = normals on the eta faces
!qvar  = conserved variable q
!x = x coordinates
!y = y coordinates
!area = cell areas
!qvar = q for the complete domain
!       delta_w = dummy, not used
!===============================================================
use flux
use inf
use st_end
use pri
implicit none

integer(kind=i2):: bk,icell,jcell
integer(kind=i2):: idir,iside,ibound
real(kind=dp):: aml(ndim),amr(ndim),amfr(ndim)
real(kind=dp):: q1(ndim+2),q2(ndim+2),q3(ndim+2),q4(ndim+2)
real(kind=dp):: qpos(ndim+2),qneg(ndim+2),fintr(ndim+2)

integer(kind=i2):: i,j,k
real(kind=dp):: gamma_m1,grad_met,amx,amy
real(kind=dp):: rho1,u1,v1,e1,p1,t1,c1
real(kind=dp):: rho2,u2,v2,e2,p2,t2,c2
real(kind=dp):: rho3,u3,v3,e3,p3,t3
real(kind=dp):: pneg,roneg,uneg,vneg,cneg,contra_neg,velt_neg,hneg
real(kind=dp):: ppos,ropos,upos,vpos,cpos,contra_pos,velt_pos,hpos
real(kind=dp):: alambda1,alambda3,riem1,riem2,riem3,contra_rie
real(kind=dp):: ut_rie,vt_rie,velt_rie,ro_rie,p_rie,vel_ries,e_rie
real(kind=dp):: u_rie,v_rie,gterm2,vel_nor,pwall,fac_two,r,q_rie,h_rie
real(kind=dp):: term1,term2,term4,term5,term6
real(kind=dp):: qfleft(ndim+2),qleft(ndim+2)
real(kind=dp):: qright(ndim+2),qfright(ndim+2)
real(kind=dp):: r1,r2,r3,r4
k=1

call start_end(bk)

gamma_m1 = gamma-one
r        = znd

grad_met = dsqrt(amr(1)*amr(1) + amr(2)*amr(2))
amx      = amr(1)/grad_met
amy      = amr(2)/grad_met

call q_conv(q1)
rho1 = rho
u1   = u
v1   = v
e1   = e
p1   = p
t1   = t
c1   = a

call q_conv(q2)
rho2 = rho
u2   = u
v2   = v
e2   = e
p2   = p
t2   = t
c2   = a

call cell_dist(bk,cx,cy,icell,jcell,idir,iside,r1,r2,r3,r4)

select case(ibound)

!ibound = 0 : block interface
case(0)

  select case(idir)

!  (idir .eq. 1)
    case(1)

    select case(iside)

!      (iside .eq. 1)
      case(1)

        q3(:) = blk(bk)%fv(icell,jcell,k)%qvar(:)
        q4(:) = blk(bk)%fv(icell-1,jcell,k)%qvar(:)

        qfleft(:)  = q4(:)
        qleft(:)   = q3(:)
        qright(:)  = q2(:)
        qfright(:) = q1(:)

        call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

!       (iside .eq. 2)
      case(2)

        q3(:) = blk(bk)%fv(icell+1,jcell,k)%qvar(:)
        q4(:) = blk(bk)%fv(icell+2,jcell,k)%qvar(:)

        qfleft(:)  = q1(:)
        qleft(:)   = q2(:)
        qright(:)  = q3(:)
        qfright(:) = q4(:)

        call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

    end select     ! iside select

!  (idir .eq. 2)
    case(2)

    select case(iside)

!      (iside .eq. 1)
      case(1)

        q3(:) = blk(bk)%fv(icell,js-1,k)%qvar(:)
        q4(:) = blk(bk)%fv(icell,js-2,k)%qvar(:)

        qfleft(:)  = q4(:)
        qleft(:)   = q3(:)
        qright(:)  = q2(:)
        qfright(:) = q1(:)

        call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

!      (iside .eq. 2)
      case(2)

        q3(:) = blk(bk)%fv(icell,je+1,k)%qvar(:)
        q4(:) = blk(bk)%fv(icell,je+2,k)%qvar(:)

        qfleft(:)  = q1(:)
        qleft(:)   = q2(:)
        qright(:)  = q3(:)
        qfright(:) = q4(:)

        call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

    end select         !iside select

  end select

!ibound = 1 : characteristic
case(1)
!
! calculation of interior flow field values
!
  select case(iside)

!  (iside .eq. 1)
  case(1)
    qpos(:)  = q2(:)
    qneg(:)  = q3(:)

!  (iside .eq. 2)
  case(2)

    qneg(:)  = q2(:)
    qpos(:)  = q3(:)
  end select
!
! the positive state is obtained from the right state
!
!
! defining the left riemann state primitives
!
  call q_conv(qneg)
  pneg  = p
  roneg = rho
  uneg  = u
  vneg  = v
  cneg  = a
  hneg  = h
  contra_neg = uneg*amx + vneg*amy
  velt_neg  = dsqrt(u*u + v*v - contra_neg*contra_neg)
!
! defining the right riemann state primitives
!
  call q_conv(qpos)
  ppos  = p
  ropos = rho
  upos  = u
  vpos  = v
  cpos  = a
  hpos  = h
  contra_pos= upos*amx + vpos*amy
  velt_pos  = dsqrt(u*u + v*v - contra_pos*contra_pos)
!
!calculating the interface eigenvalues 1d
!
  select case(iside)

!  (iside .eq. 1)
  case(1)
    alambda1   = contra_neg - cneg
    alambda3   = contra_neg + cneg
!  (iside .eq. 2)
  case(2)
    alambda1   = contra_pos - cpos
    alambda3   = contra_pos + cpos
  end select
!
! computing the interface riemann invariants
!
  if(alambda1 .le. 0.d0) then
    riem3 = contra_pos - 2.d0*cpos/gamma_m1
   else
    riem3 = contra_neg - 2.d0*cneg/gamma_m1
  endif

  if(alambda3 .le. 0.d0) then
    riem1 = contra_pos + 2.d0*cpos/gamma_m1
   else
    riem1 = contra_neg + 2.d0*cneg/gamma_m1
  endif
!
! decoupling the riemann invariants at the interface
!
  contra_rie = 0.5d0*(riem1 + riem3)

  if(contra_rie .le. 0.d0) then
    riem2    = ppos/(ropos**gamma)
    ut_rie   = upos - contra_pos*amx
    vt_rie   = vpos - contra_pos*amy
    velt_rie = velt_pos
    h_rie    = hpos
    q_rie    = contra_rie*contra_rie + velt_rie*velt_rie
   else
    riem2    = pneg/(roneg**gamma)
    ut_rie   = uneg - contra_neg*amx
    vt_rie   = vneg - contra_neg*amy
    velt_rie = velt_neg
    h_rie    = hneg
    q_rie    = contra_rie*contra_rie + velt_rie*velt_rie
  endif

  term1      = gamma_m1*(riem1 - riem3)*0.25d0
!  term1      = 0.5d0*gamma_m1*(h_rie - q_rie)
  term2      = term1*term1/gamma/riem2
  ro_rie     = term2**(1.d0/gamma_m1)
  p_rie      = riem2*(ro_rie**gamma)
  vel_ries   = contra_rie**2 + velt_rie**2
  e_rie      = p_rie/gamma_m1 + 0.5d0*ro_rie*vel_ries
  u_rie      = ut_rie + contra_rie*amx
  v_rie      = vt_rie + contra_rie*amy
!
! the interface inviscid flux
!
  contra_rie = contra_rie*grad_met

  fintr(1)   = (e_rie + p_rie)*contra_rie
  fintr(2)   = ro_rie*contra_rie
  fintr(3)   = fintr(2)*u_rie + p_rie*amr(1)
  fintr(4)   = fintr(2)*v_rie + p_rie*amr(2)

! Was commented
!--------------------------------------
!  q3(1)  = 2.d0*e_rie  - q2(1)
!  q3(2)  = 2.d0*ro_rie - q2(2)
!  q3(3)  = 2.d0*ro_rie*u_rie - q2(3)
!  q3(4)  = 2.d0*ro_rie*v_rie - q2(4)
  !q3(5)  = 2.d0*ro_rie*w_rie - q2(5)
!--------------------------------------



!ibound = 2 : wall
case(2)

  gterm2 = gamma/(gamma - 1.d0)
!
! calculation of the interior riemann state
!
!
! calculating the normal velocity
!
  vel_nor = u2*amx + v2*amy
  vel_nor = -vel_nor*(-one)**iside
!
! calculates wall pressure by solving 1d riemann problem
!
  term4 = gamma_m1/2.d0*(vel_nor/c2)
  term5 = 1.d0 - term4
  term6 = term5**(2.d0*sngl(gterm2))
  pwall = p2*term6
  if(pwall .lt. 0.d0) pwall = p2
!  pwall = p2

!
! imposes proper b.c. for fictitious cell vector for
! consistent wall inviscid flux since it is required
! for the next interface calculation
!
  fac_two = 2.d0*vel_nor
  u3      = u2 - fac_two*amx
  v3      = v2 - fac_two*amy
  e3      = e2
  rho3    = rho2
!
!  calculating the fictitious cell vector
!

  q3(1) = e3
  q3(2) = rho3
  q3(3) = rho3*u3
  q3(4) = rho3*v3
!
! the inviscid wall flux vector
!
  fintr(1)   = 0.d0
  fintr(2)   = 0.d0
  fintr(3)   = pwall*amx*grad_met
  fintr(4)   = pwall*amy*grad_met
!


!ibound = 3 : extrapolation
case(3)

  q3(:) = q2(:)
  q4(:) = q3(:)

  select case(iside)

!  (iside .eq. 1)
  case(1)

    qfleft(:)  = q4(:)
    qleft(:)   = q3(:)
    qright(:)  = q2(:)
    qfright(:) = q1(:)

    call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

!  (iside .eq. 2)
  case(2)

    qfleft(:)  = q1(:)
    qleft(:)   = q2(:)
    qright(:)  = q3(:)
    qfright(:) = q4(:)

    call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

  end select

!ibound = 4 : freestream
case(4)

  q3(:) = qinf(:)
  q4(:) = q3(:)

  select case(iside)

!  (iside .eq. 1)
  case(1)

    qfleft(:)  = q4(:)
    qleft(:)   = q3(:)
    qright(:)  = q2(:)
    qfright(:) = q1(:)

    call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

!  (iside .eq. 2)
  case(2)

    qfleft(:)  = q1(:)
    qleft(:)   = q2(:)
    qright(:)  = q3(:)
    qfright(:) = q4(:)

    call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

  end select
!

!ibound = 5 : symmetry
case(5)

 gterm2 = gamma/(gamma - 1.d0)
 vel_nor = u2*amx + v2*amy
 fac_two = 2.d0*vel_nor
 e3      = e2
 u3      = u2 - fac_two*amx
 v3      = v2 - fac_two*amy
 p3      = p2
 rho3    = rho2
!
!  calculating the fictitious cell vector
!

  q3(1) = e3
  q3(2) = rho3
  q3(3) = rho3*u3
  q3(4) = rho3*v3
!
! the inviscid wall flux vector
!
  fintr(1)   = 0.d0
  fintr(2)   = 0.d0
  fintr(3)   = p2*amx*grad_met
  fintr(4)   = p2*amy*grad_met
!


!ibound = 6 : periodic
case(6)


  select case(idir)

!  (idir .eq. 1)
    case(1)

    select case(iside)

!      (iside .eq. 1)
      case(1)

        q3(:) = blk(bk)%fv(ie  ,jcell,k)%qvar(:)
        q4(:) = blk(bk)%fv(ie-1,jcell,k)%qvar(:)

        qfleft(:)  = q4(:)
        qleft(:)   = q3(:)
        qright(:)  = q2(:)
        qfright(:) = q1(:)

        call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

!       (iside .eq. 2)
      case(2)

        q3(:) = blk(bk)%fv(is,jcell,k)%qvar(:)
        q4(:) = blk(bk)%fv(is+1,jcell,k)%qvar(:)

        qfleft(:)  = q1(:)
        qleft(:)   = q2(:)
        qright(:)  = q3(:)
        qfright(:) = q4(:)

        call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

    end select     ! iside select

!  (idir .eq. 2)
    case(2)

    select case(iside)

!      (iside .eq. 1)
      case(1)

        q3(:) = blk(bk)%fv(icell,je  ,k)%qvar(:)
        q4(:) = blk(bk)%fv(icell,je-1,k)%qvar(:)

        qfleft(:)  = q4(:)
        qleft(:)   = q3(:)
        qright(:)  = q2(:)
        qfright(:) = q1(:)

        call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

!      (iside .eq. 2)
      case(2)

        q3(:) = blk(bk)%fv(icell,js,k)%qvar(:)
        q4(:) = blk(bk)%fv(icell,js+1,k)%qvar(:)

        qfleft(:)  = q1(:)
        qleft(:)   = q2(:)
        qright(:)  = q3(:)
        qfright(:) = q4(:)

        call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

    end select         !iside select



  end select


!ibound = 7 : outflow pressure specified
case(7)

  rho3 = rho2 + (pout-p2)/(c2*c2)
  u3   = u2 + amx*(p2-pout)/(rho2*c2)
  v3   = v2 + amy*(p2-pout)/(rho2*c2)
  e3   = pout/(gamma-one) + 0.5d0*rho3*(u3*u3 + v3*v3)

  q3(1) = e3
  q3(2) = rho3
  q3(3) = rho3*u3
  q3(4) = rho3*v3

  qfleft(:)  = q1(:)
  qleft(:)   = q2(:)
  qright(:)  = q3(:)
  qfright(:) = q3(:)

  call fluxinv_int(cx,cy,idir,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)

end select

!====================
end subroutine bound
!====================




!=============================================================
subroutine bound_delta(bk,nx,ny)
!=============================================================

!not used
use flux
use st_end
use inf
use pri
implicit none

integer(kind=i2):: bk,nx,ny
integer(kind=i2):: i,j,k
real(kind=dp):: amx,amy,grad,qtemp(ndim+2)
real(kind=dp):: r1,u1,v1,p1,r2,u2,v2,p2
real(kind=dp):: delta_r,delta_u,delta_v,delta_p
k=1

do bk = one,nblk
   call start_end(bk)

i = ie-one
do j = js,je

  amx = 0.5d0*(blk(bk)%zhi(i+1,j,k)%nx+blk(bk)%zhi(i,j,k)%nx)
  amy = 0.5d0*(blk(bk)%zhi(i+1,j,k)%ny+blk(bk)%zhi(i,j,k)%ny)

  grad = dsqrt(amx*amx + amy*amy)
  amx  = amx/grad
  amy  = amy/grad

  qtemp(:)=blk(bk)%fv(i,j,k)%qold(:)
  call q_conv(qtemp)
  r2 = rho
  u2 = u
  v2 = v
  p2 = p

  qtemp(:)=blk(bk)%fv(i,j,k)%qvar(:)
  call q_conv(qtemp)
  r1 = rho
  u1 = u
  v1 = v
  p1 = p

  delta_r = r1 - r2
  delta_u = u1 - u2
  delta_v = v1 - v2
  delta_p = p1 - p2

  delta_w(1,j,2) = -delta_p/(a*a)   + delta_r
  delta_w(1,j,3) =  amy*delta_u     - amx*delta_v
  delta_w(1,j,4) =  delta_p/(rho*a) + amx*delta_u + amy*delta_v

enddo

i = ie
do j = js,je

  amx = 0.5d0*(blk(bk)%zhi(i+1,j,k)%nx+blk(bk)%zhi(i,j,k)%nx)
  amy = 0.5d0*(blk(bk)%zhi(i+1,j,k)%ny+blk(bk)%zhi(i,j,k)%ny)
  grad = dsqrt(amx*amx + amy*amy)
  amx  = amx/grad
  amy  = amy/grad

  qtemp(:)=blk(bk)%fv(i,j,k)%qold(:)
  call q_conv(qtemp)
  r2 = rho
  u2 = u
  v2 = v
  p2 = p

  qtemp(:)=blk(bk)%fv(i,j,k)%qvar(:)
  call q_conv(qtemp)
  r1 = rho
  u1 = u
  v1 = v
  p1 = p

  delta_r = r1 - r2
  delta_u = u1 - u2
  delta_v = v1 - v2
  delta_p = p1 - p2

  delta_w(is,j,2) = -delta_p/(a*a)   + delta_r
  delta_w(is,j,3) =  amy*delta_u     - amx*delta_v
  delta_w(is,j,4) =  delta_p/(rho*a) + amx*delta_u + amy*delta_v

enddo

do j = js,je
  delta_w(3,j,:) = 2*delta_w(is,j,:) - delta_w(1,j,:)
enddo

enddo

!==========================
end subroutine bound_delta
!==========================

!===================================
!subroutine walldatas(cx,cy,mbc)
subroutine walldatas(bk)
!===================================
!calculates the fictitous cell values of conserved variable
!qvar as well as the invisid boundary interface fluxes fin_zhi and
!fin_eta for the complete boundary
!
!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!cx    = same as n for this code
!cy    = same as n for this code
!mbc   = boundary condition vector
!n_zhi = normals on the xi faces
!n_eta = normals on the eta faces
!qvar  = conserved variable q
!fin_zhi = fluxes on the zhi faces (output)
!fin_eta = fluxes on the eta faces (output)
!x = x coordinates
!y = y coordinates
!area = cell areas
!fintr = flux at the given cell interface
!       delta_w = dummy, not used
!===============================================================
use data_kinds
use dims
use flux
use vis
use st_end
implicit none

integer(kind=i2):: idir,iface,iside,ibound
integer(kind=i2):: i,j,k,bk,icell,jcell
real(kind=dp):: q1(ndim+2),q2(ndim+2),q3(ndim+2)!,fintr(ndim+2)
real(kind=dp):: aml(ndim),amr(ndim),amfr(ndim)
real(kind=dp):: tauw,sf,cp,pwall


zhi_wall(:)%tauw(1) = 0.0d0 
zhi_wall(:)%tauw(2) = 0.0d0 
eta_wall(:)%tauw(1) = 0.0d0
eta_wall(:)%tauw(2) = 0.0d0

zhi_wall(:)%pwall(1) = 0.0d0 
zhi_wall(:)%pwall(2) = 0.0d0 
eta_wall(:)%pwall(1) = 0.0d0 
eta_wall(:)%pwall(2) = 0.0d0 

zhi_wall(:)%sf(1) = 0.0d0 
zhi_wall(:)%sf(2) = 0.0d0 
eta_wall(:)%sf(1) = 0.0d0
eta_wall(:)%sf(2) = 0.0d0

zhi_wall(:)%cp(1) = 0.0d0 
zhi_wall(:)%cp(2) = 0.0d0 
eta_wall(:)%cp(1) = 0.0d0
eta_wall(:)%cp(2) = 0.0d0

!xi direction; left boundary
k=1
!do bk = one,nblk
call start_end(bk)

if(blk(bk)%mbc(1)==2) then

idir   = one
iface  = one
iside  = one
ibound = blk(bk)%mbc(1)
icell  = is-one
do j = js,je

  jcell  = j
  aml(1)  = blk(bk)%zhi(icell,j,k)%nx
  aml(2)  = blk(bk)%zhi(icell,j,k)%ny
  amr(1)  = blk(bk)%zhi(icell+1,j,k)%nx
  amr(2)  = blk(bk)%zhi(icell+1,j,k)%ny
  amfr(1)  = blk(bk)%zhi(icell+2,j,k)%nx
  amfr(2)  = blk(bk)%zhi(icell+2,j,k)%ny


  q1(:) = blk(bk)%fv(icell+2,j,k)%qvar(:)
  q2(:) = blk(bk)%fv(icell+1,j,k)%qvar(:)
  q3(:) = blk(bk)%fv(icell,j,k)%qvar(:)

  call wall_data(cx,cy,idir,iside,icell,jcell,aml,amr,amfr,q1,q2,q3) 
                
  eta_wall(j)%pwall(1) = pwall
  eta_wall(j)%tauw(1) = tauw 
  eta_wall(j)%sf(1) = sf 
  eta_wall(j)%cp(1) = cp

enddo
endif

if(blk(bk)%mbc(2)==2) then
!xi direction; right boundary
idir   = one
iface  = two
iside  = two
ibound = blk(bk)%mbc(2)
icell  = ie
do j = js,je

  jcell  = j
  aml(1)  = blk(bk)%zhi(icell,j,k)%nx
  aml(2)  = blk(bk)%zhi(icell,j,k)%ny
  amr(1)  = blk(bk)%zhi(icell+1,j,k)%nx
  amr(2)  = blk(bk)%zhi(icell+1,j,k)%ny
  amfr(1)  = blk(bk)%zhi(icell+2,j,k)%nx
  amfr(2)  = blk(bk)%zhi(icell+2,j,k)%ny


  q1(:) = blk(bk)%fv(icell-1,j,k)%qvar(:)
  q2(:) = blk(bk)%fv(icell,j,k)%qvar(:)
  q3(:) = blk(bk)%fv(icell+1,j,k)%qvar(:)

  call wall_data(cx,cy,idir,iside,icell,jcell,aml,amr,amfr,q1,q2,q3) 

  eta_wall(j)%pwall(2) = pwall
  eta_wall(j)%tauw(2) = tauw 
  eta_wall(j)%sf(2) = sf 
  eta_wall(j)%cp(2) = cp

enddo
endif

if(blk(bk)%mbc(3)==2) then

!eta direction; left boundary
idir   = two
iface  = three
iside  = one
ibound = blk(bk)%mbc(3)
jcell  = js-one
do i = is,ie

  icell  = i
  aml(1)  = blk(bk)%eta(i,jcell,k)%nx
  aml(2)  = blk(bk)%eta(i,jcell,k)%ny
  amr(1)  = blk(bk)%eta(i,jcell+1,k)%nx
  amr(2)  = blk(bk)%eta(i,jcell+1,k)%ny
  amfr(1)  = blk(bk)%eta(i,jcell+2,k)%nx
  amfr(2)  = blk(bk)%eta(i,jcell+2,k)%ny


  q1(:) = blk(bk)%fv(i,jcell+2,k)%qvar(:)
  q2(:) = blk(bk)%fv(i,jcell+1,k)%qvar(:)
  q3(:) = blk(bk)%fv(i,jcell,k)%qvar(:)

  call wall_data(cx,cy,idir,iside,icell,jcell,aml,amr,amfr,q1,q2,q3) 

  zhi_wall(i)%pwall(1) = pwall
  zhi_wall(i)%tauw(1) = tauw 
  zhi_wall(i)%sf(1) = sf 
  zhi_wall(i)%cp(1) = cp

enddo
endif

if(blk(bk)%mbc(4)==2) then

!eta direction; right boundary
idir   = two
iface  = four
iside  = two
ibound = blk(bk)%mbc(4)
jcell  = je
do i = is,ie

  icell  = i
  aml(1)  = blk(bk)%eta(i,jcell,k)%nx
  aml(2)  = blk(bk)%eta(i,jcell,k)%ny
  amr(1)  = blk(bk)%eta(i,jcell+1,k)%nx
  amr(2)  = blk(bk)%eta(i,jcell+1,k)%ny
  amfr(1)  = blk(bk)%eta(i,jcell+2,k)%nx
  amfr(2)  = blk(bk)%eta(i,jcell+2,k)%ny

  q1(:) = blk(bk)%fv(i,jcell-1,k)%qvar(:)
  q2(:) = blk(bk)%fv(i,jcell,k)%qvar(:)
  q3(:) = blk(bk)%fv(i,jcell+1,k)%qvar(:)

  call wall_data(cx,cy,idir,iside,icell,jcell,aml,amr,amfr,q1,q2,q3) 

  zhi_wall(i)%pwall(2) = pwall
  zhi_wall(i)%tauw(2) = tauw 
  zhi_wall(i)%sf(2) = sf 
  zhi_wall(i)%cp(2) = cp

enddo
endif

if(ieuler == 1 ) then
zhi_wall(:)%tauw(1) = 0.0d0 
zhi_wall(:)%tauw(2) = 0.0d0 
eta_wall(:)%tauw(1) = 0.0d0
eta_wall(:)%tauw(2) = 0.0d0
zhi_wall(:)%sf(1) = 0.0d0 
zhi_wall(:)%sf(2) = 0.0d0 
eta_wall(:)%sf(1) = 0.0d0
eta_wall(:)%sf(2) = 0.0d0
endif
!enddo

contains
!============================================================================
subroutine wall_data(cx,cy,idir,iside,icell,jcell,aml,amr,amfr,q1,q2,q3) 
!subroutine wall_data(cx,idir,iside,icell,jcell,aml,amr,q1,q2,q3) 
!============================================================================
use data_kinds
use pri
use inf
use vis
implicit none

integer(kind=i2):: idir,iside,icell,jcell,cx,cy
real(kind=dp):: r1,r2,r3,r4!,area(n,m)
real(kind=dp):: q1(ndim+2),q2(ndim+2),q3(ndim+2)
real(kind=dp):: aml(ndim),amr(ndim),amfr(ndim)
real(kind=dp):: gamma_m1,r,grad_met,amx,amy 
real(kind=dp):: rho1,u1,v1,e1,p1,t1,c1 
real(kind=dp):: rho2,u2,v2,e2,p2,t2,c2 
real(kind=dp):: vel_nor,utan,vtan,veltan
!real(kind=dp):: term4,term5,term6,tauw,sf
real(kind=dp):: gterm2

q3=q3
cp =0.0d0
gterm2 = gamma/(gamma - 1.d0)
gamma_m1 = gamma-one
r        = znd

grad_met = dsqrt(amr(1)*amr(1) + amr(2)*amr(2))
amx      = amr(1)/grad_met
amy      = amr(2)/grad_met

call q_conv(q1)
rho1 = rho
u1   = u
v1   = v
e1   = e
p1   = p
t1   = t
c1   = a

call q_conv(q2)
rho2 = rho
u2   = u
v2   = v
e2   = e
p2   = p
t2   = t
c2   = a

call cell_dist(bk,cx,cy,icell,jcell,idir,iside,r1,r2,r3,r4)

! calculating the normal velocity
!
vel_nor = u2*amx + v2*amy
vel_nor = -vel_nor*(-one)**iside
!
! calculates wall pressure by solving 1d riemann problem
!

!term4 = gamma_m1/2.d0*(vel_nor/c2)
!term5 = 1.d0 - term4
!term6 = term5**(2.d0*sngl(gterm2))
!pwall = p2*term6
!if(pwall .lt. 0.d0) pwall = p2
pwall = p2
!
!Computing wall shear stress
!

utan = U2 - VEL_NOR * AMX
vtan = V2 - VEL_NOR * AMY
veltan = DSIGN(1.d0,utan) * DSQRT(utan*utan + vtan*vtan)


!tauw = mu(icell,jcell) * veltan / R3
tauw = blk(bk)%fv(icell,jcell,k)%mu * veltan / R3
!if(tauw/=0.d0) print*,icell,jcell,tauw
!sf = 0.5 * GRAD_MET * tauw / Re_l 
sf = grad_met * tauw / Re_l 
!write(506,*) X(ICELL,JCELL+1),sf
end subroutine wall_data

!===================
end subroutine walldatas
!===================


!=============================================================
subroutine blockinterface(bk)
!=============================================================
use data_kinds
implicit none

integer(kind=i2):: bk

call block_bc(bk)
return

print*,'Multi block solver not implemented....'
print*
stop
!=============================================================
end subroutine blockinterface
!=============================================================
!
!
!
!=============================================================
subroutine block_conectivity 
!=============================================================
! Storing start & end of neighouring blocks
!                   __________
!                  |f2      f4|
!                  |          |
!                  |   BK4    |
!                  |          |
!                  |          |
!                  |f3______f1|  
!     
!      __________   ____fc4___   _________
!     |f4      f1| |          | |f3     f2|                           
!     |          | |          | |         |                           
!     |   BK1    | fc1  BK   fc2|  BK2    |                            
!     |          | |j         | |         |                           
!     |f2______f3| |__i_fc3___| |f1_____f4| 
!                   __________
!                  |f1      f3|
!                  |          |
!                  |    BK3   |
!                  |          |
!                  |          |
!                  |f4______f2|  
!     
!     ___________________________________________________________
!     ___________________________________________________________
!     Block structure for right hand thumb block orientation  
!     ___________________________________________________________
!     BK is the parent block & BK1 to BK4 are the blocks next
!     to faces  fc1 to fc4 of BK.
!     The index f1 to f4 at corners of neibhoring blocks indicate
!     that if face 1 to 4 of neibhouring blocks are matching then
!     the origin of respective block should be at f1 to f4.
!     ___________________________________________________________
use data_kinds
use flux
use st_end
implicit none

integer(kind=i2) :: bc,bk,kk
integer(kind=i2) :: c1,c2,mesh

mesh=1
c1=1
c2=2
if(mesh==1) then
c1=0
c2=1
endif
!print '(9(A3,2x))','bk','bc','fc','i_1','i_2','j_1','j_2','si','sj'
!print '(45(a1))',('=',i=1,45)

do bk=1,nblk
   
   do bc=1,2*ndim
      kk=blk(bk)%bc(bc) ! block no. next to the "bc" of block bk 
      if(kk==0) cycle   
      
      call start_end(kk)

      !On Face 1 of parent block for face 1 to 4 of neibhoring block 
      if(bc==1) then
      !call blk_nbr_inf(bk,bc,fc,i1,i2,j1,j2,si,sj)
      call blk_nbr_inf(bk,bc,one,is+c2,is+c1,cy   ,one   ,-one,-one)
      call blk_nbr_inf(bk,bc,two,ie-c2,ie-c1,one   ,cy   , one, one)
      call blk_nbr_inf(bk,bc,three,one   ,cx   ,js+c2,js+c1, one,-one)
      call blk_nbr_inf(bk,bc,four,cx   ,one   ,je-c2,je-c1,-one, one)
      endif
      !On Face 2 of parent block for face 1 to 4 of neibhoring block 
      if(bc==2) then
      call blk_nbr_inf(bk,bc,one,is+c1,is+c2,one   ,cy   , one, one)
      call blk_nbr_inf(bk,bc,two,ie-c1,ie-c2,cy   ,one   ,-one,-one)
      call blk_nbr_inf(bk,bc,three,cx   ,one   ,js+c1,js+c2,-one, one)
      call blk_nbr_inf(bk,bc,four,one   ,cx   ,je-c1,je-c2, one,-one)
      endif
      !On Face 3 of parent block for face 1 to 4 of neibhoring block 
      if(bc==3) then
      call blk_nbr_inf(bk,bc,one,is+c2,is+c1,one   ,cy   ,-one, one)
      call blk_nbr_inf(bk,bc,two,ie-c2,ie-c1,cy   ,one   , one,-one)
      call blk_nbr_inf(bk,bc,three,cx   ,one   ,js+c2,js+c1,-one,-one)
      call blk_nbr_inf(bk,bc,four,one   ,cx   ,je-c2,je-c1, one, one)
      endif
      !On Face 4 of parent block for face 1 to 4 of neibhoring block 
      if(bc==4) then
      call blk_nbr_inf(bk,bc,one,is+c1,is+c2,cy   ,one   , one,-one)
      call blk_nbr_inf(bk,bc,two,ie-c1,ie-c2,one   ,cy   ,-one, one)
      call blk_nbr_inf(bk,bc,three,one   ,cx   ,js+c1,js+c2, one, one)
      call blk_nbr_inf(bk,bc,four,cx   ,one   ,je-c1,je-c2,-one,-one)
      endif
    enddo
enddo

contains

subroutine blk_nbr_inf(bk,bc,fc,i_1,i_2,j_1,j_2,si,sj)
use data_kinds
use flux
implicit none

integer(kind=i2) :: bk,bc,fc
integer(kind=i2) :: si,sj    
integer(kind=i2) :: i_1,i_2,j_1,j_2         

blk(bk)%fcs(bc,fc)%is=i_1 ; blk(bk)%fcs(bc,fc)%ie=i_2 ; blk(bk)%fcs(bc,fc)%inc=si  
blk(bk)%fcs(bc,fc)%js=j_1 ; blk(bk)%fcs(bc,fc)%je=j_2 ; blk(bk)%fcs(bc,fc)%jnc=sj
!print '(9(i3,2x))',bk,bc,fc,i_1,i_2,j_1,j_2,si,sj
!print '(6(i3,2x))',i_1,i_2,j_1,j_2,si,sj
end subroutine blk_nbr_inf

!=============================================================
end subroutine block_conectivity 
!=============================================================
!
!
!
!=============================================================
subroutine block_bc(bk)
!=============================================================
use data_kinds
use flux
use st_end
implicit none

integer(kind=i2) :: fc,bc,bk,kk,si,sj,i_1,j_1 
integer(kind=i2) :: wbc(4,4)

!do bk=1,nblk
call start_end(bk)

   ! Fictitious cell index number for copying 
   ! mesh data wrt the final mesh block
   wbc(1,1)=one     ; wbc(1,2)=is-one ; wbc(1,3)=one    ; wbc(1,4)=cy
   wbc(2,1)=ie+one  ; wbc(2,2)=cx     ; wbc(2,3)=one    ; wbc(2,4)=cy
   wbc(3,1)=one     ; wbc(3,2)=cx     ; wbc(3,3)=one    ; wbc(3,4)=js-one
   wbc(4,1)=one     ; wbc(4,2)=cx     ; wbc(4,3)=je+one ; wbc(4,4)=cy  

   do bc=1,2*ndim
      kk=blk(bk)%bc(bc) 
      fc=abs(blk(bk)%fc(bc))
      if(kk==0) cycle   
      !print*,'BK,kk,bc,fc',bk,kk,bc,fc
      !print '("bkk",1x,4(i3,1x))',bk,bc,kk,blk(bk)%fc(bc)
      !print '(4(i3,1x))',wbc(bc,1),wbc(bc,2),wbc(bc,3),wbc(bc,4) 
      !print'(2(i3,1x))', blk(bk)%fcs(bc,fc)%is,blk(bk)%fcs(bc,fc)%js
      !print'(2(i3,1x))', blk(bk)%fcs(bc,fc)%inc,blk(bk)%fcs(bc,fc)%jnc
      !print*
      call blkboundary(bk,kk,bc,fc,wbc) 
   enddo
!enddo


contains

subroutine blkboundary(bk,kk,bc,fc,wbc)
use data_kinds
use flux
implicit none

integer(kind=i2) :: bk,kk,bc,fc
integer(kind=i2) :: i,j,k,si,sj,ii,jj    
integer(kind=i2) :: i_0,i_1,j_0,j_1,is0,js0
integer(kind=i2) :: inc,jnc,wbc(4,4)
k=1
is0=blk(bk)%fcs(bc,fc)%is
js0=blk(bk)%fcs(bc,fc)%js
inc=blk(bk)%fcs(bc,fc)%inc 
jnc=blk(bk)%fcs(bc,fc)%jnc 

i_0=wbc(bc,1)
i_1=wbc(bc,2)
j_0=wbc(bc,3)
j_1=wbc(bc,4)

if(bc==1.or.bc==2) then
!write(56,'("#",8(i3,1x))'),bk,kk,bc,fc,i,j,ii,jj
!write(57,'("#",8(i3,1x))'),bk,kk,bc,fc,i,j,ii,jj
   if(fc==1.or.fc==2) then
!   print*,'#',bk,kk,bc,fc
   ii=is0
   si=inc 
   do i=i_0,i_1
      jj=js0
      sj=jnc 
      do j=j_0,j_1
         blk(bk)%fv(i,j,k)%qvar(1:ndim+2) =blk(kk)%fv(ii,jj,k)%qvar(1:ndim+2)
         blk(bk)%fv(i,j,k)%qold(1:ndim+2) =blk(kk)%fv(ii,jj,k)%qold(1:ndim+2)
         !print '(8(i3,1x))',bk,kk,bc,fc,i,j,ii,jj
         !write(56,*)blk(bk)%mesh(i,j,k)%x,blk(bk)%mesh(i,j,k)%y 
         !write(57,*)blk(kk)%mesh(ii,jj,k)%x,blk(kk)%mesh(ii,jj,k)%y 
         jj=jj+one*sj
      enddo
   ii=ii+one*si
   enddo

   elseif(fc==3.or.fc==4) then
  ! print*,bk,kk,bc,fc
   jj=js0
   sj=jnc 
   do i=i_0,i_1
      ii=is0
      si=inc 
      do j=j_0,j_1
         blk(bk)%fv(i,j,k)%qvar(1:ndim+2) =blk(kk)%fv(ii,jj,k)%qvar(1:ndim+2)
         blk(bk)%fv(i,j,k)%qold(1:ndim+2) =blk(kk)%fv(ii,jj,k)%qold(1:ndim+2)
         !print '(8(i3,1x))',bk,kk,bc,fc,i,j,ii,jj
         !write(56,*)blk(bk)%mesh(i,j,k)%x,blk(bk)%mesh(i,j,k)%y 
         !write(57,*)blk(kk)%mesh(ii,jj,k)%x,blk(kk)%mesh(ii,jj,k)%y 
         ii=ii+one*si
      enddo
   jj=jj+one*sj
   enddo
   endif
elseif(bc==3.or.bc==4) then

   if(fc==1.or.fc==2) then
   print*,bk,kk,bc,fc
   ii=is0
   si=inc
   do j=j_0,j_1   
      jj=js0
      sj=jnc
      do i=i_0,i_1
         blk(bk)%fv(i,j,k)%qvar(1:ndim+2) =blk(kk)%fv(ii,jj,k)%qvar(1:ndim+2)
         blk(bk)%fv(i,j,k)%qold(1:ndim+2) =blk(kk)%fv(ii,jj,k)%qold(1:ndim+2)
         !print '(8(i3,1x))',bk,kk,bc,fc,i,j,ii,jj
         !write(56,*)blk(bk)%mesh(i,j,k)%x,blk(bk)%mesh(i,j,k)%y 
         !write(57,*)blk(kk)%mesh(ii,jj,k)%x,blk(kk)%mesh(ii,jj,k)%y 
         jj=jj+one*sj
      enddo
      ii=ii+one*si
   enddo
   elseif(fc==3.or.fc==4) then
   print*,bk,kk,bc,fc
   jj=js0
   sj=jnc
   do j=j_0,j_1   
      ii=is0
      si=inc
      do i=i_0,i_1
         blk(bk)%fv(i,j,k)%qvar(1:ndim+2) =blk(kk)%fv(ii,jj,k)%qvar(1:ndim+2)
         blk(bk)%fv(i,j,k)%qold(1:ndim+2) =blk(kk)%fv(ii,jj,k)%qold(1:ndim+2)
         !print '(8(i3,1x))',bk,kk,bc,fc,i,j,ii,jj
         !write(56,*)blk(bk)%mesh(i,j,k)%x,blk(bk)%mesh(i,j,k)%y 
         !write(57,*)blk(kk)%mesh(ii,jj,k)%x,blk(kk)%mesh(ii,jj,k)%y 
         ii=ii+one*si
      enddo
      jj=jj+one*sj
   enddo
   endif
else
   print*
   print*,'face no. :',bc
   print*,'Face number should be between 1 to 4'
   print*
   stop
endif


end subroutine blkboundary 
!=============================================================
end subroutine block_bc
!=============================================================
