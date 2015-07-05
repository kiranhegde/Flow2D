! im_fl5.f.old
!===================================================
subroutine implicit_solver(bk)
!===================================================

!implicit time-stepping using matrix-free lu-ssor (lu-sgs)

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!nx    = same as n for this code
!ny    = same as n for this code
!qvar  = conserved variable q
!n_zhi = normals on the xi faces
!n_eta = normals on the eta faces
!dt    = delta-t (local time stepping)
!area = cell areas
!x = x coordinates
!y = y coordinates
!del_flux = total flux in a cell (sum of fin_zhi and fin_eta)
!del_q = incremental change in qvar
!===============================================================
use st_end
use flux
use pri
implicit none

integer(kind=i2):: bk,kstart,kend
real(kind=dp):: q_p

integer(kind=i2):: inner2,iter_inner2,jterm,i,j,k,ii
real(kind=dp):: alpha_lu,beta_lu,kappa_lu,epsilon
real(kind=dp):: amr(ndim),ev(ndim+2),qpos(ndim+2),qneg(ndim+2)
real(kind=dp):: fpos(ndim+2),fneg(ndim+2)
real(kind=dp):: grad_met,contra_vel,acoustic,spec_rad,dtbarea1
real(kind=dp):: lower,upper

type lussors
     real(kind=dp):: ap,ra,rb,dt_area
end type lussors     

type fvms
       real(kind=dp),dimension(ndim+2):: qold,dq!,qtmp
       real(kind=dp),dimension(ndim+2,ndim):: qflux
end type fvms

type(lussors),allocatable,dimension(:,:,:) :: lu 
type(fvms),allocatable,dimension(:,:,:):: fvm

allocatable :: q_p(:)


allocate(q_p(ndim+2))

k=1
call start_end(bk)

allocate(lu(cx,cy,cz),fvm(cx,cy,cz))

!call cpu_time(time_begin)
iter_inner2 = 10
alpha_lu    = 1.5_dp
beta_lu     = 0.5_dp
kappa_lu    = 2.01_dp
epsilon     = 1.d-4

!!$OMP PARALLEL

!!$OMP DO

do j=1,cy
do i=1,cx
   lu(i,j,k)%ap = 0.0_dp
   lu(i,j,k)%ra = 0.0_dp
   lu(i,j,k)%rb = 0.0_dp
   lu(i,j,k)%dt_area = 0.0_dp
   fvm(i,j,k)%dq(1:ndim+2) = 0.0_dp
   fvm(i,j,k)%qflux(1:ndim+2,:) = 0.0_dp
!   fvm(i,j,k)%qtmp(1:ndim+2) =blk(bk)%fv(i,j,k)%qvar(1:ndim+2) 
   fvm(i,j,k)%qold(1:ndim+2) = blk(bk)%fv(i,j,k)%qvar(1:ndim+2)
enddo
enddo
!!$OMP END DO



!calculate the spectral radius for zhi direction

!!$OMP DO PRIVATE(amr,q_p)
do j =js-one,je+one
  do i =is-one,ie+one
     amr(1) = 0.5_dp*(blk(bk)%zhi(i,j,k)%nx + blk(bk)%zhi(i+1,j,k)%nx)
     amr(2) = 0.5_dp*(blk(bk)%zhi(i,j,k)%ny + blk(bk)%zhi(i+1,j,k)%ny)
     q_p(:) = blk(bk)%fv(i,j,k)%qvar(:)
     !call cell_dist(bk,cx,cy,i,j,1,2,r1,r2,r3,r4) 
     lu(i,j,k)%ra = Spec_radius(amr,q_p)!,r2+r3)
  enddo
enddo
!!$OMP END DO

!calculate the spectral radius for eta direction

!!$OMP DO PRIVATE(amr,q_p)
do i =is-one,ie+one
  do j =js-one,je+one
     amr(1) = 0.5_dp*(blk(bk)%eta(i,j,k)%nx + blk(bk)%eta(i,j+1,k)%nx)
     amr(2) = 0.5_dp*(blk(bk)%eta(i,j,k)%ny + blk(bk)%eta(i,j+1,k)%ny)
     q_p(:) = blk(bk)%fv(i,j,k)%qvar(:)
     !call cell_dist(bk,cx,cy,i,j,2,2,r1,r2,r3,r4)
     lu(i,j,k)%rb = Spec_radius(amr,q_p)!,r2+r3)
  enddo
enddo
!!$OMP END DO

k=1
!calculate the diagonal element

!!$OMP DO 
do j =js,je
  do i =is,ie
     lu(i,j,k)%dt_area = blk(bk)%fv(i,j,k)%dt/blk(bk)%fv(i,j,k)%vol
     lu(i,j,k)%ap = 1 + alpha_lu*lu(i,j,k)%dt_area*(lu(i,j,k)%ra + lu(i,j,k)%rb)
  enddo
enddo
!!$OMP END DO

!!$OMP END PARALLEL


do inner2 = 1,iter_inner2


!!$OMP PARALLEL

!!$OMP DO PRIVATE(dtbarea1,amr,qpos,qneg)

!forward sweep
!(l+d)*del_q = -r

  do j = js,je
    do i = is,ie

       !calculate the diagonal element
       !lu(i,j,k)%dt_area = blk(bk)%fv(i,j,k)%dt/blk(bk)%fv(i,j,k)%vol
       !lu(i,j,k)%ap = 1 + alpha_lu*lu(i,j,k)%dt_area*(lu(i,j,k)%ra + lu(i,j,k)%rb)
       dtbarea1 = 0.5_dp*lu(i,j,k)%dt_area*alpha_lu

      do jterm = 1,ndim+2
         fvm(i,j,k)%dq(jterm) = rhs(i,j,jterm) 
      enddo

      !blk(bk)%fv(i,j,k)%qvar(:) = fvm(i,j,k)%qtmp(:) + fvm(i,j,k)%dq(:)
      blk(bk)%fv(i,j,k)%qvar(:) = fvm(i,j,k)%qold(:) + fvm(i,j,k)%dq(:)

      amr(1) = blk(bk)%zhi(i,j,k)%nx
      amr(2) = blk(bk)%zhi(i,j,k)%ny

      qpos(:) = blk(bk)%fv(i,j,k)%qvar(:)
      qneg(:) = fvm(i,j,k)%qold(:)
      call flux_pn
      fvm(i,j,k)%qflux(:,1) = fpos(:) - fneg(:)

      amr(1) = blk(bk)%eta(i,j,k)%nx
      amr(2) = blk(bk)%eta(i,j,k)%ny

      qpos(:) = blk(bk)%fv(i,j,k)%qvar(:) 
      qneg(:) = fvm(i,j,k)%qold(:)
      call flux_pn
      fvm(i,j,k)%qflux(:,2) = fpos(:) - fneg(:)

    enddo
  enddo
!  enddo
!!$OMP END DO

!backward sweep
!(u+d)*del_q = d*del_q


!!$OMP DO PRIVATE(dtbarea1,amr,qpos,qneg)
  do j = je,js,-one
    do i = ie,is,-one

      dtbarea1 = 0.5_dp*lu(i,j,k)%dt_area*alpha_lu

      do jterm = 1,ndim+2
         fvm(i,j,k)%dq(jterm) = rhs(i,j,jterm) 
      enddo

      !blk(bk)%fv(i,j,k)%qvar(:) = fvm(i,j,k)%qtmp(:) + fvm(i,j,k)%dq(:) 
      blk(bk)%fv(i,j,k)%qvar(:) = fvm(i,j,k)%qold(:) + fvm(i,j,k)%dq(:) 

      amr(1) = blk(bk)%zhi(i,j,k)%nx
      amr(2) = blk(bk)%zhi(i,j,k)%ny

      qpos(:) = blk(bk)%fv(i,j,k)%qvar(:)
      qneg(:) = fvm(i,j,k)%qold(:) 
      call flux_pn
      fvm(i,j,k)%qflux(:,1) = fpos(:) - fneg(:)

      amr(1) = blk(bk)%eta(i,j,k)%nx
      amr(2) = blk(bk)%eta(i,j,k)%ny

      qpos(:) = blk(bk)%fv(i,j,k)%qvar(:)
      qneg(:) = fvm(i,j,k)%qold(:)
      call flux_pn
      fvm(i,j,k)%qflux(:,2) = fpos(:) - fneg(:)

    enddo
  enddo
!  enddo
!!$OMP END DO

!!$OMP END PARALLEL

enddo ! inner2

!!$OMP PARALLEL

!!$OMP DO
do j=1,cy
do i=1,cx
   blk(bk)%fv(i,j,k)%qvar(:) = fvm(i,j,k)%qold(:)
enddo
enddo
!!$OMP END DO

!...........updating q (conserved variable)
!!$OMP DO
do j = js,je
   do i = is,ie
       blk(bk)%fv(i,j,k)%qvar(:) = blk(bk)%fv(i,j,k)%qold(:) + fvm(i,j,k)%dq(:)
   enddo
enddo
!!$OMP END DO

!!$OMP END PARALLEL

deallocate(lu,fvm)
call flush(6)

!enddo
deallocate(q_p)

contains

!================================================
subroutine flux_pn
!================================================

!flux for a given q

!===============================================================
!ndim  = space dimension (2)
!amr   = normals
!qpos,qneg = q values
!fpos,fneg = corresponing analytical value of flux
!===============================================================
real(kind=dp):: contra_vel

call q_conv(qpos)
contra_vel = u*amr(1) + v*amr(2)
fpos(1)    = contra_vel*(e + p)
fpos(2)    = contra_vel*rho
fpos(3)    = fpos(2)*u + p*amr(1)
fpos(4)    = fpos(2)*v + p*amr(2)

call q_conv(qneg)
contra_vel = u*amr(1) + v*amr(2)
fneg(1)    = contra_vel*(e + p)
fneg(2)    = contra_vel*rho
fneg(3)    = fneg(2)*u + p*amr(1)
fneg(4)    = fneg(2)*v + p*amr(2)

!======================
end subroutine flux_pn
!======================

!====================================
real(kind=DP) function rhs(i,j,jterm)
!====================================
implicit none

integer(kind=i2):: i,j,jterm

lower =  dtbarea1*(fvm(i-1,j,k)%qflux(jterm,1)+fvm(i,j-1,k)%qflux(jterm,2) &
              + lu(i-1,j,k)%ra*fvm(i-1,j,k)%dq(jterm) &
              + lu(i,j-1,k)%rb*fvm(i,j-1,k)%dq(jterm))

upper = -dtbarea1*(fvm(i+1,j,k)%qflux(jterm,1)+fvm(i,j+1,k)%qflux(jterm,2) &
              - lu(i+1,j,k)%ra*fvm(i+1,j,k)%dq(jterm) &
              - lu(i,j+1,k)%rb*fvm(i,j+1,k)%dq(jterm))
q_p(jterm) = -beta_lu*lu(i,j,k)%dt_area*blk(bk)%fv(i,j,k)%dflux(jterm)+lower+upper 
       rhs =  q_p(jterm)/lu(i,j,k)%ap

!===============
end function rhs 
!===============

!==========================================
real(kind=DP) function Spec_radius(amr,q_p)
!==========================================
implicit none

real(kind=dp):: q_p(ndim+2),amr(ndim)

    call q_conv(q_p)
    grad_met = dsqrt(amr(1)*amr(1) + amr(2)*amr(2))
    contra_vel = u*amr(1) + v*amr(2)
    acoustic   = a*grad_met
    ev(1)      = contra_vel - acoustic
    ev(2)      = contra_vel
    ev(3)      = contra_vel
    ev(4)      = contra_vel + acoustic
    spec_rad        = dabs(contra_vel) + acoustic + epsilon
    Spec_radius  = kappa_lu*spec_rad

!=======================
end function Spec_radius
!=======================



!==============================
end subroutine implicit_solver
!==============================
