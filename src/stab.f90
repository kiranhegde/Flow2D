! stab.f.old
!=================================================================
subroutine stability(bk,cfl)
!=================================================================

!calculates the delta-t (local time-stepping)

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!nx    = same as n for this code
!ny    = same as n for this code
!cfl   = cfl no.
!n_zhi = normals on the xi faces
!n_eta = normals on the eta faces
!qvar  = q in the domain
!dt    = delta-t, local time-stepping
!area  = area
!===============================================================
use st_end
use pri
use flux
use vis
use inf
implicit none

integer(kind=i2):: bk
integer(kind=i2):: i,j,k,ii
real(kind=dp):: cfl
real(kind=dp):: cfl_zhi,cfl_eta,dt_min,area,dt
real(kind=dp):: zhixav,zhiyav
real(kind=dp):: etaxav,etayav,conv_spec,vis_spec
real(kind=dp):: cvelx,cvely,term1,term2,dt_conv
real(kind=dp):: q1(ndim+2)
real(kind=dp):: dsi,dsj,vol 
!real(kind=dp):: c_fac,dtv,mu,mut,f1,f2,fmue,fac

cfl_zhi  = 0.0_dp
cfl_eta  = 0.0_dp
dt_min   = 1.d6
!c_fac=2.0

k=1
!do bk=1,nblk
   call start_end(bk)

!!$OMP PARALLEL DO DEFAULT(NONE) &
!!$OMP FIRSTPRIVATE(bk,k) &
!!$OMP PRIVATE(i,j,term_conv,zhixav,zhiyav,etaxav,etayav, &
!!$OMP cvelx,cvely,u,v,a) &
!!$OMP SHARED(blk,is,ie,js,je,cfl,q1)
!!$omp parallel do
!do j = js,je 
!  do i = is,ie
do j =1,cy 
  do i = 1,cx
    zhixav = 0.5_dp*(blk(bk)%zhi(i+1,j,k)%nx+blk(bk)%zhi(i,j  ,k)%nx)
    zhiyav = 0.5_dp*(blk(bk)%zhi(i+1,j,k)%ny+blk(bk)%zhi(i,j  ,k)%ny)
    etaxav = 0.5_dp*(blk(bk)%eta(i,j+1,k)%nx+blk(bk)%eta(i,j  ,k)%nx)
    etayav = 0.5_dp*(blk(bk)%eta(i,j+1,k)%ny+blk(bk)%eta(i,j  ,k)%ny)

    do ii=1,ndim+2
       q1(ii) = blk(bk)%fv(i,j,k)%qvar(ii) 
    enddo
    call q_conv(q1)

    cvelx   = u*zhixav  + v*zhiyav
    cvely   = u*etaxav  + v*etayav

    dsi= zhixav*zhixav+zhiyav*zhiyav 
    dsj= etaxav*etaxav+etayav*etayav

    conv_spec = dabs(cvelx) + a*dsqrt(dsi) + &
                dabs(cvely) + a*dsqrt(dsj)

    vol=blk(bk)%fv(i,j,k)%vol
   
    !mu=blk(bk)%fv(i,j,k)%mu
    !mut=blk(bk)%fv(i,j,k)%mut
    !fmue = (mu+mut)
    !f1   = 4.0/3.0
    !f2   = gamma/(Pr_l+Pr_t)
    !f2   = gamma/Pr_l
    !fac  = dMAX1(f1,f2)/rho
    !vis_spec=fac*fmue*(dsi+dsj)*M_inf/Re_L/vol
    !print*,conv_spec,c_fac*vis_spec
    !blk(bk)%fv(i,j,k)%dt  = cfl*vol/(conv_spec+c_fac*vis_spec)
    
    blk(bk)%fv(i,j,k)%dt  = cfl*vol/conv_spec

!    area=blk(bk)%fv(i,j,k)%vol
!    dt=blk(bk)%fv(i,j,k)%dt

!    cfl_zhi  = dmax1(cfl_zhi,dt*term1/area)
!    cfl_eta  = dmax1(cfl_eta,dt*term2/area)
!    dt_min   = dmin1(dt_min,dt)
    !blk(bk)%fv(i,j,k)%dt  = dt_min 
  enddo
enddo

!!$OMP END PARALLEL DO

!write(*,*) 'cfl_zhi  = ', cfl_zhi
!write(*,*) 'cfl_eta  = ', cfl_eta

!========================
end subroutine stability
