!====================================
subroutine TurbulenceModel(bk)
!====================================
use dims
use flux
use vis
use st_end
implicit none

integer(kind=i2):: mbc(2*ndim),bk,nx,ny,k
k=1
if ( turb == 0) then
  ! visturb(:,:)=0.0d0
  blk(bk)%fv(:,:,:)%mut=0.0d0
  return
else
   call start_end(bk)
   call BLomax(bk,cy)
endif
!=============================
end subroutine TurbulenceModel
!=============================

!=============================
subroutine BLomax(bk,ny)
!=============================
use dims
use vis
use flux
use pri
use st_end
implicit none

integer(kind=i2):: i,j,k,bk,nx,ny,jj
real(kind=dp):: qtmp(ndim+2),tmp,amr(ndim+2),aml(ndim+2),qleft(ndim+2),qright(ndim+2)
real(kind=dp):: yy(ny),rh(ny),uu(ny),vo(ny),vi(ny),tauw,tv(ny)
real(kind=dp) :: u_x,u_y,v_x,v_y,t_x,t_y


k=1
call start_end(bk)
do i = is,ie!-1
   do j = js,je

      aml(1)  = blk(bk)%zhi(i,j,k)%nx
      aml(2)  = blk(bk)%zhi(i,j,k)%ny
      amr(1)  = blk(bk)%zhi(i+1,j,k)%nx
      amr(2)  = blk(bk)%zhi(i+1,j,k)%ny

      qleft(:)   = blk(bk)%fv(i,j,k)%qvar(:)
      qright(:)  = blk(bk)%fv(i+1,j,k)%qvar(:)
      !call deriv_int_xi (bk,i,j,1,qleft,qright,amr,aml) 
      u_x=blk(bk)%zhi(i,j,k)%u_x
      u_y=blk(bk)%zhi(i,j,k)%u_y
      v_x=blk(bk)%zhi(i,j,k)%v_x
      v_y=blk(bk)%zhi(i,j,k)%v_y
      t_x=blk(bk)%zhi(i,j,k)%t_x
      t_y=blk(bk)%zhi(i,j,k)%t_y 

      Blmx(i,j,1)%vo = dabs(v_x-u_y)
   enddo      
enddo      

do j = js,je-1
   do i = is,ie
      aml(1)  = blk(bk)%eta(i,j,k)%nx
      aml(2)  = blk(bk)%eta(i,j,k)%ny
      amr(1)  = blk(bk)%eta(i,j+1,k)%nx
      amr(2)  = blk(bk)%eta(i,j+1,k)%ny

      qleft(:)   = blk(bk)%fv(i,j,k)%qvar(:)
      qright(:)  = blk(bk)%fv(i,j+1,k)%qvar(:)
      !call deriv_int_eta(bk,i,j,2,qleft,qright,amr,aml)
      u_x=blk(bk)%eta(i,j,k)%u_x
      u_y=blk(bk)%eta(i,j,k)%u_y
      v_x=blk(bk)%eta(i,j,k)%v_x
      v_y=blk(bk)%eta(i,j,k)%v_y
      t_x=blk(bk)%eta(i,j,k)%t_x
      t_y=blk(bk)%eta(i,j,k)%t_y

      Blmx(i,j,2)%vo = dabs(v_x-u_y) 
   enddo
enddo

do i = is,ie
   do j = js,je

      qtmp(:) = blk(bk)%fv(i,j,k)%qvar(:) 
      call q_conv(qtmp) 
      Blmx(i,j,:)%rh = rho 
      Blmx(i,j,:)%uu = dsqrt(u*u+v*v) 
                 tmp =dsqrt(Blmx(i,j,1)%vo*Blmx(i,j,1)%vo+Blmx(i,j,2)%vo*Blmx(i,j,2)%vo) 
      Blmx(i,j,:)%vo = tmp
      Blmx(i,j,:)%vi = blk(bk)%fv(i,j,k)%mu 
      Blmx(i,j,:)%tauw = blk(bk)%fv(1,j,k)%mu*tmp 
         
      yy(j) = Blmx(i,j,3)%yy
      rh(j) = rho
      uu(j) = dsqrt(u*u+v*v) 
      vo(j) = tmp 
      !vi(j) = mu(i,j)
      vi(j) = blk(bk)%fv(i,j,k)%mu 
   enddo
      jj=js
      !tauw  =mu(i,jj)*vo(jj)! blk(bk)%fv(i,j,k)%mu
      !tauw  = zhi_wall(i)%tauw(1)
      call BaldwinLomax(ny,yy,rh,uu,vo,vi,tauw,tv)
!      print*,"i,            mu,     vo            tv"
!      print*,i,mu(i,jj),vo(jj),tv(jj+4) 
!      pause
      !visturb(i,:)=tv(:)
      blk(bk)%fv(i,:,k)%mut=tv(:)
      !visturb(i,:)=0.0d0
enddo

!do i=1,cx
!   write(21,*)i, blk(bk)%fv(i,4,k)%mut, blk(bk)%fv(i,5,k)%mut, blk(bk)%fv(i,6,k)%mut
!enddo
!stop

!====================
end subroutine BLomax
!====================

!======================================
subroutine walldistance(bk)
!======================================
use flux
use dims
use vis
use st_end

implicit none

integer(kind=i2):: i,j,k,bk
real(kind=dp):: d1,d2,dy!,area(n,m)
real(kind=dp) :: x1,x2,y1,y2,dx1,dy1,dd,ds,nx1,ny1

!.......Computing  distnce from the  wall

!.......zhi direction
k=1
! face 1
if(Blmx(1,1,1)%bc==2) then
!do bk=1,nblk
!call start_end(bk)
do i=1,cx
do j=1,cy
   dx1 = blk(bk)%eta(i,j,k)%nx
   dy1 = blk(bk)%eta(i,j,k)%ny
   d1 = dsqrt(dx1*dx1+dy1*dy1)

   dx1 = blk(bk)%eta(i,j+1,k)%nx
   dy1 = blk(bk)%eta(i,j+1,k)%ny
   d2 = dsqrt(dx1*dx1+dy1*dy1)

   dy = 2.0d0*blk(bk)%fv(i,j,k)%vol/(d1+d2)

   Blmx(i,j,1)%yy=Blmx(i,j,1)%yy+0.5d0*dy  
   Blmx(i+1,j,1)%yy=Blmx(i,j,1)%yy+0.5d0*dy  
enddo
enddo
!enddo
endif



! face 2
if(Blmx(1,1,2)%bc==2) then
!do bk=1,nblk
!call start_end(bk)
do i=cx,1,-1
do j=cy,1,-1
   dx1 = blk(bk)%eta(i,j,k)%nx
   dy1 = blk(bk)%eta(i,j,k)%ny
   d1 = dsqrt(dx1*dx1+dy1*dy1)

   dx1 = blk(bk)%eta(i,j+1,k)%nx
   dy1 = blk(bk)%eta(i,j+1,k)%ny
   d2 = dsqrt(dx1*dx1+dy1*dy1)
   dy = 2.0d0*blk(bk)%fv(i,j,k)%vol/(d1+d2)

   Blmx(i-1,j,2)%yy=Blmx(i,j,2)%yy+0.5d0*dy  
   Blmx(i,j,2)%yy=Blmx(i,j,2)%yy+0.5d0*dy  
enddo
enddo
!enddo
endif

! face 3
if(Blmx(1,1,3)%bc==2) then
!do bk=1,nblk
!call start_end(bk)
do j=1,cy
do i=1,cx
   dx1 = blk(bk)%zhi(i,j,k)%nx
   dy1 = blk(bk)%zhi(i,j,k)%ny
   d1 = dsqrt(dx1*dx1+dy1*dy1)

   dx1 = blk(bk)%zhi(i+1,j,k)%nx
   dy1 = blk(bk)%zhi(i+1,j,k)%ny
   d2 = dsqrt(dx1*dx1+dy1*dy1)
   dy = 2.0d0*blk(bk)%fv(i,j,k)%vol/(d1+d2)
   Blmx(i,j,3)%yy=Blmx(i,j,3)%yy+0.5d0*dy  
   Blmx(i,j+1,3)%yy=Blmx(i,j,3)%yy+0.5d0*dy  
enddo
enddo
!enddo
endif

! face 4
if(Blmx(1,1,4)%bc==2) then
!do bk=1,nblk
!call start_end(bk)
do j=cy,1,-1
do i=cx,1,-1
   dx1 = blk(bk)%zhi(i,j,k)%nx
   dy1 = blk(bk)%zhi(i,j,k)%ny
   d1 = dsqrt(dx1*dx1+dy1*dy1)

   dx1 = blk(bk)%zhi(i+1,j,k)%nx
   dy1 = blk(bk)%zhi(i+1,j,k)%ny
   d2 = dsqrt(dx1*dx1+dy1*dy1)
   dy = 2.0d0*blk(bk)%fv(i,j,k)%vol/(d1+d2)

   Blmx(i,j-1,4)%yy=Blmx(i,j,4)%yy+0.5d0*dy  
   Blmx(i,j,4)%yy=Blmx(i,j,4)%yy+0.5d0*dy  
enddo
enddo
!enddo
endif


!do j=1,m
!do i=1,n
!   write(13,*)x(i,j),y(i,j)
!enddo
!   write(13,*)
!enddo
!
!do i=1,n
!do j=1,m
!   write(13,*)x(i,j),y(i,j)
!enddo
!   write(13,*)
!enddo
!
!print*,nx,ny
!print*,n,m
!
!!do i=1,n
!do i=is,ie
!do j=js,js+1
!   x1=x(i,j) ; y1 = y(i,j) 
!   x2=x(i+1,j) ; y2 = y(i+1,j) 
!   x1=0.5d0*(x1+x2)
!   y1=0.5d0*(y1+y2)
!   write(14,*)x1,y1
!   dx1=x2-x1 ; dy1 = y2-y1
!   ds=dsqrt(dx1*dx1+dy1*dy1)
!   nx1=-dy1 ; ny1 = dx1
!   dd = Blmx(i,j,3)%yy  
!
!   !write(14,*)x1+nx1*dd/ds,y1+ny1*dd/ds
!   write(14,*)x1+nx1*ds/dd,y1+ny1*ds/dd
!   write(14,*)
!enddo
!   write(14,*)
!enddo

!==========================
end subroutine walldistance 
!==========================


!=====================================================================
! Subroutine to compute turbulence viscosity using Baldwin-Lomax model
!=====================================================================
subroutine BaldwinLomax(m,yy,rh,uu,vo,vi,tauw,tv) 
!=====================================================================
!
!-----Purpose:  This subroutine computes the turbulent viscosity
!               coefficient using the Baldwin-Lomax model.
!

!-----Use modules
Use dims   ! kind parameters
use vis
use st_end
!-----Require explicit typing of variables
Implicit none

!-----Parameter statements
!Include 'mxdim.par'

!-----Common blocks
!Include 'test.inc'

!-----Input arguments
!Integer(kind=i2), intent(in) :: jedge       ! Index of bl edge
Integer(kind=i2), intent(in) :: m       ! Index of bl edge

!Real(kind=dp), intent(in) :: re_l          ! Reference Reynolds number
Real(kind=dp), intent(in) :: yy(M)   ! Distance from wall
Real(kind=dp), intent(in) :: rh(M)   ! Static density
Real(kind=dp), intent(in) :: uu(m)   ! Velocity
Real(kind=dp), intent(in) :: vo(M)   ! Vorticity
Real(kind=dp), intent(in) :: vi(M)   ! Laminar viscosity coeff
!Real(kind=dp), intent(in) :: tauw        ! Shear stress at the wall
Real(kind=dp) :: tauw        ! Shear stress at the wall

!-----Output arguments
!Real(kind=dp), intent(out) :: tv(jedge)  ! Turbulent viscosity coeff
Real(kind=dp), intent(out) :: tv(m)  ! Turbulent viscosity coeff

!-----Local variables
Integer(kind=i2) icross   ! Flag for inner/outer region
Integer(kind=i2) j        ! Do loop index
Integer(kind=i2) jedg     ! Index limit for Fmax search
Integer(kind=i2) itest1,itest2,itest3                   

Real(kind=dp) al       ! Mixing length
Real(kind=dp) aplus    ! Van Driest damping constant
Real(kind=dp) bigk     ! Clauser constant
Real(kind=dp) ccp      ! Constant in outer region model
Real(kind=dp) ckleb    ! Constant in Klebanoff intermittency factor
Real(kind=dp) cwk      ! Constant in outer region model
Real(kind=dp) fkleb    ! Klebanoff intermittency factor
Real(kind=dp) fl       ! Baldwin-Lomax F parameter
Real(kind=dp) fmax     ! Baldwin-Lomax Fmax parameter
Real(kind=dp) frac     ! Fractional decrease in F req'd for peak
Real(kind=dp) fwake    ! Baldwin-Lomax Fwake parameter
Real(kind=dp) rdum     ! Ratio of distance from wall to ymax
Real(kind=dp) smlk     ! Von Karman constant
Real(kind=dp) tvi      ! Inner region turbulent viscosity coeff
Real(kind=dp) tvo      ! Outer region turbulent viscosity coeff
Real(kind=dp) udif     ! Max velocity difference
Real(kind=dp) umax     ! Max velocity
Real(kind=dp) umin     ! Min velocity
Real(kind=dp) ymax     ! Distance from wall to location of Fmax
Real(kind=dp) yp       ! y+
Real(kind=dp) ypcon    ! Coeff term for y+, based on wall values
Real(kind=dp) ypconl   ! Coeff term for y+
Real(kind=dp) yyj      ! Distance from wall

!-----Set constants

aplus = 26.d0
ccp   = 1.6d0
ckleb = 0.3d0
cwk   = 0.25d0
smlk  = 0.4d0
bigk  = 0.0168d0

itest1 = 0
itest2 = 1
itest3 = 0

If (itest1 == 1) bigk = 0.0180   ! Comp. correction (cfl3de)

!-----Compute stuff needed in model

!-----Get coefficient term for y+
If (itest2 == 1) then             ! Using wall vorticity as in cfl3de
   ypcon = Sqrt(re_l*rh(js)*vo(js)/vi(js))
Else                                 ! Using wall shear stress
   If (tauw <= 1.d-9) tauw = 1.d-9
      ypcon = Sqrt(re_l*rh(js)*tauw)/vi(js)
   End if

!-----Set index limit for Fmax search, and fractional decrease needed to
!-----   qualify as first peak
jedg = je

frac = .70d0
If (itest3 > 0) then           ! User-spec. frac. decrease
frac = dble(itest3)/1000.d0
Else if (itest3 < 0) then      ! Reset search range, use max
jedg = Min(jedg,-itest3)   !    value, not first peak
frac = 0.0d0
Endif

!-----Get max velocity and max velocity difference
umin = 0.d0
umax = 0.d0
!Do j = 2,m
Do j = js,jedg
umax = dMax1(umax,uu(j))
umin = dMin1(umin,uu(j))
End do
udif = umax - umin

!-----Get Fmax by searching for first peak in F

ymax = 0.d0
fmax = 0.d0
ypconl = ypcon
!Do j = 2,jedg
Do j = js,jedg
yyj = yy(j)
If (itest2 == 1) then               ! Use local values in y+
  ypconl = ypcon*Sqrt(rh(j)/rh(js))*vi(js)/vi(j)
End if
yp = ypconl*yyj                        ! y+
fl = yyj*vo(j)*(1.d0 - Exp(-yp/aplus))   ! B-L F parameter
If (fl > fmax) then                    ! Set new Fmax
  fmax = fl
  ymax = yyj
Else if (fl > frac*fmax) then          ! Keep searching
  Cycle
Else                                   ! Found Fmax, so get out
  Exit
End if
End do

!-----Reset ymax and Fmax if necessary, to avoid overflows
If (ymax < 1.d-6) ymax = 5.d-5
If (fmax < 1.d-6) fmax = 5.d-5

!-----Compute turbulent viscosity

icross = 0
ypconl = ypcon
!Do j = 2,jedge
!Do j = 2,m
Do j = js,jedg
yyj = yy(j)
tvi = 1.d10

!--------Inner region value, if we're still there
If (icross == 0) then
  If (itest2 == 1) then             ! Use local values in y+
     ypconl = ypcon*Sqrt(rh(j)/rh(js))*vi(js)/vi(j)
  End if
  yp = ypconl*yyj                      ! y+
  al = smlk*yyj*(1.d0 - Exp(-yp/aplus))  ! Mixing length
  tvi = rh(j)*al*al*vo(j)
End if

!--------Outer region value
rdum = yyj/ymax
If (rdum >= 1.d5) then                   ! Prevent overflow
 fkleb = 0.0d0
Else                                     ! Klebanoff intermittency factor
 fkleb = 1.d0/(1.d0 + 5.5d0*(ckleb*rdum)**6)
End if
fwake = dMin1(ymax*fmax,cwk*ymax*udif*udif/fmax)
tvo = bigk*ccp*rh(j)*fwake*fkleb

!--------Set turbulent viscosity, plus flag if we're in outer region
tv(j) = tvi
If (tvo < tvi) then
  icross = 1
  tv(j) = tvo
End if

!--------Non-dimensionalize
tv(j) = re_l*tv(j)
End do   ! Do j = 2,jedge

!-----Zero out turbulent viscosity at wall
tv(js) = 0.0d0

end

