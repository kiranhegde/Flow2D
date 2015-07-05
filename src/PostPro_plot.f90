!================
program  flow_2d_post
!================
use flux
use vis
use st_end
use pri
use inf
use omp_lib
implicit none

!2-d finite volume solver, cell centred approach

!===============================================================
!n    = total no. of cells in xi direction
!m    = total no. of cells in eta direction
!ndim = space dimension (2 for this code)
!===============================================================
integer::core
integer(kind=i2):: i,j,k,bk,jterm,niter,iter,irks,ii,io!,n,m
integer(kind=i2):: isize(ndim),saveinterval
integer(kind=i2):: nx,ny,nz,istart
real(kind=dp):: time_begin,time_end
real(kind=dp):: xx,yy,x0,y0,delta,q_tmp(ndim+two),r1,r2
real(kind=dp):: cfl_max,cl,cd,xc,cp
real(kind=dp):: x1(ndim),x2(ndim),x3(ndim),x4(ndim),a1
real(kind=dp):: cfl,rk_fac(3)
real(kind=dp):: flux_diff(ndim+two),resi,qtemp(ndim+two)

logical :: file_exists



pi = 4*datan(1.d0)
gamma     = 1.4
gas_const = 287.06

!.......Ficticious cells based on Boundary conditions
bcc(0)=2 ; bcc(1)=1
bcc(2)=1 ; bcc(3)=1
bcc(4)=1 ; bcc(5)=1
bcc(6)=1 ; bcc(7)=1
bcc(8)=0

!........read simulation and  solver input
call read_input(niter,cfl_max,istart)

!.......read grid data
call read_mesh

call normal_zhi
call normal_eta

if ( turb > 0) then
do bk=1,nblk
   call start_end(bk)
   allocate(blmx(cx+1,cy+1,2*ndim))
   do j=1,cy
   do i=1,cx
      blmx(i,j,:)%yy = 0.0d0
      blmx(i,j,:)%rh = 0.0d0
      blmx(i,j,:)%uu = 0.0d0
      blmx(i,j,:)%vo = 0.0d0
      blmx(i,j,:)%vi = 0.0d0
      blmx(i,j,:)%tauw = 0.0d0
      blmx(i,j,:)%tv = 0.0d0
      blmx(i,j,1)%bc  = blk(bk)%mbc(1)
      blmx(i,j,2)%bc  = blk(bk)%mbc(2)
      blmx(i,j,3)%bc  = blk(bk)%mbc(3)
      blmx(i,j,4)%bc  = blk(bk)%mbc(4)
    enddo
    enddo

   call walldistance(bk)
enddo
endif


!.......intialise solution
call intialize


call vertex_data
!call tecplot


k=1

  open(unit=12,file='p_tar')
     bk=1
     j=2
  do i = 1,cx 
    qtemp(:)=blk(bk)%fv(i,j,k)%qvar(:)
    call q_conv(qtemp)
    cp=(p-pinf)/(0.5*qinf(2)*(qinf(3)*qinf(3)+qinf(4)*qinf(4)))
    xc=blk(bk)%mesh(i,j,k)%x
    write(12,*)i,xc,-cp
  enddo
  close(12)

  !101 format(1x,i6,2x,2(f10.6,1x))

!call plot3d
!====================
end program  flow_2d_post
!====================


!====================================================
subroutine read_input(niter,cfl_max,istart)
!====================================================
use flux
use pri
use inf
use vis
implicit none

integer(kind=i2):: niter,istart
real(kind=dp)::  cfl_max

pi = 4.0d0*datan(1.d0)

open(unit=12,file='solver.inp')
!  mach no., freestream pressure, freestream temparature, aoa, exit pressure
read(12,*) 
read(12,*) m_inf,p_inf,t_inf,aoa,pout
!  euler(1)/ns(0) simulation
read(12,*) 
read(12,*) ieuler,turb
!  adiabatic(0)/isothermal(1) wall, wall temperature  
read(12,*) 
read(12,*) iwall,t_w
!  reynolds number     prandtl number
read(12,*) 
read(12,*) re_l,pr_l,pr_t
!  no of iteration, max cfl,  cfl_fraction(dummy here)  beginning(0)/restart(1)
read(12,*) 
read(12,*) niter,cfl_max,cfl_fact,istart
close(12)

aoa = aoa*pi/180.d0
t_w = t_w/t_inf


!========================
end subroutine read_input
!========================


!==========================
subroutine read_mesh
!==========================
use flux
implicit none

integer(kind=i2):: bk,i,j,k,ok,nx,ny,nz,cx,cy,cz
integer(kind=i2):: b1,b2,b3,b4,f1,f2,f3,f4
integer(kind=i2):: bc1,bc2,bc3,bc4
integer(kind=i2):: is,ie,js,je
integer(kind=i2):: m1,n1

ok =100

!grid coordinates include fictitous cells
open(unit=12,file='grid.dat')
read(12,*) nblk
allocate(blk(nblk))

do bk=1,nblk
   read(12,*) nx,ny,nz
   blk(bk)%imax= nx ; blk(bk)%jmax = ny ; blk(bk)%kmax = nz
   blk(bk)%cx= nx-one ; blk(bk)%cy   = ny-one; blk(bk)%cz = nz!-one   
   allocate(blk(bk)%mesh(nx,ny,nz+1),stat=ok)
   if(ok/=0) stop 'allocation failure : read_mesh '
enddo

do bk=1,nblk
   nx=blk(bk)%imax ; ny=blk(bk)%jmax ; nz=blk(bk)%kmax
   read(12,*)(((blk(bk)%mesh(i,j,k)%x,i=1,nx),j=1,ny),k=1,nz), &
             (((blk(bk)%mesh(i,j,k)%y,i=1,nx),j=1,ny),k=1,nz), &
             (((blk(bk)%mesh(i,j,k)%z,i=1,nx),j=1,ny),k=1,nz)
   do j=1,ny
   do i=1,nx
      blk(bk)%mesh(i,j,2)%z=0.001
      blk(bk)%mesh(i,j,2)%x=blk(bk)%mesh(i,j,1)%x
      blk(bk)%mesh(i,j,2)%y=blk(bk)%mesh(i,j,1)%y
   enddo
   enddo
enddo

close(12)
! Reading block connectivity &
! boundary condition at each block boundary
open(3,file='block.inf')
read(3,*)
read(3,*)
read(3,*)
read(3,*)
read(3,*)
do i=1,nblk
read(3,*)bk,b1,b2,b3,b4,f1,f2,f3,f4,bc1,bc2,bc3,bc4
blk(bk)%bc(1)=b1 ; blk(bk)%fc(1)=f1 ; blk(bk)%mbc(1)=bc1
blk(bk)%bc(2)=b2 ; blk(bk)%fc(2)=f2 ; blk(bk)%mbc(2)=bc2
blk(bk)%bc(3)=b3 ; blk(bk)%fc(3)=f3 ; blk(bk)%mbc(3)=bc3
blk(bk)%bc(4)=b4 ; blk(bk)%fc(4)=f4 ; blk(bk)%mbc(4)=bc4

nx=blk(bk)%imax ; ny=blk(bk)%jmax ; nz=blk(bk)%kmax
cx=blk(bk)%cx   ; cy=blk(bk)%cy   ; cz=blk(bk)%cz  

! Start and end of computational cell values
! along i,j & k direction ( based on ghost cells
! at the respective boundary )
bcc(0)=2 ; bcc(1)=1
bcc(2)=1 ; bcc(3)=1
bcc(4)=1 ; bcc(5)=1
bcc(6)=1 ; bcc(7)=1
bcc(8)=0

! i-direction

m1=1
! BC=1
blk(bk)%is=bcc(bc1)+m1

! BC=2
blk(bk)%ie=cx-bcc(bc2)

! j-direction
! BC=3
n1=1
blk(bk)%js=bcc(bc3)+n1
! BC=4
blk(bk)%je=cy-bcc(bc4)

! k-direction
if(ndim==2) then
blk(bk)%ks=1
blk(bk)%ke=nz!-one
!if(b5==0) blk(bk)%ks=gc_blk+1
!if(b6==0) blk(bk)%ke=nz-gc_blk-one
elseif(ndim==3) then
print*
print*,'Wrong   dimension : ndim > 2  not allowed'
print*
stop
endif
nx=blk(bk)%imax ; ny=blk(bk)%jmax ; nz=blk(bk)%kmax
is = blk(bk)%is
ie = blk(bk)%ie
js = blk(bk)%js
je = blk(bk)%je

print *,"==>",bk
print '(" I ","1 -->",i3,"-->",i3,"-->",i3)', is,ie+1,nx
print '(" J ","1 -->",i3,"-->",i3,"-->",i3)', js,je+1,ny

enddo
close(3)

!.......Allocating variables  
do bk=1,nblk

   n1=blk(bk)%cx
   m1=blk(bk)%cy

   nx=blk(bk)%imax
   ny=blk(bk)%jmax
   nz=blk(bk)%kmax

   allocate(blk(bk)%fv(n1,m1,nz),stat=ok)
   if(ok/=0) stop 'allocation failure : fv'
   allocate(blk(bk)%zhi(nx,m1,nz),stat=ok)
   if(ok/=0) stop 'allocation failure : zhi'
   allocate(blk(bk)%eta(n1,ny,nz),stat=ok)
   if(ok/=0) stop 'allocation failure : eta'

enddo

!allocate(delta_w(3,m,ndim),mu(n,m),visturb(n,m))
allocate(zhi_wall(n1),eta_wall(m1))
!bk=1
!nx=blk(bk)%imax
!ny=blk(bk)%jmax
!allocate(x(nx,ny),y(nx,ny))
!n = nx-one !no. of cells
!m = ny-one

!=======================
end subroutine read_mesh
!=======================

!================================================
subroutine intialize
!================================================
use flux
use inf
use st_end
use vis
implicit none

integer(kind=i2):: i,j,k,bk,istart,jterm
logical :: file_exists
k=1
!define intial condition vector
qinf(2) = 1.d0
qinf(3) = dcos(aoa)
qinf(4) = dsin(aoa)
pinf    = 1.d0/(gamma*m_inf*m_inf)
qinf(1) = pinf/(gamma-1.d0) + 0.5d0
pout    = (pout/p_inf)*pinf
rho_inf = qinf(2)
ent_inf = pinf/(rho_inf**gamma)
a_inf   = dsqrt(gamma*pinf/rho_inf)
pdyna   = 0.5d0*rho_inf*(qinf(3)*qinf(3)+qinf(4)*qinf(4))
mu_inf  = rho_inf/re_l

!print*, 'p_dyna =',pdyna

print*
if(ieuler == 0 .and. turb == 0 )  print*,'navier-stokes computation for laminar flow'
if(ieuler == 0 .and. turb > 0 )  then 
   print*,'navier-stokes computation for turbulent flow'
   print*,  ' baldwin-lomax model ' 
endif
if(ieuler == 1 ) print*,'euler computation'
!if(iflow.eq.turbulent)print*,'turbulent navier-stokes computation'
print*
print*,'free-stream values:'
print*  
write(*,'(5x, " mach number =", f8.4)')m_inf
write(*,'(5x, " aoa         =",f8.4)')aoa*180.0/PI
write(*,'(5x, " u velocity  =",f8.4)')qinf(3)
write(*,'(5x, " v velocity  =",f8.4)')qinf(4)
write(*,'(5x, " pressure=", f8.4)')  pinf
print*

inquire(file="q_i.dat", exist=file_exists) ! file_exists will be true if the file
if(.not.file_exists) then
     print*
     print*,'Copy q_o.dat to q_i.dat for  postprocessing...' 
     !call system('cp -fr q_i.dat  q_i.dat.old')
     !call system('cp -fr q_o.dat  q_i.dat')
 stop
endif   

open(unit=12,file='q_i.dat',form='unformatted')
read(12) initer
do bk = 1,nblk
   call start_end(bk)
   read(12) ((( blk(bk)%fv(i,j,k)%qvar(jterm),i=1,cx),j=1,cy),jterm=1,ndim+two)
enddo
close(12)


do bk = 1,nblk
   call start_end(bk)
   if(ieuler == 0 .and. turb == 0 )  then 
          call transpt(bk,cx,cy)
          call walldatas(bk)
   endif
   if(ieuler == 0 .and. turb >0  )  then 
         print*,'Turb'  
         call transpt(bk,cx,cy)
         call walldatas(bk)
         call turbulencemodel(bk)
   endif
enddo



print*

!=======================
end subroutine intialize
!=======================

!========================
subroutine tecplot
!========================
use flux
use st_end
use pri
use vis
implicit none

integer(kind=i2):: i,j,k,bk,ib_m,jb_m,si,ei,sj,ej
integer(kind=i2):: bc1,bc2,bc3,bc4,ii,mm
real(kind=dp):: q_tmp(2+ndim),pm,pr
real(kind=dp):: p_p,p_m,x1,y1,mu,mut,visc
character(len=100)::var,var1
CHARACTER(LEN=30) :: FMT1 
CHARACTER(LEN=3) :: str4 
allocatable  :: p_p(:,:,:) 

k=1
!call blockinterface
var="variables = x,y,rho,u,v,p,t,mach"
mm=8
if(ieuler == 0 .and. turb == 0 )then 
   var1=trim(var)//",mu"
   var=trim(var1)
   mm=mm+one
elseif(ieuler == 0 .and. turb > 0 ) then
   var1=trim(var)//",mu,mut"
   var=trim(var1)
   mm=mm+two
endif
print*,'mm=',mm
print*,var
!var=var(1:len(var))//"'"

write(str4,'(I2)') mm 
str4=adjustr(str4)
FMT1=trim(str4)//"(f15.6,1x)"
fmt1=trim(fmt1)

open(unit=23,file='Tec.dat') 
write(23,*) 'title = "2D CFD"'
write(23,*) var(1:len(var))

do bk = 1,nblk

   call start_end(bk)
   print*,bk,'cx=',cx,'cy=',cy
  
   allocate(p_p(cx,cy,mm)) 

   p_p(:,:,:)=0.d0
   
   do j=1,cy+one   
   do i=1,cx+one   
      p_p(i,j,1)=blk(bk)%mesh(i,j,k)%x
      p_p(i,j,2)=blk(bk)%mesh(i,j,k)%y
   enddo 
   enddo 

   do j=js,je      
   do i=is,ie      

      q_tmp(1:ndim+two) = blk(bk)%fv(i,j,k)%qvar(1:ndim+two)
      call q_conv(q_tmp)
 
      !pr   = p/znd*p_inf
      pr   = 2*(p-znd)
      pm  = dsqrt(u*u+v*v)/a
      p_p(i,j,3)=rho
      p_p(i,j,4)=u
      p_p(i,j,5)=v
      p_p(i,j,6)=pr
      p_p(i,j,7)=t
      p_p(i,j,8)=pm
      if(ieuler == 0 .and. turb == 0 )  p_p(i,j,9)=blk(bk)%fv(i,j,k)%mu
      if(ieuler == 0 .and. turb >  0 )  p_p(i,j,10)=blk(bk)%fv(i,j,k)%mut
   enddo 
   enddo 


   ib_m = ie-is+two
   jb_m = je-js+two

   write(23,*) 'zone  t="',bk,'"  i=',ib_m,'j=',jb_m !,'DATAPACKING=POINT'
   do j =js,je+one 
   do i =is,ie+one 
        write(23,"(1x"//FMT1//")")  (p_p(i,j,ii),ii=1,mm) 
  enddo
  enddo
  deallocate(p_p)
enddo

close(23)



!220 format(1x,i3,3x,3(2x,i4,2x,f15.6))
!201 format(1x,2(f15.6,1x),6(f15.6,1x))
!202 format(1x,2(f15.6,1x),7(f15.6,1x))
!203 format(1x,<mm>(f15.6,1x))
!=====================
end subroutine tecplot
!=====================


!subroutine start_end(bk)
!use st_end
!use flux
!implicit none
!
!integer (kind=i2):: bk,nx,ny,nz
!
!!ci1=0
!!cj1=0
!!ci2=0
!!cj2=0
!! if bc is wall then
!!if(bcc(blk(bk)%mbc(1))==2) ci1 = 1
!!if(bcc(blk(bk)%mbc(2))==2) ci2 = 1
!!if(bcc(blk(bk)%mbc(3))==2) cj1 = 1
!!if(bcc(blk(bk)%mbc(4))==2) cj2 = 1
!is=blk(bk)%is
!ie=blk(bk)%ie
!js=blk(bk)%js
!je=blk(bk)%je
!ks=blk(bk)%ks
!ke=blk(bk)%ke
!
!!nx=blk(bk)%imax  
!!ny=blk(bk)%jmax  
!!nz=blk(bk)%kmax  
!cx=blk(bk)%cx
!cy=blk(bk)%cy
!cz=blk(bk)%cz
!
!end subroutine start_end

!============================================
subroutine plot3d
use flux
use st_end
use pri
use vis
use inf
implicit none

integer(kind=i2):: i,j,k,bk,nx,ny,nz,l,ii,jj,kk
integer(kind=i2):: var,si1,si2,sj1,sj2
real(kind=dp):: q_tmp(2+ndim),visc
real(kind=dp):: p_p,p_m,mu,mut,sol
allocatable :: sol(:,:,:,:)


!open(23,file='FLOW.xyz',form='unformatted')
open(23,file='FLOW.xyz')
!rewind 23
write(23,*)nblk
do bk = 1,nblk
   call start_end(bk)

   si1=is!-one
   si2=ie+one
   sj1=js!-one
   sj2=je+one

   nx = si2-si1+one
   ny = sj2-sj1+one

   write(23,*)nx,ny,2!ke-ks+one
enddo

do bk = 1,nblk
   call start_end(bk)

   si1=is!-one
   si2=ie+one
   sj1=js!-one
   sj2=je+one

   write(23,*)(((blk(bk)%mesh(i,j,k)%x,i=si1,si2),j=sj1,sj2),k=1,2), &
              (((blk(bk)%mesh(i,j,k)%y,i=si1,si2),j=sj1,sj2),k=1,2), &
              (((blk(bk)%mesh(i,j,k)%z,i=si1,si2),j=sj1,sj2),k=1,2)
enddo
close(23)

!Write out 3D variable name Plot3D file:
open(23,file='FLOW.nam')
write(23,*)'u-velocity ; velocity'
write(23,*)'v-velocity'
write(23,*)'w-velocity'
write(23,*)'Density'
write(23,*)'Temparature'
write(23,*)'Presssure'
write(23,*)'Mach No.'
if(ieuler == 0 .and. turb == 0 .or.turb>0)write(23,*)'Mu'
if(ieuler == 0 .and. turb == 1 )write(23,*)'Mut'
close(23)

if(ieuler == 1 ) var = 7
if(ieuler == 0 .and. turb == 0 ) var = 8
if(ieuler == 0 .and. turb > 1 ) var = 9
!var=var+one
!Write out 3D solution Plot3D file:
!open(23,file='FLOW.f',form='unformatted')
open(23,file='FLOW.f')
write(23,*)nblk
do bk=1,nblk
   call start_end(bk)

   si1=is!-one
   si2=ie+one
   sj1=js!-one
   sj2=je+one

   nx = si2-si1+one
   ny = sj2-sj1+one
   nz=ke-ks+one

   if(ieuler == 1 ) write(23,*) nx,ny,2,var
   if(ieuler == 0 .and. turb == 0 ) write(23,*) nx,ny,2 ,var
   if(ieuler == 0 .and. turb > 0 ) write(23,*) nx,ny,2 ,var
enddo

do bk = 1,nblk
   call start_end(bk)

   si1=is!-one
   si2=ie+one
   sj1=js!-one
   sj2=je+one

   nx = si2-si1+one
   ny = sj2-sj1+one
   nz=ke-ks+one

   if(ieuler == 1 )allocate(sol(nx,ny, 2,var))
   if(ieuler == 0 .and. turb == 0 )allocate(sol(nx,ny,2 ,var))
   if(ieuler == 0 .and. turb  > 0 )allocate(sol(nx,ny,2 ,var))

   sol=0.0d0
    do k=1,2
  ! k=1 
   do j=1,ny
   do i=1,nx
      ii=i+is-one
      jj=j+js-one
      kk=k!+ks-one
      !print '(6(i3,1x))',i,j,k,ii,jj,kk  
      q_tmp(:) = blk(bk)%fv(ii,jj,1 )%qvar(:)
      call q_conv(q_tmp)
      !call viscous(t,vis)
      !p_p(i,j) = p/znd*p_inf
      p_p   = 2*(p-znd)
      p_m   = dsqrt(u*u+v*v)/a
      mu = blk(bk)%fv(ii,jj, 1)%mu
      mut = blk(bk)%fv(ii,jj, 1)%mut
      sol(i,j,k,1)=u
      sol(i,j,k,2)=v
      sol(i,j,k,3)=0.0d0
      sol(i,j,k,4)=rho
      sol(i,j,k,5)=t
      sol(i,j,k,6)=p_p
      sol(i,j,k,7)=p_m
      if(ieuler == 0 .and. turb == 0 .or. turb == 1)sol(i,j,k,8)=mu
      if(ieuler == 0 .and. turb >  0 )sol(i,j,k,9)=mut
   enddo
   enddo
   enddo

   do j=1,ny
   do i=1,nx
        sol(i,j,2,:)=sol(i,j,1,:)
   enddo
   enddo

   if(ieuler == 1) write(23,*)((((sol(i,j,k,l),i=1,nx),j=1,ny),k=1, 2),l=1,var)
   if(ieuler == 0 .and. turb == 0 )write(23,*)((((sol(i,j,k,l),i=1,nx),j=1,ny),k=1, 2),l=1,var)
   if(ieuler == 0 .and. turb == 1 )write(23,*)((((sol(i,j,k,l),i=1,nx),j=1,ny),k=1, 2),l=1,var)
   deallocate(sol)
enddo

close(23)
end subroutine plot3d

subroutine vertex_data
use flux
use st_end
use pri
use vis
implicit none

integer(kind=i2):: i,j,k,bk,ib_m,jb_m,si,ei,sj,ej
integer(kind=i2):: bc1,bc2,bc3,bc4,ii,mm,sub
real(kind=dp):: q_tmp(2+ndim),pm,pr,xi,yi,p_p
real(kind=dp):: p_c,p_v,x1,y1,mu,mut,visc,vari
character(len=100)::var,var1
allocatable  :: p_p(:,:,:),vari(:) 
CHARACTER(LEN=3) :: str4
CHARACTER(LEN=30) :: fmt1 

k=1

var="variables = x,y,rho,u,v,p,t,mach"
mm=8
sub=2
if(ieuler == 0 .and. turb == 0 )then 
   var1=trim(var)//",mu"
   var=trim(var1)
   sub=one
   mm=mm+sub
elseif(ieuler == 0 .and. turb == 1 ) then
   var1=trim(var)//",mu,mut"
   var=trim(var1)
   sub=two
   mm=mm+sub
endif
print*,'mm   =',mm
print*,var

write(str4,'(I2)') mm 
str4=adjustr(str4)
FMT1=trim(str4)//"(f15.8,1x)"
fmt1=trim(fmt1)

open(unit=23,file='Tec_vrt.dat') 
write(23,*) 'title = "2D CFD"'
write(23,*) var(1:len(var))

sub=sub-2

do bk = 1,nblk

   call start_end(bk)
   print*,bk,'cx=',cx,'cy=',cy
   
   allocate(p_p(cx+one,cy+one,mm),vari(mm-sub))
   print*, ":",mm,sub,mm-sub
   p_p(:,:,:)=0.d0

   do j=1,cy+one   
   do i=1,cx+one   
      p_p(i,j,1)=blk(bk)%mesh(i,j,k)%x
      p_p(i,j,2)=blk(bk)%mesh(i,j,k)%y
   enddo 
   enddo 


   do j=js,je+one
   do i=is,ie+one
      !call interpol(bk,i,j,vari,mm-two)
      call pseudo_Laplacian(bk,i,j,vari,mm-sub)
      print*,'hi',i,j,mm
      p_p(i,j,3:mm)=vari(1:mm-sub)
   enddo
   enddo
  print*,'kkr' 
   ib_m = ie-is+two
   jb_m = je-js+two 
   
   write(23,*) 'zone  t="',bk,'"  i=',ib_m,'j=',jb_m !,'DATAPACKING=POINT'
   do j =js,je+one 
   do i =is,ie+one 
        !write(23,FMT1)(p_p(i,j,ii),ii=1,mm) 
        write(23,"(1x"//FMT1//")")  (p_p(i,j,ii),ii=1,mm) 
        !write(23,'(1x,  9(f15.6,1x))')(p_p(i,j,ii),ii=1,mm) 
   enddo
   enddo
   deallocate(p_p,vari)
enddo

close(23)

! WRITE(FMT,*) N+1

!           WRITE(6,"(I" // ADJUSTL(FMT) // ")") INT1


!220 format(1x,i3,3x,3(2x,i4,2x,f15.6))
!201 format(1x,2(f15.6,1x),6(f15.6,1x))
!202 format(1x,2(f15.6,1x),7(f15.6,1x))
!203 format(1x,<mm>(f15.6,1x))


end subroutine vertex_data


subroutine primitive(bk,i,j,var,mm)
use dims
use inf
use pri
use vis
use flux
implicit none
integer(kind=i2):: i,j,k,bk,mm
real(kind=dp)::var(mm),q_tmp(ndim+two),pr,pm 

k=1
      q_tmp(1:ndim+two) = blk(bk)%fv(i,j,k)%qvar(1:ndim+two)
      call q_conv(q_tmp)
      !pr   = p/znd*p_inf
      pr   = 2*(p-znd)
      pm  = dsqrt(u*u+v*v)/a
      var(1)=rho
      var(2)=u
      var(3)=v
      var(4)=pr
      var(5)=t
      var(6)=pm
      if(ieuler == 0 .and. turb == 0 ) var(7)=blk(bk)%fv(i,j,k)%mu
      if(ieuler == 0 .and. turb >  0 ) then 
          var(7)=blk(bk)%fv(i,j,k)%mu  
          var(8)=blk(bk)%fv(i,j,k)%mut
      endif
end subroutine primitive


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
real(kind=dp):: qleft(ndim+two),vis
k=1

do j = 1,cy
  do i = 1,cx
    qleft(:) = blk(bk)%fv(i,j,k)%qvar(:)
    call q_conv(qleft)
    call viscous(t,vis)
    blk(bk)%fv(i,j,k)%mu = vis
  enddo
enddo


print *, 'mu min',minval(blk(bk)%fv(:,:,:)%mu)
print *, 'mu max',maxval(blk(bk)%fv(:,:,:)%mu)

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
use vis
implicit none

real(kind=dp):: ts,mu
real(kind=dp):: sut,term1,term2,T_infd,num,den

!sut   = 110.4d0/t_inf
!term1 = dsqrt(ts*ts*ts)
!term2 = (1.d0 + sut)/(ts + sut)
!mu   = term1*term2/re_l

!T_infd  = 300.0d0
!Sut  = 110.4d0*T_inf/T_infd

sut=110.4d0

term1=(ts/t_inf)**1.5
term2=(t_inf+sut)/(ts+sut)
mu  = term1*term2!*mu_inf


!======================
end subroutine viscous
!======================


!====================================
subroutine TurbulenceModel(bk)
!====================================
use dims
use flux
use vis
use st_end
!use pri_deriv
implicit none

integer(kind=i2):: mbc(2*ndim),bk,nx,ny,k
k=1
if ( turb == 0) then
  ! visturb(:,:)=0.0d0
  blk(bk)%fv(:,:,:)%mut=0.0d0
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
real(kind=dp):: qtmp(ndim+two),tmp,amr(ndim+two),aml(ndim+two),qleft(ndim+two),qright(ndim+two)
real(kind=dp):: yy(ny),rh(ny),uu(ny),vo(ny),vi(ny),tv(ny),tauw
real(kind=dp):: u_y,v_x 

!pause 'Tu'


k=1
call start_end(bk)
do i = is,ie!-one
   do j = js,je

      aml(1)  = blk(bk)%zhi(i,j,k)%nx
      aml(2)  = blk(bk)%zhi(i,j,k)%ny
      amr(1)  = blk(bk)%zhi(i+one,j,k)%nx
      amr(2)  = blk(bk)%zhi(i+one,j,k)%ny

      qleft(:)   = blk(bk)%fv(i,j,k)%qvar(:)
      qright(:)  = blk(bk)%fv(i+one,j,k)%qvar(:)
      !call deriv_int_xi (bk,i,j,one,qleft,qright,amr,aml) 
      v_x=blk(bk)%zhi(i,j,k)%v_x
      u_y=blk(bk)%zhi(i,j,k)%u_y
      Blmx(i,j,1)%vo = dabs(v_x-u_y)
   enddo      
enddo      

do j = js,je-one
   do i = is,ie
      aml(1)  = blk(bk)%eta(i,j,k)%nx
      aml(2)  = blk(bk)%eta(i,j,k)%ny
      amr(1)  = blk(bk)%eta(i,j+one,k)%nx
      amr(2)  = blk(bk)%eta(i,j+one,k)%ny

      qleft(:)   = blk(bk)%fv(i,j,k)%qvar(:)
      qright(:)  = blk(bk)%fv(i,j+one,k)%qvar(:)
      !call deriv_int_eta(bk,i,j,two,qleft,qright,amr,aml)
      v_x=blk(bk)%zhi(i,j,k)%v_x
      u_y=blk(bk)%zhi(i,j,k)%u_y
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
      tauw  = zhi_wall(i)%tauw(1)
      call BaldwinLomax(ny,yy,rh,uu,vo,vi,tv,tauw)
      !print*,"i,,     vo            tv"
      !print*,i,vo(jj),tv(jj+4),tauw 
      !pause
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

   dx1 = blk(bk)%eta(i,j+one,k)%nx
   dy1 = blk(bk)%eta(i,j+one,k)%ny
   d2 = dsqrt(dx1*dx1+dy1*dy1)

   dy = 2.0d0*blk(bk)%fv(i,j,k)%vol/(d1+d2)

   Blmx(i,j,1)%yy=Blmx(i,j,1)%yy+0.5d0*dy  
   Blmx(i+one,j,1)%yy=Blmx(i,j,1)%yy+0.5d0*dy  
enddo
enddo
!enddo
endif



! face 2
if(Blmx(1,1,2)%bc==2) then
!do bk=1,nblk
!call start_end(bk)
do i=cx,1,-one
do j=cy,1,-one
   dx1 = blk(bk)%eta(i,j,k)%nx
   dy1 = blk(bk)%eta(i,j,k)%ny
   d1 = dsqrt(dx1*dx1+dy1*dy1)

   dx1 = blk(bk)%eta(i,j+one,k)%nx
   dy1 = blk(bk)%eta(i,j+one,k)%ny
   d2 = dsqrt(dx1*dx1+dy1*dy1)
   dy = 2.0d0*blk(bk)%fv(i,j,k)%vol/(d1+d2)

   Blmx(i-one,j,2)%yy=Blmx(i,j,2)%yy+0.5d0*dy  
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

   dx1 = blk(bk)%zhi(i+one,j,k)%nx
   dy1 = blk(bk)%zhi(i+one,j,k)%ny
   d2 = dsqrt(dx1*dx1+dy1*dy1)
   dy = 2.0d0*blk(bk)%fv(i,j,k)%vol/(d1+d2)
   Blmx(i,j,3)%yy=Blmx(i,j,3)%yy+0.5d0*dy  
   Blmx(i,j+one,3)%yy=Blmx(i,j,3)%yy+0.5d0*dy  
enddo
enddo
!enddo
endif

! face 4
if(Blmx(1,1,4)%bc==2) then
!do bk=1,nblk
!call start_end(bk)
do j=cy,1,-one
do i=cx,1,-one
   dx1 = blk(bk)%zhi(i,j,k)%nx
   dy1 = blk(bk)%zhi(i,j,k)%ny
   d1 = dsqrt(dx1*dx1+dy1*dy1)

   dx1 = blk(bk)%zhi(i+one,j,k)%nx
   dy1 = blk(bk)%zhi(i+one,j,k)%ny
   d2 = dsqrt(dx1*dx1+dy1*dy1)
   dy = 2.0d0*blk(bk)%fv(i,j,k)%vol/(d1+d2)

   Blmx(i,j-one,4)%yy=Blmx(i,j,4)%yy+0.5d0*dy  
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
!do j=js,js+one
!   x1=x(i,j) ; y1 = y(i,j) 
!   x2=x(i+one,j) ; y2 = y(i+one,j) 
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
subroutine BaldwinLomax(m,yy,rh,uu,vo,vi,tv,tauw) 
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
!Real(kind=dp) :: tauw        ! Shear stress at the wall

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
Real(kind=dp) tauw 

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


!=================================================
subroutine normal_zhi
!=================================================

!calculates normals on zhi faces for all the computational cells

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!cx    = same as n for this code
!cy    = same as n for this code
!x1    = x coordinates
!y1    = y coordinates
!n_zhi = normals on the xi faces
!===============================================================
use flux
use st_end
implicit none

integer(kind=i2):: i,j,k,bk
k=1
do bk=1,nblk
   call start_end(bk)
i = 1
do j = 1,cy
do i = 1,is-one
  blk(bk)%zhi(i,j,k)%nx = -(blk(bk)%mesh(i,j+one,k)%y - blk(bk)%mesh(i,j,k)%y)
  blk(bk)%zhi(i,j,k)%ny =  (blk(bk)%mesh(i,j+one,k)%x - blk(bk)%mesh(i,j,k)%x)
enddo
enddo

do j = 1,cy
  do i = is,cx+one
    blk(bk)%zhi(i,j,k)%nx =  (blk(bk)%mesh(i,j+one,k)%y - blk(bk)%mesh(i,j,k)%y)
    blk(bk)%zhi(i,j,k)%ny = -(blk(bk)%mesh(i,j+one,k)%x - blk(bk)%mesh(i,j,k)%x)
  enddo
enddo
enddo

!=========================
end subroutine normal_zhi
!=========================




!=================================================
subroutine normal_eta
!=================================================

!calculates normals on eta faces for all the computational cells

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!cx    = same as n for this code
!cy    = same as n for this code
!x1    = x coordinates
!y1    = y coordinates
!n_eta = normals on the eta faces
!===============================================================
use st_end
use flux
implicit none

integer(kind=i2):: i,j,k,bk

k=1
do bk=1,nblk
   call start_end(bk)
j = 1
do j = 1,js-one
do i = 1,cx
    blk(bk)%eta(i,j,k)%nx =  (blk(bk)%mesh(i+one,j,k)%y - blk(bk)%mesh(i,j,k)%y)
    blk(bk)%eta(i,j,k)%ny = -(blk(bk)%mesh(i+one,j,k)%x - blk(bk)%mesh(i,j,k)%x)
enddo
enddo

do j = js,cy+one
  do i = 1,cx
    blk(bk)%eta(i,j,k)%nx = -(blk(bk)%mesh(i+one,j,k)%y - blk(bk)%mesh(i,j,k)%y)
    blk(bk)%eta(i,j,k)%ny =  (blk(bk)%mesh(i+one,j,k)%x - blk(bk)%mesh(i,j,k)%x)
  enddo
enddo
enddo


!=========================
end subroutine normal_eta
!=========================

subroutine  interpol(bk,i,j,vari,mm)
use dims
use flux
use st_end
implicit none
integer(kind=i2):: i,j,k,bk,mm
real(kind=dp)::xi,yi,vari(mm),var(mm),r1,r2,r3,r4,rsum 
real(kind=dp)::cell_centroid_dist,x1,y1

k=1
      vari(:)=0.d0
      xi=blk(bk)%mesh(i,j,k)%x
      yi=blk(bk)%mesh(i,j,k)%y

      r1=cell_centroid_dist(bk,i-one,j-one,xi,yi)
      if(i==2.or.j==2) r1=0.0
      call primitive(bk,i-one,j-one,var,mm)
      vari(:)=var(:)*r1

      r2=cell_centroid_dist(bk,i-one,j  ,xi,yi)
      if(i==2.or.j==cy-one) r2=0.0
      call primitive(bk,i-one,j  ,var,mm)
      vari(:)=vari(:)+var(:)*r2

      r3=cell_centroid_dist(bk,i  ,j  ,xi,yi)
      if(i==cx-one.or.j==cy-one) r3=0.0
      call primitive(bk,i  ,j  ,var,mm)
      vari(:)=vari(:)+var(:)*r3

      r4=cell_centroid_dist(bk,i  ,j-one,xi,yi)
      if(i==cx-one.or.j==2) r4=0.0
      call primitive(bk,i  ,j-one,var,mm)
      vari(:)=vari(:)+var(:)*r4

      rsum=r1+r2+r3+r4
 
      vari(:)=vari(:)/rsum

end subroutine  interpol

real(kind=dp)  function cell_centroid_dist(bk,i,j,xi,yi)
use dims
use flux
implicit none

integer(kind=i2):: i,j,bk
real(kind=dp)::xi,yi,x1,y1,dist

     call cell_centroid(bk,i,j,x1,y1)
     dist=dsqrt( (xi-x1)**2+(yi-y1)**2)
     dist=1.0/dist
     cell_centroid_dist=dist

end function  cell_centroid_dist

subroutine  cell_centroid(bk,i,j,xi,yi)
use dims
use flux
implicit none

integer(kind=i2):: i,j,k,bk,ii
real(kind=dp)::xi,yi 
real(kind=dp)::x(5),y(5),area,tmp

k=1
     x(1)=blk(bk)%mesh(i,j,k)%x
     x(2)=blk(bk)%mesh(i+one,j,k)%x
     x(3)=blk(bk)%mesh(i+one,j+one,k)%x
     x(4)=blk(bk)%mesh(i,j+one,k)%x
     x(5)=x(1)
     y(1)=blk(bk)%mesh(i,j,k)%y
     y(2)=blk(bk)%mesh(i+one,j,k)%y
     y(3)=blk(bk)%mesh(i+one,j+one,k)%y
     y(4)=blk(bk)%mesh(i,j+one,k)%y
     y(5)=y(1)

     area=0.d0
     xi=0.d0
     yi=0.d0
     do ii=1,4
        tmp=x(ii)*y(ii+one)-x(ii+one)*y(ii)
        area=area+tmp
        xi=xi+tmp*(x(ii)+x(ii+one))
        yi=yi+tmp*(y(ii)+y(ii+one))
     enddo
     area=0.5d0*area
     xi=xi/6.d0/area
     yi=yi/6.d0/area
     !x1=0.25*(x(1)+x(2)+x(3)+ x(4))
     !y1=0.25*(y(1)+y(2)+y(3)+ y(4))

end subroutine cell_centroid


subroutine  pseudo_Laplacian(bk,i,j,vari,mm)
use dims
use st_end
use flux
implicit none

integer(kind=i2):: i,j,k,bk,ii,mm
real(kind=dp)::x(4),y(4),area,tmp,x1,y1,cell_centroid_dist
real(kind=dp)::xi,yi,vari(mm),var(mm),r1,r2,r3,r4,rsum 
real(kind=dp)::ixx,iyy,ixy,rx,ry,gx,gy,dx,dy,det 

k=1
vari(:)=0.d0
xi=blk(bk)%mesh(i,j,k)%x
yi=blk(bk)%mesh(i,j,k)%y
ixx=0.d0 ; iyy=0.d0
ixy=0.d0
rx=0.d0 ; ry=0.d0
gx=0.d0 ; gy=0.d0


   !print*,'hi',i,j
   call cell_centroid(bk,i-one,j-one,x1,y1)
   dx=x1-xi ; dy=y1-yi
   rx=dx    ; ry=dy   
   ixx=ixx+dx*dx
   iyy=iyy+dy*dy
   ixy=ixy+dx*dy
   x(1)=dx  ; y(1)=dy

   call cell_centroid(bk,i-one,j  ,x1,y1)
   dx=x1-xi ; dy=y1-yi
   rx=rx+dx    
   ry=ry+dy   
   ixx=ixx+dx*dx
   iyy=iyy+dy*dy
   ixy=ixy+dx*dy
   x(2)=dx  ; y(2)=dy

   call cell_centroid(bk,i  ,j  ,x1,y1)
   dx=x1-xi ; dy=y1-yi 
   rx=rx+dx    
   ry=ry+dy   
   ixx=ixx+dx*dx
   iyy=iyy+dy*dy
   ixy=ixy+dx*dy
   x(3)=dx  ; y(3)=dy

   call cell_centroid(bk,i  ,j-one,x1,y1)
   dx=x1-xi ; dy=y1-yi
   rx=rx+dx    
   ry=ry+dy   
   ixx=ixx+dx*dx
   iyy=iyy+dy*dy
   ixy=ixy+dx*dy
   x(4)=dx  ; y(4)=dy

   det=ixx*iyy-ixy*ixy 
   if(det==0.d0) then
      print*,'det=',det
   endif
   gx=(ixy*ry-iyy*rx)/det
   gy=(ixy*rx-ixx*ry)/det
 
   call start_end(bk)

   r1=1+gx*x(1)+gy*y(1)
   call primitive(bk,i-one,j-one,var,mm)
   vari(:)=var(:)*r1

   r2=1+gx*x(2)+gy*y(2)
   call primitive(bk,i-one,j  ,var,mm)
   vari(:)=vari(:)+var(:)*r2

   r3=1+gx*x(3)+gy*y(3)
   call primitive(bk,i  ,j  ,var,mm)
   vari(:)=vari(:)+var(:)*r3

   r4=1+gx*x(4)+gy*y(4)
   call primitive(bk,i  ,j-one,var,mm)
   vari(:)=vari(:)+var(:)*r4

   rsum=r1+r2+r3+r4
   vari(:)=vari(:)/rsum

end subroutine  pseudo_Laplacian

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
real(kind=dp):: qleft(ndim+two),qright(ndim+two)
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
eta_x2 = blk(bk)%eta(icell  ,jcell+one,k)%nx
eta_x3 = blk(bk)%eta(icell+one,jcell+one,k)%nx
eta_x4 = blk(bk)%eta(icell+one,jcell  ,k)%nx
eta_y1 = blk(bk)%eta(icell  ,jcell  ,k)%ny
eta_y2 = blk(bk)%eta(icell  ,jcell+one,k)%ny
eta_y3 = blk(bk)%eta(icell+one,jcell+one,k)%ny
eta_y4 = blk(bk)%eta(icell+one,jcell  ,k)%ny

u_zhi   = uright - uleft
v_zhi   = vright - vleft
t_zhi   = tright - tleft

do jterm = 1,ndim+two
  qleft(jterm)  = blk(bk)%fv(icell,jcell-one,k)%qvar(jterm)
  qright(jterm) = blk(bk)%fv(icell+one,jcell-one,k)%qvar(jterm)
enddo

call q_conv(qleft)
ubotl = u
vbotl = v
tbotl = t

call q_conv(qright)
ubotr = u
vbotr = v
tbotr = t

do jterm = 1,ndim+two
  qleft(jterm)  = blk(bk)%fv(icell,jcell+one,k)%qvar(jterm)
  qright(jterm) = blk(bk)%fv(icell+one,jcell+one,k)%qvar(jterm)
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
real(kind=dp):: qleft(ndim+two),qright(ndim+two)
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
zhi_x2 = blk(bk)%zhi(icell  ,jcell+one,k)%nx
zhi_x3 = blk(bk)%zhi(icell+one,jcell+one,k)%nx
zhi_x4 = blk(bk)%zhi(icell+one,jcell  ,k)%nx
zhi_y1 = blk(bk)%zhi(icell  ,jcell  ,k)%ny
zhi_y2 = blk(bk)%zhi(icell  ,jcell+one,k)%ny
zhi_y3 = blk(bk)%zhi(icell+one,jcell+one,k)%ny
zhi_y4 = blk(bk)%zhi(icell+one,jcell  ,k)%ny

u_eta   = uright - uleft
v_eta   = vright - vleft
t_eta   = tright - tleft

do jterm = 1,ndim+two
  qleft(jterm)  = blk(bk)%fv(icell-one,jcell,k)%qvar(jterm)
  qright(jterm) = blk(bk)%fv(icell+one,jcell,k)%qvar(jterm)
enddo

call q_conv(qleft)
ubotl = u
vbotl = v
tbotl = t

call q_conv(qright)
ubotr = u
vbotr = v
tbotr = t

do jterm = 1,ndim+two
  qleft(jterm)  = blk(bk)%fv(icell-one,jcell+one,k)%qvar(jterm)
  qright(jterm) = blk(bk)%fv(icell+one,jcell+one,k)%qvar(jterm)

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


!==================================================
subroutine cent_fac(bk,i,j,idir,xfc,yfc)
!==================================================

!calculates the face centres of a given cell face

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!i     = xi index
!j     = eta index
!idir  = 1 for xi direction, 2 for eta direction
!x     = x coordinates
!y     = y coordinates
!xfc   = x coordinate for the face centre
!yfc   = y coordinate for the face centre
!===============================================================
use flux

implicit none

integer(kind=i2):: bk,i,j,k,idir
real(kind=dp):: xfc,yfc

k=1

select case(idir)

case(1)
  xfc = 0.5d0*(blk(bk)%mesh(i+one,j,k)%x+blk(bk)%mesh(i+one,j+one,k)%x)
  yfc = 0.5d0*(blk(bk)%mesh(i+one,j,k)%y+blk(bk)%mesh(i+one,j+one,k)%y)
case(2)
  xfc = 0.5d0*(blk(bk)%mesh(i,j+one,k)%x+blk(bk)%mesh(i+one,j+one,k)%x)
  yfc = 0.5d0*(blk(bk)%mesh(i,j+one,k)%y+blk(bk)%mesh(i+one,j+one,k)%y)
end select


!=======================
end subroutine cent_fac
!=======================


!===================================================
subroutine centroid1(bk,i,j,idir,xfc,yfc)
!===================================================

!calculates the centroid of a given cell

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!i     = xi index
!j     = eta index
!idir  = 1 for xi direction, 2 for eta direction
!xfc   = x coordinate for the centroid
!yfc   = y coordinate for the centroid
!x     = x coordinates
!y     = y coordinates
!===============================================================
use dims
use flux
implicit none
integer(kind=i2):: i,j,k,bk,idir
real(kind=dp):: xfc,yfc

integer(kind=i2):: ip1,jp1
real(kind=dp):: xa,xb,xc,xd,ya,yb,yc,yd
real(kind=dp):: third
real(kind=dp):: xc1,yc1,termx,termy,dab,dbc,dca,speri,area1
real(kind=dp):: xc2,yc2,dda,dac,dcd,area2,area
k=1
third = 1.d0/3.d0
!
!initializing the vertices for zhi direction
!
ip1 = i+one
jp1 = j+one

xa  = blk(bk)%mesh(i  ,j,k)%x
xb  = blk(bk)%mesh(i  ,jp1,k)%x
xc  = blk(bk)%mesh(ip1,jp1,k)%x
xd  = blk(bk)%mesh(ip1,j,k)%x

ya  = blk(bk)%mesh(i  ,j,k)%y
yb  = blk(bk)%mesh(i  ,jp1,k)%y
yc  = blk(bk)%mesh(ip1,jp1,k)%y
yd  = blk(bk)%mesh(ip1,j,k)%y

!
!dividing into two triangles abc  dac
!
!first triangle
!
xc1   = (xa + xb + xc)*third
yc1   = (ya + yb + yc)*third

termx = xb - xa
termy = yb - ya
dab   = dsqrt(termx*termx +  termy*termy)

termx = xc - xb
termy = yc - yb
dbc   = dsqrt(termx*termx + termy*termy)

termx = xa - xc
termy = ya - yc
dca   = dsqrt(termx*termx + termy*termy)

speri = 0.5d0*(dab + dbc + dca)
area1 = dsqrt(speri*(speri - dab)*(speri - dbc)*(speri - dca))

!
!second triangle
!
xc2   = (xd + xa + xc)*third
yc2   = (yd + ya + yc)*third

termx = xa - xd
termy = ya - yd
dda   = dsqrt(termx*termx + termy*termy)

dac   = dca

termx = xd - xc
termy = yd - yc
dcd   = dsqrt(termx*termx + termy*termy)

speri = 0.5d0*(dda + dac + dcd)
area2 = dsqrt(speri*(speri - dda)*(speri - dac)*(speri - dcd))
!
!calculation of cell face centroids
!
area = area1 + area2

if(area .ne. 0.d0) then
  xfc = (xc1*area1 + xc2*area2)/area
  yfc = (yc1*area1 + yc2*area2)/area
else
  xfc = 0.5d0*(xc1 + xc2)
  yfc = 0.5d0*(yc1 + yc2)
endif

!========================
end subroutine centroid1
!========================




!===========================================================
subroutine cell_dist(bk,cx,cy,icell,jcell,idir,iside,r1,r2,r3,r4)!,area)
!===========================================================

!for a given interface, this subroutine computes the cell center
!distances of two cells on either sides from the given interface

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!cx    = same as n for this code
!cy    = same as n for this code
!icell = xi index
!jcell = eta index
!idir  = 1 for xi direction, 2 for eta direction
!iside = 1 for left interface, 2 for right interface
!r1    = distance of the interface to the cell centre of 2nd cell
!        to the left
!r2    = distance of the interface to the cell centre of 1st cell
!        to the left
!r3    = distance of the interface to the cell centre of 1st cell
!        to the right
!r4    = distance of the interface to the cell centre of 2nd cell
!        to the right
!n_zhi = normals on the xi faces
!n_eta = normals on the eta faces
!x     = x coordinates
!y     = y coordinates
!area  = area of the cell
!===============================================================
use flux
implicit none
integer(kind=i2):: bk,k,icell,jcell,idir,iside,cx,cy
real(kind=dp):: r1,r2,r3,r4

real(kind=dp):: amx,amy,grad_met,xfc,yfc
real(kind=dp):: x1,x2,x3,x4,y1,y2,y3,y4,dx,dy
k=1

select case(idir)

case(1)
  call cent_fac(bk,icell,jcell,idir,xfc,yfc)
  amx = blk(bk)%zhi(icell+one,jcell,k)%nx
  amy = blk(bk)%zhi(icell+one,jcell,k)%ny

  grad_met = dsqrt(amx**2 + amy**2)
  amx = amx/grad_met
  amy = amy/grad_met

  select case(iside)

  case(1)

    call centroid(bk,icell+two,jcell,idir,x1,y1)
    !r1 = dabs((x1-xfc)*amx + (y1-yfc)*amy)
    dx=(x1-xfc)*amx  ; dy=(y1-yfc)*amy
    r1 = dsqrt(dx*dx+dy*dy) 

    call centroid(bk,icell+one,jcell,idir,x2,y2)
    !r2 = dabs((x2-xfc)*amx + (y2-yfc)*amy)
    dx=(x2-xfc)*amx  ; dy=(y2-yfc)*amy
    r2 = dsqrt(dx*dx+dy*dy) 

    call centroid(bk,icell,jcell,idir,x3,y3)
    !r3 = dabs((x3-xfc)*amx + (y3-yfc)*amy)
    dx=(x3-xfc)*amx  ; dy=(y3-yfc)*amy
    r3 = dsqrt(dx*dx+dy*dy) 

    if (icell .ne. 1) then
      call centroid(bk,icell-one,jcell,idir,x4,y4)
      !r4 = dabs((x4-xfc)*amx + (y4-yfc)*amy)
      dx=(x4-xfc)*amx  ; dy=(y4-yfc)*amy
      r4 = dsqrt(dx*dx+dy*dy) 
    else
      !r4 = 3.d0*(r2+r3) - r1
      r4 = 3.d0*(r2+r3) - r1
    endif

  case(2)

    call centroid(bk,icell-one,jcell,idir,x1,y1)
    !r1 = dabs((x1-xfc)*amx + (y1-yfc)*amy)
    dx=(x1-xfc)*amx  ; dy=(y1-yfc)*amy
    r1 = dsqrt(dx*dx+dy*dy) 

    call centroid(bk,icell,jcell,idir,x2,y2)
    !r2 = dabs((x2-xfc)*amx + (y2-yfc)*amy)
    dx=(x2-xfc)*amx  ; dy=(y2-yfc)*amy
    r2 = dsqrt(dx*dx+dy*dy) 

    call centroid(bk,icell+one,jcell,idir,x3,y3)
    !r3 = dabs((x3-xfc)*amx + (y3-yfc)*amy)
    dx=(x3-xfc)*amx  ; dy=(y3-yfc)*amy
    r3 = dsqrt(dx*dx+dy*dy) 

    if (icell .ne. cx-one) then
      call centroid(bk,icell+two,jcell,idir,x4,y4)
      !r4 = dabs((x4-xfc)*amx + (y4-yfc)*amy)
      dx=(x4-xfc)*amx  ; dy=(y4-yfc)*amy
      r4 = dsqrt(dx*dx+dy*dy) 
     else
      r4 = 3.d0*(r2+r3) - r1
    endif

  end select

case(2)

  call cent_fac(bk,icell,jcell,idir,xfc,yfc)
  amx = blk(bk)%eta(icell,jcell+one,k)%nx
  amy = blk(bk)%eta(icell,jcell+one,k)%ny

  grad_met = dsqrt(amx**2 + amy**2)
  amx = amx/grad_met
  amy = amy/grad_met

  select case(iside)

  case(1)

    call centroid(bk,icell,jcell+two,idir,x1,y1)
    !r1 = dabs((x1-xfc)*amx + (y1-yfc)*amy)
    dx=(x1-xfc)*amx  ; dy=(y1-yfc)*amy
    r1 = dsqrt(dx*dx+dy*dy) 

    call centroid(bk,icell,jcell+one,idir,x2,y2)
    !r2 = dabs((x2-xfc)*amx + (y2-yfc)*amy)
    dx=(x2-xfc)*amx  ; dy=(y2-yfc)*amy
    r2 = dsqrt(dx*dx+dy*dy) 

    call centroid(bk,icell,jcell,idir,x3,y3)
    !r3 = dabs((x3-xfc)*amx + (y3-yfc)*amy)
    dx=(x3-xfc)*amx  ; dy=(y3-yfc)*amy
    r3 = dsqrt(dx*dx+dy*dy) 

    if(jcell .ne. 1) then
      call centroid(bk,icell,jcell-one,idir,x4,y4)
      !r4 = dabs((x4-xfc)*amx + (y4-yfc)*amy)
      dx=(x4-xfc)*amx  ; dy=(y4-yfc)*amy
      r4 = dsqrt(dx*dx+dy*dy) 
     else
      r4 = 3.d0*(r2+r3) - r1
    endif

  case(2)

    call centroid(bk,icell,jcell-one,idir,x1,y1)
    !r1 = dabs((x1-xfc)*amx + (y1-yfc)*amy)
    dx=(x1-xfc)*amx  ; dy=(y1-yfc)*amy
    r1 = dsqrt(dx*dx+dy*dy) 

    call centroid(bk,icell,jcell,idir,x2,y2)
    !r2 = dabs((x2-xfc)*amx + (y2-yfc)*amy)
    dx=(x2-xfc)*amx  ; dy=(y2-yfc)*amy
    r2 = dsqrt(dx*dx+dy*dy) 

    call centroid(bk,icell,jcell+one,idir,x3,y3)
    !r3 = dabs((x3-xfc)*amx + (y3-yfc)*amy)
    dx=(x3-xfc)*amx  ; dy=(y3-yfc)*amy
    r3 = dsqrt(dx*dx+dy*dy) 

    if(jcell .ne. cy-one) then
      call centroid(bk,icell,jcell+two,idir,x4,y4)
      !r4 = dabs((x4-xfc)*amx + (y4-yfc)*amy)
      dx=(x4-xfc)*amx  ; dy=(y4-yfc)*amy
      r4 = dsqrt(dx*dx+dy*dy) 
     else
      r4 = 3.d0*(r2+r3) - r1
    endif

  end select

end select

!========================
end subroutine cell_dist
!========================




!==================================================
subroutine centroid(bk,i,j,idir,xfc,yfc)
!==================================================

!simple centroid calculation

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!i     = xi index
!j     = eta index
!idir  = 1 for xi direction, 2 for eta direction
!xfc   = x coordinate for the centroid
!yfc   = y coordinate for the centroid
!x     = x coordinates
!y     = y coordinates
!===============================================================
use flux

implicit none
integer(kind=i2):: i,j,k,bk,idir
real(kind=dp):: xfc,yfc

integer(kind=i2):: ip1,jp1
real(kind=dp):: xa,xb,xc,xd,ya,yb,yc,yd
real(kind=dp):: xac,yac,xbd,ybd
k=1
ip1 = i+one
jp1 = j+one

xa  = blk(bk)%mesh(i  ,j,k)%x
xb  = blk(bk)%mesh(i  ,jp1,k)%x
xc  = blk(bk)%mesh(ip1,jp1,k)%x
xd  = blk(bk)%mesh(ip1,j,k)%x

ya  = blk(bk)%mesh(i  ,j,k)%y
yb  = blk(bk)%mesh(i  ,jp1,k)%y
yc  = blk(bk)%mesh(ip1,jp1,k)%y
yd  = blk(bk)%mesh(ip1,j,k)%y

xac = 0.5d0*(xa + xc)
yac = 0.5d0*(ya + yc)

xbd = 0.5d0*(xb + xd)
ybd = 0.5d0*(yb + yd)

xfc = 0.5d0*(xac + xbd)
yfc = 0.5d0*(yac + ybd)

!=======================
end subroutine centroid
!=======================


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
real(kind=dp):: q1(ndim+two),q2(ndim+two),q3(ndim+two)!,fintr(ndim+two)
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
  amr(1)  = blk(bk)%zhi(icell+one,j,k)%nx
  amr(2)  = blk(bk)%zhi(icell+one,j,k)%ny
  amfr(1)  = blk(bk)%zhi(icell+two,j,k)%nx
  amfr(2)  = blk(bk)%zhi(icell+two,j,k)%ny


  q1(:) = blk(bk)%fv(icell+two,j,k)%qvar(:)
  q2(:) = blk(bk)%fv(icell+one,j,k)%qvar(:)
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
  amr(1)  = blk(bk)%zhi(icell+one,j,k)%nx
  amr(2)  = blk(bk)%zhi(icell+one,j,k)%ny
  amfr(1)  = blk(bk)%zhi(icell+two,j,k)%nx
  amfr(2)  = blk(bk)%zhi(icell+two,j,k)%ny


  q1(:) = blk(bk)%fv(icell-one,j,k)%qvar(:)
  q2(:) = blk(bk)%fv(icell,j,k)%qvar(:)
  q3(:) = blk(bk)%fv(icell+one,j,k)%qvar(:)

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
  amr(1)  = blk(bk)%eta(i,jcell+one,k)%nx
  amr(2)  = blk(bk)%eta(i,jcell+one,k)%ny
  amfr(1)  = blk(bk)%eta(i,jcell+two,k)%nx
  amfr(2)  = blk(bk)%eta(i,jcell+two,k)%ny


  q1(:) = blk(bk)%fv(i,jcell+two,k)%qvar(:)
  q2(:) = blk(bk)%fv(i,jcell+one,k)%qvar(:)
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
  amr(1)  = blk(bk)%eta(i,jcell+one,k)%nx
  amr(2)  = blk(bk)%eta(i,jcell+one,k)%ny
  amfr(1)  = blk(bk)%eta(i,jcell+two,k)%nx
  amfr(2)  = blk(bk)%eta(i,jcell+two,k)%ny

  q1(:) = blk(bk)%fv(i,jcell-one,k)%qvar(:)
  q2(:) = blk(bk)%fv(i,jcell,k)%qvar(:)
  q3(:) = blk(bk)%fv(i,jcell+one,k)%qvar(:)

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
real(kind=dp):: q1(ndim+two),q2(ndim+two),q3(ndim+two)
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
!write(506,*) X(ICELL,JCELL+one),sf
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
   wbc(1,1)=one        ; wbc(1,2)=is-one ; wbc(1,3)=one        ; wbc(1,4)=cy
   wbc(2,1)=ie+one     ; wbc(2,2)=cx   ; wbc(2,3)=one        ; wbc(2,4)=cy
   wbc(3,1)=one        ; wbc(3,2)=cx   ; wbc(3,3)=one        ; wbc(3,4)=js-one
   wbc(4,1)=one        ; wbc(4,2)=cx   ; wbc(4,3)=je+one     ; wbc(4,4)=cy  
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
         blk(bk)%fv(i,j,k)%qvar(1:ndim+two) =blk(kk)%fv(ii,jj,k)%qvar(1:ndim+two)
         blk(bk)%fv(i,j,k)%qold(1:ndim+two) =blk(kk)%fv(ii,jj,k)%qold(1:ndim+two)
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
         blk(bk)%fv(i,j,k)%qvar(1:ndim+two) =blk(kk)%fv(ii,jj,k)%qvar(1:ndim+two)
         blk(bk)%fv(i,j,k)%qold(1:ndim+two) =blk(kk)%fv(ii,jj,k)%qold(1:ndim+two)
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
         blk(bk)%fv(i,j,k)%qvar(1:ndim+two) =blk(kk)%fv(ii,jj,k)%qvar(1:ndim+two)
         blk(bk)%fv(i,j,k)%qold(1:ndim+two) =blk(kk)%fv(ii,jj,k)%qold(1:ndim+two)
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
         blk(bk)%fv(i,j,k)%qvar(1:ndim+two) =blk(kk)%fv(ii,jj,k)%qvar(1:ndim+two)
         blk(bk)%fv(i,j,k)%qold(1:ndim+two) =blk(kk)%fv(ii,jj,k)%qold(1:ndim+two)
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
