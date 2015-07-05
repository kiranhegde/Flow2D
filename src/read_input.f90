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

!pi = 4.0d0*datan(1.d0)

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
   blk(bk)%cx= nx-one ; blk(bk)%cy   = ny-one; blk(bk)%cz = nz!-1   
   allocate(blk(bk)%mesh(nx,ny,nz+1),stat=ok)
   if(ok/=0) stop 'allocation failure : read_mesh '
enddo

do bk=1,nblk
   nx=blk(bk)%imax ; ny=blk(bk)%jmax ; nz=blk(bk)%kmax
   read(12,*)(((blk(bk)%mesh(i,j,k)%x,i=1,nx),j=1,ny),k=1,nz), &
             (((blk(bk)%mesh(i,j,k)%y,i=1,nx),j=1,ny),k=1,nz), &
             (((blk(bk)%mesh(i,j,k)%z,i=1,nx),j=1,ny),k=1,nz)
!open(unit=12,file=fname)
!read(12,*) NX,NY ! no. of nodes
!allocate(x(nx,ny),y(nx,ny))
!read(12,*) ((x(i,j),y(i,j),i=1,nx),j=1,ny)
!   read(12,*)(((blk(bk)%mesh(i,j,k)%x, blk(bk)%mesh(i,j,k)%y, i=1,nx),j=1,ny),k=1,nz)
!             blk(bk)%mesh(:,:,:)%z=0.d0
!close(12)

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
!boundary types:
!bcc(0) : block interface
!bcc(1) : characteristic
!bcc(2) : wall
!bcc(3) : extrapolation
!bcc(4) : freestream
!bcc(5) : symmetry
!bcc(6) : periodic
!bcc(7) : outflow pressure specified

bcc(0)=2 
bcc(1)=1
bcc(2)=1 
bcc(3)=1
bcc(4)=1 
bcc(5)=1
bcc(6)=1 
bcc(7)=1
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
blk(bk)%ke=nz!-1
!if(b5==0) blk(bk)%ks=gc_blk+1
!if(b6==0) blk(bk)%ke=nz-gc_blk-1
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
print '(i3,"  I  ","1 -->",i3,"-->",i3,"-->",i3,"  J " , & 
                   "1 -->",i3,"-->",i3,"-->",i3)',bk,is,ie+1,nx,js,je+1,ny
               

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
!allocate(zhi_wall(n1),eta_wall(m1))
!bk=1
!nx=blk(bk)%imax
!ny=blk(bk)%jmax
!allocate(x(nx,ny),y(nx,ny))
!n = nx-1 !no. of cells
!m = ny-1

!=======================
end subroutine read_mesh
!=======================

!================================================
subroutine intialize(istart)
!================================================
use flux
use inf
use st_end
use vis
implicit none

integer(kind=i2):: i,j,k,bk,istart,jterm
real(kind=dp):: aoa1
logical :: file_exists
k=1
aoa1= aoa*pi/180.d0
!define intial condition vector
qinf(2) = 1.d0
qinf(3) = dcos(aoa1)
qinf(4) = dsin(aoa1)
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
if(ieuler == 0 .and. turb == 1 )  then 
   print*,'navier-stokes computation for turbulent flow'
   print*,  ' baldwin-lomax model ' 
endif
if(ieuler == 1 ) print*,'euler computation'
!if(iflow.eq.turbulent)print*,'turbulent navier-stokes computation'
print*
print*,'free-stream values:'
print*  
write(*,'(5x, " mach number =", f8.4)')m_inf
write(*,'(5x, " aoa         =",f8.4)')aoa
write(*,'(5x, " u velocity  =",f8.4)')qinf(3)
write(*,'(5x, " v velocity  =",f8.4)')qinf(4)
write(*,'(5x, " pressure=", f8.4)')  pinf
print*

!.......initialize qvar
if(istart == 0) then

  initer  = 0


  do bk=1,nblk
     call start_end(bk)
     do j = 1,cy
        do i = 1,cx
           blk(bk)%fv(i,j,k)%qvar(:)  = qinf(:)
        enddo
        if(blk(bk)%mbc(2) .eq. 7) then
           blk(bk)%fv(ie+1,j,k)%qold(1)  = pout/(gamma-1.d0) + 0.5d0
           !blk(bk)%fv(ie,j,k)%qold(1)  = pout/(gamma-1.d0) + 0.5d0
        endif
      enddo
   enddo

!  residue file
  open(unit=15,file='res.dat')
  write(15,*) '#residue'
  if(ieuler == 0 .and. turb == 0 )  write(15,*)'# navier-stokes computation for laminar flow'
  if(ieuler == 0 .and. turb == 1 )  then
     write(15,*)'# navier-stokes computation for turbulent flow'
     write(15,*)'# Baldwin-Lomax model '
  endif

  if(ieuler == 1 ) write(15,*)'# euler computation'
  write(15,'(1x, "# mach number =", f8.4)')m_inf
  write(15,'(1x, "# aoa         =",f8.4)')aoa
  write(15,'(1x, "# p_inf,t_inf =",2(f8.4,2x))')p_inf,t_inf
  write(15,'(1x, "# Re,Pr,Pr_t  =",e15.2,2x,2(f8.4,2x))')re_l,pr_l,pr_t
  close(15)

else

!  restart file (copy q_o.dat to q_i.dat for restart)
  inquire(file="q_i.dat", exist=file_exists) ! file_exists will be true if the file
  if(.not.file_exists) then
     print*
     print*,'Copying q_o.dat to q_i.dat for restart...'
     call system('cp -fr q_i.dat  q_i.dat.old')
     call system('cp -fr q_o.dat  q_i.dat')
  endif   

  open(unit=12,file='q_i.dat',form='unformatted')
     read(12) initer
  do bk = 1,nblk
     call start_end(bk)
     read(12) ((( blk(bk)%fv(i,j,k)%qvar(jterm),i=1,cx),j=1,cy),jterm=1,ndim+2)
  enddo
  close(12)
  print*
  print*,'Restaring from iteration',initer

endif

!visturb(:,:)=0.0d0
  do bk = 1,nblk
     call start_end(bk)
     do j = js,je
        do i = is,ie
           blk(bk)%fv(i,j,k)%mu=0.0d0
           blk(bk)%fv(i,j,k)%mut=0.0d0
        enddo 
     enddo 
  enddo 

  do bk = 1,nblk
     call start_end(bk)
     do j = js,je
        do i = is,ie
           blk(bk)%fv(i,j,k)%qold(:)=blk(bk)%fv(i,j,k)%qvar(:)
        enddo
     enddo
  enddo


!=======================
end subroutine intialize
!=======================
