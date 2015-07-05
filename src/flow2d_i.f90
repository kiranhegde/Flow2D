!================
program  flow_2d
!================
use flux
use vis
use st_end
use pri
use inf
!use ifport ! for secnds() function
!$ use omp_lib
implicit none

!2-d finite volume solver, cell centred approach

!===============================================================
!n    = total no. of cells in xi direction
!m    = total no. of cells in eta direction
!ndim = space dimension (2 for this code)
!===============================================================
integer::core
integer(kind=i2):: i,j,k,bk,niter,iter,irks,ii
integer(kind=i2):: saveinterval
integer(kind=i2):: istart
!integer(kind=i2):: jterm,io 
!real(kind=dp):: time_begin,time_end
!real(kind=dp):: xx,yy,x0,y0,delta,r1,r2
real(kind=dp):: cfl_max,cl,cd
real(kind=dp):: cfl,rk_fac(3)
real(kind=dp):: resi
real(kind=sp):: t1,t2
!real(kind=dp):: ts,te
real(kind=dp):: time_s,time_e
real(kind=dp):: ds,nx,ny,x1,y1


logical :: file_exists
!CHARACTER*100 BUFFER

irks=1
core=1
!CALL GETARG(1,BUFFER,io)
!if(io==-1) then
!   print*
!   print*,"The input format is ..."
!   print*, " ./mis no_of_cores"
!   print*, "Exiting....."
!   print*
!   stop
!endif
!READ(BUFFER,*)core 

!!$OMP PARALLEL      &
!!$OMP DEFAULT(none) &
!!$OMP SHARED(Nthreads)
!!$  Nthreads = OMP_GET_NUM_THREADS()
!!$OMP END PARALLEL


saveinterval=10
initer=0
pi = 4_dp*datan(1.d0)

gamma     = 1.4_dp
gas_const = 287.06_dp

rk_fac(1) = 0.25d0
rk_fac(2) = 0.50d0
rk_fac(3) = 1.00d0
!.......Ficticious cells based on Boundary conditions
bcc(0)=2 ; bcc(1)=1
bcc(2)=1 ; bcc(3)=1
bcc(4)=1 ; bcc(5)=1
bcc(6)=1 ; bcc(7)=1
bcc(8)=0

!$ write ( *, '(a,i8)' ) '  The number of processors available = ', omp_get_num_procs ( )
!$ write ( *, '(a,i8)' ) '  The number of threads available    = ', omp_get_max_threads ( )
!$ CALL OMP_SET_NUM_THREADS(core)
!$ write ( *, '(a,i8)' ) '  The number of processors available = ', omp_get_num_procs ( )
!$ write ( *, '(a,i8)' ) '  The number of threads available    = ', omp_get_max_threads ( )

!....read simulation and  solver input
call read_input(niter,cfl_max,istart)

!....read grid data
call read_mesh

!....intialise solution
call intialize(istart)

!....calculate cell areas 
call cell_area

!....calculate face normals
call normal_zhi
call normal_eta
k=1  



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

!if (nblk>1) call block_conectivity
call block_conectivity


!....BC locations
bk=1
call start_end(bk)
i  = is
open(3,file='bc_1.dat')
do j = js,je+one
   write(3,*)blk(bk)%mesh(i,j,k)%x,blk(bk)%mesh(i,j,k)%y
enddo
close(3)

i  = ie+one
open(3,file='bc_2.dat')
do j = js,je+one
   write(3,*)blk(bk)%mesh(i,j,k)%x,blk(bk)%mesh(i,j,k)%y
enddo
close(3)

j  = js-one
open(3,file='bc_3.dat')
do i = is,ie+one
   write(3,*)blk(bk)%mesh(i,j,k)%x,blk(bk)%mesh(i,j,k)%y
enddo
close(3)

j  = je+one
open(3,file='bc_4.dat')
do i = is,ie+one
   write(3,*)blk(bk)%mesh(i,j,k)%x,blk(bk)%mesh(i,j,k)%y
enddo
close(3)



!=====================================================
call cpu_time(time_s)
t1 = secnds(0.0)


k=1
open(unit=15,file='res.dat',position='append')

!$ ts = omp_get_wtime()
!..start of iterations
do iter = initer+one,niter
   iter2=iter
   cfl = 0.5d0*cfl_max*(1+dtanh(iter*7.d0/500.0 - 3.5))
   !cfl = 0.5d0*cfl_max*(1+dtanh(iter*7.d0/500.0 - 9.0))


!   if(nblk>1) then
      do bk=1,nblk
         call blockinterface(bk)
         call blockinterface(bk)
      enddo
!   endif


   do bk=1,nblk

      call start_end(bk)

!.....dt calculation
      call stability(bk,cfl)

!.....save old  

!$OMP PARALLEL  DEFAULT(NONE) & 
!$OMP FIRSTPRIVATE(bk,k) &
!$OMP PRIVATE(i,j,ii) &
!$OMP SHARED(cx,cy,blk)
!$OMP DO  
      do j = 1,cy
         do i = 1,cx
         do ii=1,ndim+2   
            blk(bk)%fv(i,j,k)%qold(ii) = blk(bk)%fv(i,j,k)%qvar(ii) 
            blk(bk)%fv(i,j,k)%dflux(ii)= 0.0d0
         enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL


!  do irks = 1,3

      call BC_del_flux(bk) 
      call implicit_solver(bk)
!     call walldatas(bk)

!  enddo !r-k

   enddo  


!..residue calculation
   !call residue(iter,cfl,resi)
   call residue(resi)
   write(15,*) iter, dlog10(resi),cfl

   if(mod(iter,saveinterval) .eq. 0)then
    t2 = secnds(t1)
    write(*,*) 'Time passed [s]',t2
    t1 = secnds(0.0)
    call clcd(cl,cd)
    call write_result(iter,cfl,cl,cd,resi)
    inquire(file="stop.run", exist=file_exists) ! file_exists will be true if the file

    if(file_exists) then
       print*," Stopping  at", iter, "iterations..."
       print*," Please remove the file stop.run to restart"
       print*
       print*
       call clcd(cl,cd)
       call write_result(iter,cfl,cl,cd,resi)
       stop
    endif
   endif

enddo

close(15)

!$ write(6,*) 'Time taken (parallel) ', omp_get_wtime()-ts
call cpu_time(time_e)
t2 = secnds(t1)
write(*,*) 'Time taken',time_e-time_s
write(*,*) 'passed [s]',t2


call clcd(cl,cd)
call write_result(iter-1,cfl,cl,cd,resi)

open(3,file='BC02.nrm')
do bk=1,nblk
   call start_end(bk)
     j=js
  do i = is,cx
    x1 = 0.5d0*(blk(bk)%mesh(i+1,j,k)%x + blk(bk)%mesh(i,j,k)%x)
    y1 = 0.5d0*(blk(bk)%mesh(i+1,j,k)%y + blk(bk)%mesh(i,j,k)%y)

    nx=blk(bk)%eta(i,j,k)%nx
    ny=blk(bk)%eta(i,j,k)%ny 
    ds=dsqrt(nx*nx+ny*ny)
    write(3,100)x1,y1,nx,ny
  enddo
enddo
close(3)

open(3,file='BC04.nrm')
do bk=1,nblk
   call start_end(bk)
     j=je+one
  do i = is,cx
    x1 = 0.5d0*(blk(bk)%mesh(i+1,j,k)%x + blk(bk)%mesh(i,j,k)%x)
    y1 = 0.5d0*(blk(bk)%mesh(i+1,j,k)%y + blk(bk)%mesh(i,j,k)%y)

    nx=blk(bk)%eta(i,j,k)%nx
    ny=blk(bk)%eta(i,j,k)%ny 
    ds=dsqrt(nx*nx+ny*ny)
    write(3,100)x1,y1,nx,ny
  enddo
enddo
close(3)

open(3,file='BC01.nrm')
do bk=1,nblk
   call start_end(bk)
     i=is
  do j = js,cy
    x1 = 0.5d0*(blk(bk)%mesh(i,j+1,k)%x + blk(bk)%mesh(i,j,k)%x)
    y1 = 0.5d0*(blk(bk)%mesh(i,j+1,k)%y + blk(bk)%mesh(i,j,k)%y)

    nx=blk(bk)%zhi(i,j,k)%nx
    ny=blk(bk)%zhi(i,j,k)%ny 
    ds=dsqrt(nx*nx+ny*ny)
    write(3,100)x1,y1,nx,ny
  enddo
enddo
close(3)


open(3,file='BC03.nrm')
do bk=1,nblk
   call start_end(bk)
     i=ie+one
  do j = js,cy
    x1 = 0.5d0*(blk(bk)%mesh(i,j+1,k)%x + blk(bk)%mesh(i,j,k)%x)
    y1 = 0.5d0*(blk(bk)%mesh(i,j+1,k)%y + blk(bk)%mesh(i,j,k)%y)

    nx=blk(bk)%zhi(i,j,k)%nx
    ny=blk(bk)%zhi(i,j,k)%ny 
    ds=dsqrt(nx*nx+ny*ny)
    write(3,100)x1,y1,nx,ny
  enddo
enddo
close(3)
100 format(1x,4(f15.6,1x))
!101 format(1x,2(f15.6,1x))





!====================
end program  flow_2d
!====================


!=============================================
subroutine residue(resi1)
!=============================================
use flux
use inf
use st_end
implicit none

integer(kind=i2):: i,j,k,bk
real(kind=dp):: resi1,resi2,dt,total  
  k=1
  resi1 = 0.d0
  total = 0.d0

  do bk = 1,nblk
     call start_end(bk)

!$OMP PARALLEL  DEFAULT(NONE) & 
!$OMP FIRSTPRIVATE(bk,k) &
!$OMP PRIVATE(i,j,dt,resi2) &
!$OMP SHARED(is,ie,js,je,blk) &
!$OMP REDUCTION(+:resi1) 
!$OMP DO
  do j = js,je
    do i = is,ie
         ! Check density
      if (blk(bk)%fv(i,j,k)%qvar(2)<=0.d0) then
          write(*,107)i,j,blk(bk)%fv(i,j,k)%qvar(2)
          blk(bk)%fv(i,j,k)%qvar(2) = dmin1(0.1_dp*M_inf**2, 1.0e-04_dp)
          !call tecplot
          !  pause
          !stop
      endif
      resi2 = blk(bk)%fv(i,j,k)%qvar(2)/blk(bk)%fv(i,j,k)%qold(2) - 1.d0
      dt=blk(bk)%fv(i,j,k)%dt
      resi1 = resi1 + resi2*resi2/(dt*dt)
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

      total=total+(ie-is+1)*(je-js+1) 
  enddo

  resi1 = dsqrt(resi1)/total

!.........write the residue
!print *, iter, dlog10(resi1),cfl
107 format(1x,2(i4,1x),2x,f10.5)

!=====================
end subroutine residue 
!=====================
