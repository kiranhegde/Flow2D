!========================================================
subroutine BC_del_flux(bk)
!========================================================
!  Applying BC and computing interface fluxes
!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!nx    = same as n for this code
!ny    = same as n for this code
!n_zhi = normals on the xi faces
!n_eta = normals on the eta faces
!qvar  = conserved variable q
!fin_zhi = fluxes on the zhi faces (output)
!fin_eta = fluxes on the eta faces (output)
!x = x coordinates
!y = y coordinates
!area = cell areas
!fintr = flux at the given cell interface
!===============================================================
use flux
use vis
use st_end

implicit none
integer(kind=i2):: bk,i,j,k,idir,ii
real(kind=dp):: fintr(ndim+2)

real(kind=dp):: aml(ndim),amr(ndim),amfr(ndim)
real(kind=dp):: qfleft(ndim+2),qleft(ndim+2),qright(ndim+2),qfright(ndim+2)
real(kind=dp):: r1,r2,r3,r4
k=1

call start_end(bk)


if(ieuler .eq. 0) then

!.........invicid boundary condition
  call boundary(bk)

!.........viscous boundary condition
  call boundary_vis(bk)

!.........invicid & viscous flux calculation 

     idir = 1
     !       xi direction fluxes
     do i = is,ie-one
       do j = js,je!

          it=i ;jt=j 
          aml(1)  = blk(bk)%zhi(i  ,j,k)%nx
          aml(2)  = blk(bk)%zhi(i  ,j,k)%ny
          amr(1)  = blk(bk)%zhi(i+1,j,k)%nx
          amr(2)  = blk(bk)%zhi(i+1,j,k)%ny
          amfr(1) = blk(bk)%zhi(i+2,j,k)%nx
          amfr(2) = blk(bk)%zhi(i+2,j,k)%ny
          
          qfleft(1:ndim+2)  = blk(bk)%fv(i-1,j,k)%qvar(1:ndim+2)
          qleft(1:ndim+2)   = blk(bk)%fv(i  ,j,k)%qvar(1:ndim+2)
          qright(1:ndim+2)  = blk(bk)%fv(i+1,j,k)%qvar(1:ndim+2)
          qfright(1:ndim+2) = blk(bk)%fv(i+2,j,k)%qvar(1:ndim+2)
          
          fintr(:)=0.d0  
          call cell_dist(bk,cx,cy,i,j,one,two,r1,r2,r3,r4)
          call fluxinv_int(cx,cy,one,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)
            
          blk(bk)%fv(i,j,k)%dflux(1:ndim+2) = blk(bk)%fv(i,j,k)%dflux(1:ndim+2) + fintr(1:ndim+2)
          blk(bk)%fv(i+1,j,k)%dflux(1:ndim+2) = blk(bk)%fv(i+1,j,k)%dflux(1:ndim+2) - fintr(1:ndim+2)
           
          fintr(:)=0.d0  
          call fluxvis_int_xi(bk,i,j,idir,qleft,qright,amr,aml,fintr)
          blk(bk)%fv(i,j,k)%dflux(1:ndim+2) = blk(bk)%fv(i,j,k)%dflux(1:ndim+2) - fintr(1:ndim+2)
          blk(bk)%fv(i+1,j,k)%dflux(1:ndim+2) = blk(bk)%fv(i+1,j,k)%dflux(1:ndim+2) + fintr(1:ndim+2)
     
       enddo
     enddo
     
     idir = 2
     !eta direction fluxes
     
     do j = js,je-one
       do i = is,ie
       aml(1)  = blk(bk)%eta(i  ,j,k)%nx
       aml(2)  = blk(bk)%eta(i  ,j,k)%ny
       amr(1)  = blk(bk)%eta(i,j+1,k)%nx
       amr(2)  = blk(bk)%eta(i,j+1,k)%ny
       amfr(1)  = blk(bk)%eta(i,j+2,k)%nx
       amfr(2)  = blk(bk)%eta(i,j+2,k)%ny
       
       qfleft(1:ndim+2)  = blk(bk)%fv(i,j-1,k)%qvar(1:ndim+2)
       qleft(1:ndim+2)   = blk(bk)%fv(i  ,j,k)%qvar(1:ndim+2)
       qright(1:ndim+2)  = blk(bk)%fv(i,j+1,k)%qvar(1:ndim+2)
       qfright(1:ndim+2) = blk(bk)%fv(i,j+2,k)%qvar(1:ndim+2)
       
       fintr(:)=0.d0  
       call cell_dist(bk,cx,cy,i,j,two,two,r1,r2,r3,r4)
       call fluxinv_int(cx,cy,two,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)
       
       blk(bk)%fv(i,j,k)%dflux(1:ndim+2) = blk(bk)%fv(i,j,k)%dflux(1:ndim+2) + fintr(1:ndim+2)
       blk(bk)%fv(i,j+1,k)%dflux(1:ndim+2) = blk(bk)%fv(i,j+1,k)%dflux(1:ndim+2) - fintr(1:ndim+2)
     
       fintr(:)=0.d0  
       call fluxvis_int_eta(bk,i,j,idir,qleft,qright,amr,aml,fintr) 
       blk(bk)%fv(i,j,k)%dflux(1:ndim+2) = blk(bk)%fv(i,j,k)%dflux(1:ndim+2) - fintr(1:ndim+2)
       blk(bk)%fv(i,j+1,k)%dflux(1:ndim+2) = blk(bk)%fv(i,j+1,k)%dflux(1:ndim+2) + fintr(1:ndim+2)
       enddo
     enddo

!.........computing turbulent viscosity 
     if(turb>0)  call turbulencemodel(bk) 

else
!.........invicid boundary condition
  call boundary(bk)

!.........invicid flux calculation 
!       xi  direction fluxes
     do i = is,ie-one
       do j = js,je!

          it=i ;jt=j 
          aml(1)  = blk(bk)%zhi(i  ,j,k)%nx
          aml(2)  = blk(bk)%zhi(i  ,j,k)%ny
          amr(1)  = blk(bk)%zhi(i+1,j,k)%nx
          amr(2)  = blk(bk)%zhi(i+1,j,k)%ny
          amfr(1) = blk(bk)%zhi(i+2,j,k)%nx
          amfr(2) = blk(bk)%zhi(i+2,j,k)%ny
          
          qfleft( 1:ndim+2) = blk(bk)%fv(i-1,j,k)%qvar(1:ndim+2)
          qleft(  1:ndim+2) = blk(bk)%fv(i  ,j,k)%qvar(1:ndim+2)
          qright( 1:ndim+2) = blk(bk)%fv(i+1,j,k)%qvar(1:ndim+2)
          qfright(1:ndim+2) = blk(bk)%fv(i+2,j,k)%qvar(1:ndim+2)
          do ii=1,4
             if(isnan(qfleft(ii)).or.isnan(qleft(ii)).or.isnan(qright(ii)).or.isnan(qfright(ii))) then
               print*,'Xi flux',iter2
               print '(2(i3,1x),4(f15.8,1x))',i,j,fintr(1:4)
               print '(2(i3,1x),4(f15.8,1x))',i,j,qfleft(1:4)
               print '(2(i3,1x),4(f15.8,1x))',i,j,qfleft(1:4)
               print '(2(i3,1x),4(f15.8,1x))',i,j,qleft(1:4)
               print '(2(i3,1x),4(f15.8,1x))',i,j,qright(1:4)
               print '(2(i3,1x),4(f15.8,1x))',i,j,qfright(1:4)
               !print '(2(i3,1x),4(f15.8,1x))',it,jt,pr-pl,fa,c12,ro12
               !print '(1x,"Q-",6(f15.8,2x))',rol,ul,vl,hl,pl,cl 
               !print '(1x,"Q+",6(f15.8,2x))',ror,ur,vr,hr,pr,cr 
               !call get_cell(bk)
               stop
             endif
          enddo


          
          fintr(:)=0.d0  
          call cell_dist(bk,cx,cy,i,j,one,two,r1,r2,r3,r4)
          call fluxinv_int(cx,cy,one,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)
          
          blk(bk)%fv(i  ,j,k)%dflux(1:ndim+2) = blk(bk)%fv(i  ,j,k)%dflux(1:ndim+2) + fintr(1:ndim+2)
          blk(bk)%fv(i+1,j,k)%dflux(1:ndim+2) = blk(bk)%fv(i+1,j,k)%dflux(1:ndim+2) - fintr(1:ndim+2)
          

       enddo
     enddo

!       eta direction fluxes
     do j = js,je-one
       do i = is,ie
          aml(1)  = blk(bk)%eta(i,j  ,k)%nx
          aml(2)  = blk(bk)%eta(i,j  ,k)%ny
          amr(1)  = blk(bk)%eta(i,j+1,k)%nx
          amr(2)  = blk(bk)%eta(i,j+1,k)%ny
          amfr(1) = blk(bk)%eta(i,j+2,k)%nx
          amfr(2) = blk(bk)%eta(i,j+2,k)%ny
          
          qfleft( 1:ndim+2) = blk(bk)%fv(i,j-1,k)%qvar(1:ndim+2)
          qleft(  1:ndim+2) = blk(bk)%fv(i,j  ,k)%qvar(1:ndim+2)
          qright( 1:ndim+2) = blk(bk)%fv(i,j+1,k)%qvar(1:ndim+2)
          qfright(1:ndim+2) = blk(bk)%fv(i,j+2,k)%qvar(1:ndim+2)
             if(isnan(qfleft(ii)).or.isnan(qleft(ii)).or.isnan(qright(ii)).or.isnan(qfright(ii))) then
               print*,'Eta flux',iter2
               print '(2(i3,1x),4(f15.8,1x))',i,j,fintr(1:4)
               print '(2(i3,1x),4(f15.8,1x))',i,j,qfleft(1:4)
               print '(2(i3,1x),4(f15.8,1x))',i,j,qfleft(1:4)
               print '(2(i3,1x),4(f15.8,1x))',i,j,qleft(1:4)
               print '(2(i3,1x),4(f15.8,1x))',i,j,qright(1:4)
               print '(2(i3,1x),4(f15.8,1x))',i,j,qfright(1:4)
               !print '(2(i3,1x),4(f15.8,1x))',it,jt,pr-pl,fa,c12,ro12
               !print '(1x,"Q-",6(f15.8,2x))',rol,ul,vl,hl,pl,cl 
               !print '(1x,"Q+",6(f15.8,2x))',ror,ur,vr,hr,pr,cr 
               !call get_cell(bk)
               stop
             endif
          
          fintr(:)=0.d0  
          call cell_dist(bk,cx,cy,i,j,two,two,r1,r2,r3,r4)
          call fluxinv_int(cx,cy,two,r1,r2,r3,r4,aml,amr,qleft,qright,qfleft,qfright,fintr)
          
          blk(bk)%fv(i,j  ,k)%dflux(1:ndim+2) = blk(bk)%fv(i,j  ,k)%dflux(1:ndim+2) + fintr(1:ndim+2)
          blk(bk)%fv(i,j+1,k)%dflux(1:ndim+2) = blk(bk)%fv(i,j+1,k)%dflux(1:ndim+2) - fintr(1:ndim+2)
       enddo
     enddo
endif

!===========================
end subroutine BC_del_flux
!===========================
