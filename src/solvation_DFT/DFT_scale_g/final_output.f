      SUBROUTINE final_output
      use rho
	use solute_kind
	use control
      IMPLICIT NONE

      include 'param.h'
      include 'system.h'
      include 'numbers.h'
      include 'fft.h'
      include 'c.h'
      include 'rho.h'
      include 'iteration.h'
      include 'bridge.h'


      INTEGER,parameter:: n_gr=400
      integer i,j,k,nr,n
      real*8 rnr,delta_gr
      REAL*8 X, Y, Z, RR
      real*8 gro(n_gr),grh(n_gr),ngr(n_gr)
      real*8 ftot,xc,yc,zc,rrc,sr,ftemp,chg,ftotx,ftoty,ftotz
	  real*8 dft_peak,sc_rate1,sc_rate2


	  if(one_step)then
	    out_dens=.false.
	  else
	    if(scale_md)then
	      final=.false.
	    endif
	  endif
	  
      delta_gr=0.2

      do i=1,n_gr
        gro(i)=0.0D0
        grh(i)=0.0D0
        ngr(i)=0.0D0
      enddo
      

      if(out_dens.and.final)then
      open(300,file='./density_output.dat')      
      open(302,file='./gr_O.dat')
      open(303,file='./gr_H.dat')
      else
      open(302,file='./gr_O.dat')
      open(303,file='./gr_H.dat')
      endif
     
      ftotx=0.0D0
      ftoty=0.0D0
      ftotz=0.0D0

      DO i = 1, nfft1
        DO j = 1, nfft2
          DO k = 1, nfft3
       
           x = (i-1)*Lg/FLOAT(nfft1)-lg/2.0
           y = (j-1)*Lg/FLOAT(nfft2)-lg/2.0
           z = (k-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = DSQRT(x**2 + y**2 + z**2)
       
           if(out_dens.and.final)then
      write(300,*)x,y,z,rho_O_old(i,j,k),rho_H_old(i,j,k)
           endif

           nr=int(rr/delta_gr)+1
           rnr=nr*1.0D0-rr/delta_gr   

           gro(nr)=gro(nr)+rnr*rho_O_old(i,j,k)
           gro(nr+1)=gro(nr+1)+(1.0D0-rnr)*rho_O_old(i,j,k)
           grh(nr)=grh(nr)+rnr*rho_H_old(i,j,k)
           grh(nr+1)=grh(nr+1)+(1.0D0-rnr)*rho_H_old(i,j,k)
           ngr(nr)=ngr(nr)+rnr 
           ngr(nr+1)=ngr(nr+1)+(1.0D0-rnr)       

           do n=1,nb_solute_sites
             xc=x - x_mol(n)+lg/2.0
             yc=y - y_mol(n)+lg/2.0
             zc=z - z_mol(n)+lg/2.0
             rc=dsqrt(xc*xc+yc*yc+zc*zc)
             if(rc.gt.1.0)then
               sr=sig_mol(n)/rc
               sr=sr**6
               ftemp=24*eps_mol(n)*(2*sr**2-sr)*rho_O_old(i,j,k)
               chg=(rho_O_old(i,j,k)-rho_H_old(i,j,k))*chg_solv(1)
               ftemp=ftemp+QFACT*chg_mol(n)*chg/rc
               ftemp=ftemp/rc
               ftotx=ftotx+ftemp*xc/rc
               ftoty=ftoty+ftemp*yc/rc
               ftotz=ftotz+ftemp*zc/rc
             endif
             
           enddo           

          enddo
        enddo
      enddo

      ftotx=ftotx*DeltaV
      ftoty=ftoty*DeltaV
      ftotz=ftotz*DeltaV
      ftot=dsqrt(ftotx**2+ftoty**2+ftotz**2)
      write(*,*)"ftot=",ftot


      dft_peak=0.0
	  do i=1,n_gr
        rr=(i-1)*delta_gr
        if(rr.gt.(lg/2.0-4.0))then
          gro(i)=1.0D0
          grh(i)=1.0D0
          write(302,*)rr,gro(i)
          write(303,*)rr,grh(i)
        elseif(ngr(i).gt.0.0D0)then
          gro(i)=gro(i)/ngr(i)
          grh(i)=grh(i)/ngr(i)
          write(302,*)rr,gro(i)
          write(303,*)rr,grh(i)
		  if(dft_peak.lt.gro(i))then
		    dft_peak=gro(i)
		  endif
        endif
      enddo

	  write(*,*)"peak of gO(r)", dft_peak
	  
      if(out_dens.and.final)then
      close(300)
      endif

      close(302)
      close(303)
 
      if(final)then
	open(400,ACCESS='APPEND',file='../tip3p.dat')
	write(400,*) solute_name, re_sfe*8.31*0.3/4.183,fmx
	close(400)
	endif

      if(final)then
      write(*,*)'final output finished'
      else
      write(*,*)'output finished'
      endif

	  if(scale_md.and.(.not.one_step))then
	    write(*,*)"Scaling"
	    if(dft_peak.lt.1.0)then
		   write(*,*)"Scale is not available",dft_peak
		else
		  sc_rate1=(md_peak-1.0)/(dft_peak-1.0)
	      sc_rate2=log((md_peak-1.0)/(dft_peak-1.0))/dft_peak
		  DO i = 1, nfft1
            DO j = 1, nfft2
              DO k = 1, nfft3
			    x=rho_O_old(i,j,k)
				if(x.lt.1.0)then
				  rho_O_old(i,j,k)=1.0+(x-1.0)*exp(sc_rate2*x)
				else
				  rho_O_old(i,j,k)=1.0+(x-1.0)*sc_rate1
				endif
				x=rho_H_old(i,j,k)
				if(x.lt.1.0)then
				  rho_H_old(i,j,k)=1.0+(x-1.0)*exp(sc_rate2*x)
				else
				  rho_H_old(i,j,k)=1.0+(x-1.0)*sc_rate1
				endif
			  enddo
			enddo
		  enddo
		endif
	  endif
	  
	  
	  
      return
      end     
     
	 
	 
	 
	
	  
	  
	  Subroutine free_space
      use rho
      use vext
      use pe



      DEALLOCATE(rho_O_old)      
      DEALLOCATE(rho_H_old)
      DEALLOCATE(rho_O_new)
      DEALLOCATE(rho_H_new)


      DEALLOCATE(fai_pe)
      DEALLOCATE(bd_ax)
      DEALLOCATE(bd_bx)
      DEALLOCATE(bd_ay)
      DEALLOCATE(bd_by)
      DEALLOCATE(bd_az)
      DEALLOCATE(bd_bz)
      DEALLOCATE(dpar)
      DEALLOCATE(fai_pe1D)
      DEALLOCATE(ele_rho)


      DEALLOCATE(VextO)
      DEALLOCATE(VextH)
      DEALLOCATE(VLJO)
      DEALLOCATE(VintO)
      DEALLOCATE(VintH)

      return
      end
