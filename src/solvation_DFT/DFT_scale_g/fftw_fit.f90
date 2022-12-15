         subroutine fftw_3d_for(nx,ny,nz,input,output)
         use, intrinsic::iso_c_binding
         include '/home/jiafu/bin/fftw333/include/fftw3.f'
         integer nx,ny,nz
         complex(C_DOUBLE_COMPLEX) input(nx,ny,nz),output(nx,ny,nz)
         type(C_PTR)::plan         
         call dfftw_plan_dft_3d(plan,nz,ny,nx,input,output,FFTW_FORWARD,FFTW_ESTIMATE)
         call dfftw_execute_dft(plan,input,output)         
         call dfftw_destroy_plan(plan)
         return
         end
         
         subroutine fftw_3d_bak(nx,ny,nz,input,output)
         use, intrinsic::iso_c_binding
         include '/home/jiafu/bin/fftw333/include/fftw3.f'
         integer nx,ny,nz,i,j,k
         complex(C_DOUBLE_COMPLEX) input(nx,ny,nz),output(nx,ny,nz)
         type(C_PTR)::plan
         real*8 vxyz
         vxyz=nx*ny*nz*1.0D0         
         call dfftw_plan_dft_3d(plan,nz,ny,nx,input,output,FFTW_BACKWARD,FFTW_ESTIMATE)
         call dfftw_execute_dft(plan,input,output)
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  output(i,j,k)=output(i,j,k)/vxyz
               enddo
            enddo
         enddo         
         call dfftw_destroy_plan(plan)
         return
         end

