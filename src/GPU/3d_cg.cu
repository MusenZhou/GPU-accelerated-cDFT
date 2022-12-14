#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cufft.h>
#include<time.h>

#define PI 3.141592653589793
#define running_block_size 128



int *k_gas;
int *na, *nb, *nc;
int *na_device, *nb_device, *nc_device;
double *la_expanded, *lb_expanded, *lc_expanded;

double *grand_pontential, *grand_pontential_device;
double *Vext_device;
double *chem_device;
double *rho_bulk_device;

double *sqrt_dens_device;
double *gradient_device;
double *diameter_HS_device;
int *k_gas_device;
double *w_n0_FFT_device;
double *w_n3_FFT_device;
double *w_nv2x_FFT_device;
double *w_nv2y_FFT_device;
double *w_nv2z_FFT_device;
double *uatt_device;



cufftHandle fft_plan_device;
cufftDoubleComplex *FFT_in_device;
cufftDoubleComplex *n0_FFT_in_device, *n1_FFT_in_device, *n2_FFT_in_device, *n3_FFT_in_device;
cufftDoubleComplex *nv1x_FFT_in_device, *nv1y_FFT_in_device, *nv1z_FFT_in_device;
cufftDoubleComplex *nv2x_FFT_in_device, *nv2y_FFT_in_device, *nv2z_FFT_in_device;





cufftDoubleComplex *FFT_out_device;

cufftDoubleComplex *n0_FFT_out_device, *n1_FFT_out_device, *n2_FFT_out_device, *n3_FFT_out_device;
cufftDoubleComplex *nv1x_FFT_out_device, *nv1y_FFT_out_device, *nv1z_FFT_out_device;
cufftDoubleComplex *nv2x_FFT_out_device, *nv2y_FFT_out_device, *nv2z_FFT_out_device;


double *denskr_device, *denski_device;





// __global__ defines the funciton that can be called from the host (CPU) and executed in the device (GPU)
__global__
void check_int(int n, int *x)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i, j;
	for (j=0; j<n; j++)
	{
		for (i=index; i<n; i+=stride)
		{
			if (i==j)
			{
				//printf("index: %d\t%d\n", i, x[i]);
				printf("%d %d\n", i, x[i]);
			}
		}
	}
}

__global__
void check_double(int n, double *x)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i, j;
	for (j=0; j<n; j++)
	{
		for (i=index; i<n; i+=stride)
		{
			if (i==j)
			{
				printf("index: %d\t%lf\n", i, x[i]);
			}
		}
	}
}




__global__
void check_complex(int n, cufftDoubleComplex *x)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i, j;
	for (j=0; j<n; j++)
	{
		for (i=index; i<n; i+=stride)
		{
			if (i==j)
			{
				printf("index: %d\t%lf\t%lf\n", i, x[i].x, x[i].y);
			}
		}
	}
}



























__global__
void cal_weight_density(int *na_device, int *nb_device, int *nc_device, 
	double *la_expanded_device, double *lb_expanded_device, double *lc_expanded_device,
	int *k_gas_device, double *diameter_HS_device,
	double *w_n0_FFT_device, double *w_n3_FFT_device, double *w_nv2x_FFT_device, double *w_nv2y_FFT_device, double *w_nv2z_FFT_device)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	int index_gas, index_x, index_y, index_z;
	double fft_xk,  fft_yk, fft_zk;
	double dx, dy, dz;
	double fft_k;
	double fft_abs_kR;

	for (i=index; i<(na_device[0]*nb_device[0]*nc_device[0]*k_gas_device[0]); i+=stride)
	{
		dx = 1.0*la_expanded_device[0]/na_device[0];
		dy = 1.0*lb_expanded_device[0]/nb_device[0];
		dz = 1.0*lc_expanded_device[0]/nc_device[0];

		// decompose the row, column information from the overal index
		index_gas = (int) ( (i)/(na_device[0]*nb_device[0]*nc_device[0]) );
		index_x = (int) ( (i-index_gas*na_device[0]*nb_device[0]*nc_device[0])/(nb_device[0]*nc_device[0]) );
		index_y = (int) ( (i-index_gas*na_device[0]*nb_device[0]*nc_device[0]-index_x*nb_device[0]*nc_device[0])/(nc_device[0]) );
		index_z = (int) (i-index_gas*na_device[0]*nb_device[0]*nc_device[0]-index_x*nb_device[0]*nc_device[0]-index_y*nc_device[0]);

		// the first point (origin) is unique
		if ((i%(na_device[0]*nb_device[0]*nc_device[0])) == 0)
		{
			w_n0_FFT_device[i] = 1.0 / (na_device[0]*nb_device[0]*nc_device[0]);
			w_n3_FFT_device[i] = 1.0*PI*pow(diameter_HS_device[index_gas],3)/6/(na_device[0]*nb_device[0]*nc_device[0]);
			w_nv2x_FFT_device[i] = 0;
			w_nv2y_FFT_device[i] = 0;
			w_nv2z_FFT_device[i] = 0;
		}
		else
		{
			// reformat the index into fft space with fftw fashion
			if (index_x<=0.5*na_device[0])
			{
				fft_xk = 2.0*PI*index_x/na_device[0]/dx;
			}
			else
			{
				fft_xk = -2.0*PI*(na_device[0]-index_x)/na_device[0]/dx;
			}
			if (index_y<=0.5*nb_device[0])
			{
				fft_yk = 2.0*PI*index_y/nb_device[0]/dy;
			}
			else
			{
				fft_yk = -2.0*PI*(nb_device[0]-index_y)/nb_device[0]/dy;
			}
			if (index_z<=0.5*nc_device[0])
			{
				fft_zk = 2.0*PI*index_z/nc_device[0]/dz;
			}
			else
			{
				fft_zk = -2.0*PI*(nc_device[0]-index_z)/nc_device[0]/dz;
			}

			fft_k = sqrt(pow(fft_xk,2)+pow(fft_yk,2)+pow(fft_zk,2));
			fft_abs_kR = 1.0*fabs(fft_k)*diameter_HS_device[index_gas]/2.0;

			w_n0_FFT_device[i] = sin(fft_abs_kR)/fft_abs_kR / (na_device[0]*nb_device[0]*nc_device[0]);
			w_n3_FFT_device[i] = 4.0*PI*(sin(fft_abs_kR)-fft_abs_kR*cos(fft_abs_kR))/pow(fabs(fft_k), 3) / 
			(na_device[0]*nb_device[0]*nc_device[0]);
			w_nv2x_FFT_device[i] = -1.0*w_n3_FFT_device[i]*fft_xk;
			w_nv2y_FFT_device[i] = -1.0*w_n3_FFT_device[i]*fft_yk;
			w_nv2z_FFT_device[i] = -1.0*w_n3_FFT_device[i]*fft_zk;
		}
	}
}












__global__
void cal_Vext(int *na_device, int *nb_device, int *nc_device, 
	double *la_expanded_device, double *lb_expanded_device, double *lc_expanded_device, 
	int *N_atoms_expanded_device, double *epsilon_host_star_expanded_device, double *sigma_host_expanded_device,
	double *x_host_expanded_device, double *y_host_expanded_device, double *z_host_expanded_device,
	int *k_gas_device, double *epsilon_star_device, double *sigma_device, double *rho_bulk_device, 
	double *diameter_HS_device, double *sqrt_dens_device, double *cutoff_device, double *Vext_device)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i, ii;
	int index_gas, index_x, index_y, index_z;
	double pos_x, pos_y, pos_z;
	double dis_x, dis_y, dis_z;
	double dis;
	double epsilon_mix, sigma_mix;
	double pot, pot_cutoff;
	double dens;
	// printf("%d\n", na_device[0]*nb_device[0]*nc_device[0]*k_gas_device[0]);
	for (i=index; i<(na_device[0]*nb_device[0]*nc_device[0]*k_gas_device[0]); i+=stride)
	{
		// decompose the row, column information from the overal index
		index_gas = (int) ( (i)/(na_device[0]*nb_device[0]*nc_device[0]) );
		index_x = (int) ( (i-index_gas*na_device[0]*nb_device[0]*nc_device[0])/(nb_device[0]*nc_device[0]) );
		index_y = (int) ( (i-index_gas*na_device[0]*nb_device[0]*nc_device[0]-index_x*nb_device[0]*nc_device[0])/(nc_device[0]) );
		index_z = (int) (i-index_gas*na_device[0]*nb_device[0]*nc_device[0]-index_x*nb_device[0]*nc_device[0]-index_y*nc_device[0]);

		Vext_device[i] = 0;

		pos_x = 1.0*index_x*la_expanded_device[0]/na_device[0];
		pos_y = 1.0*index_y*lb_expanded_device[0]/nb_device[0];
		pos_z = 1.0*index_z*lc_expanded_device[0]/nc_device[0];

		for (ii=0; ii<N_atoms_expanded_device[0]; ii++)
		{
			dis_x = pos_x - x_host_expanded_device[ii];
			dis_y = pos_y - y_host_expanded_device[ii];
			dis_z = pos_z - z_host_expanded_device[ii];
			if (dis_x > (0.5*la_expanded_device[0]))
			{
				dis_x = dis_x - la_expanded_device[0];
			}
			else if (dis_x < (-0.5*la_expanded_device[0]))
			{
				dis_x = dis_x + la_expanded_device[0];
			}

			if (dis_y > (0.5*lb_expanded_device[0]))
			{
				dis_y = dis_y - lb_expanded_device[0];
			}
			else if (dis_y < (-0.5*lb_expanded_device[0]))
			{
				dis_y = dis_y + lb_expanded_device[0];
			}

			if (dis_z > (0.5*lc_expanded_device[0]))
			{
				dis_z = dis_z - lc_expanded_device[0];
			}
			else if (dis_z < (-0.5*lc_expanded_device[0]))
			{
				dis_z = dis_z + lc_expanded_device[0];
			}
			dis = sqrt(pow(dis_x,2)+pow(dis_y,2)+pow(dis_z,2));
			if (dis<cutoff_device[0])
			{
				sigma_mix = 1.0*(sigma_device[index_gas]+sigma_host_expanded_device[ii])/2.0;
				epsilon_mix = sqrt(epsilon_star_device[index_gas]*epsilon_host_star_expanded_device[ii]);
				if (dis < 0.1*sigma_mix)
				{
					dis = 0.1*sigma_mix;
				}
				pot = 4.0*epsilon_mix*(pow((1.0*sigma_mix/dis),12)-pow((1.0*sigma_mix/dis),6));
				pot_cutoff = 4.0*epsilon_mix*(pow((1.0*sigma_mix/cutoff_device[0]),12)-pow((1.0*sigma_mix/cutoff_device[0]),6));
				Vext_device[i] = Vext_device[i] + (pot - pot_cutoff);
			}
		}

		dens = rho_bulk_device[index_gas]*exp(-Vext_device[i])/100;
		// dens = rho_bulk_device[index_gas]/100;


		if (dens > 2.0/pow(diameter_HS_device[index_gas],3)/k_gas_device[0])
		{
			dens = 2.0/pow(diameter_HS_device[index_gas],3)/k_gas_device[0];
		}

		sqrt_dens_device[i] = sqrt(dens);
	}
}











__global__
void ini_device_value_double(int n, double *x)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	for (i=index; i<n; i+=stride)
	{
		x[i] = 0;
	}
}



__global__
void ini_device_value_complex(int n, cufftDoubleComplex *x)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	for (i=index; i<n; i+=stride)
	{
		x[i].x = 0;
		x[i].y = 0;
	}
}



__global__
void real_split_into_complex_batch(int n, double *in, cufftDoubleComplex *out, int num_batch)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	int index_in;

	for (i=index; i<n; i+=stride)
	{
		index_in = num_batch*n + i;
		out[i].x = in[index_in]*in[index_in];
		out[i].y = 0;
	}
}

__global__
void real_split_into_complex_separate_batch(int n, double *in_r, double *in_i, cufftDoubleComplex *out, int num_batch)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	int index_in;

	for (i=index; i<n; i+=stride)
	{
		index_in = num_batch*n + i;
		out[i].x = in_r[index_in];
		out[i].y = in_i[index_in];
	}
}

__global__
void complext_split_into_real_separate_batch(int n, double *out_r, double *out_i, cufftDoubleComplex *in, int num_batch)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	int index_in;

	for (i=index; i<n; i+=stride)
	{
		index_in = num_batch*n + i;
		out_r[index_in] = in[i].x;
		out_i[index_in] = in[i].y;
	}
}






__global__
void add_weighted_density(int *na_device, int *nb_device, int *nc_device, double *diameter_HS_device, int *k_gas_device, 
	cufftDoubleComplex *FFT_out_device, cufftDoubleComplex *n0_FFT_in_device, cufftDoubleComplex *n1_FFT_in_device, 
	cufftDoubleComplex *n2_FFT_in_device, cufftDoubleComplex *n3_FFT_in_device, 
	cufftDoubleComplex *nv1x_FFT_in_device, cufftDoubleComplex *nv1y_FFT_in_device, cufftDoubleComplex *nv1z_FFT_in_device, 
	cufftDoubleComplex *nv2x_FFT_in_device, cufftDoubleComplex *nv2y_FFT_in_device, cufftDoubleComplex *nv2z_FFT_in_device, 
	double *w_n0_FFT_device, double *w_n3_FFT_device, double *w_nv2x_FFT_device, double *w_nv2y_FFT_device, double *w_nv2z_FFT_device, 
	double *uatt_device, double *denskr_device, double  *denski_device,	int j)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i, jj;
	int batch_i;
	int temp_index1, temp_index2, temp_index;
	for (i=index; i<(na_device[0]*nb_device[0]*nc_device[0]); i+=stride)
	{
		batch_i = i + j*na_device[0]*nb_device[0]*nc_device[0];

		n0_FFT_in_device[i].x = n0_FFT_in_device[i].x + FFT_out_device[i].x*w_n0_FFT_device[batch_i];
		n0_FFT_in_device[i].y = n0_FFT_in_device[i].y + FFT_out_device[i].y*w_n0_FFT_device[batch_i];
		n1_FFT_in_device[i].x = n1_FFT_in_device[i].x + FFT_out_device[i].x*w_n0_FFT_device[batch_i]*diameter_HS_device[j]*0.5;
		n1_FFT_in_device[i].y = n1_FFT_in_device[i].y + FFT_out_device[i].y*w_n0_FFT_device[batch_i]*diameter_HS_device[j]*0.5;
		n2_FFT_in_device[i].x = n2_FFT_in_device[i].x + FFT_out_device[i].x*w_n0_FFT_device[batch_i]*PI*diameter_HS_device[j]*diameter_HS_device[j];
		n2_FFT_in_device[i].y = n2_FFT_in_device[i].y + FFT_out_device[i].y*w_n0_FFT_device[batch_i]*PI*diameter_HS_device[j]*diameter_HS_device[j];
		n3_FFT_in_device[i].x = n3_FFT_in_device[i].x + FFT_out_device[i].x*w_n3_FFT_device[batch_i];
		n3_FFT_in_device[i].y = n3_FFT_in_device[i].y + FFT_out_device[i].y*w_n3_FFT_device[batch_i];

		nv1x_FFT_in_device[i].x = nv1x_FFT_in_device[i].x + -FFT_out_device[i].y*w_nv2x_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);
		nv1x_FFT_in_device[i].y = nv1x_FFT_in_device[i].y + FFT_out_device[i].x*w_nv2x_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);
		nv1y_FFT_in_device[i].x = nv1y_FFT_in_device[i].x + -FFT_out_device[i].y*w_nv2y_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);
		nv1y_FFT_in_device[i].y = nv1y_FFT_in_device[i].y + FFT_out_device[i].x*w_nv2y_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);
		nv1z_FFT_in_device[i].x = nv1z_FFT_in_device[i].x + -FFT_out_device[i].y*w_nv2z_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);
		nv1z_FFT_in_device[i].y = nv1z_FFT_in_device[i].y + FFT_out_device[i].x*w_nv2z_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);

		nv2x_FFT_in_device[i].x = nv2x_FFT_in_device[i].x + -FFT_out_device[i].y*w_nv2x_FFT_device[batch_i];
		nv2x_FFT_in_device[i].y = nv2x_FFT_in_device[i].y + FFT_out_device[i].x*w_nv2x_FFT_device[batch_i];
		nv2y_FFT_in_device[i].x = nv2y_FFT_in_device[i].x + -FFT_out_device[i].y*w_nv2y_FFT_device[batch_i];
		nv2y_FFT_in_device[i].y = nv2y_FFT_in_device[i].y + FFT_out_device[i].x*w_nv2y_FFT_device[batch_i];
		nv2z_FFT_in_device[i].x = nv2z_FFT_in_device[i].x + -FFT_out_device[i].y*w_nv2z_FFT_device[batch_i];
		nv2z_FFT_in_device[i].y = nv2z_FFT_in_device[i].y + FFT_out_device[i].x*w_nv2z_FFT_device[batch_i];

		for (jj=0; jj<k_gas_device[0]; jj++)
		{
			if (jj>j)
			{
				temp_index1 = jj;
				temp_index2 = j;
			}
			else
			{
				temp_index1 = j;
				temp_index2 = jj;
			}
			temp_index = (int) floor((temp_index1+0)*(temp_index1+1)*0.5)*na_device[0]*nb_device[0]*nc_device[0] 
			+ temp_index2*na_device[0]*nb_device[0]*nc_device[0] + i;

			denskr_device[i+jj*na_device[0]*nb_device[0]*nc_device[0]] += FFT_out_device[i].x*uatt_device[temp_index];
			denski_device[i+jj*na_device[0]*nb_device[0]*nc_device[0]] += FFT_out_device[i].y*uatt_device[temp_index];
		}
	}
}





__global__
void cal_F_deri(int *na_device, int *nb_device, int *nc_device, 
	cufftDoubleComplex *n0_FFT_out_device,
	cufftDoubleComplex *n1_FFT_out_device, cufftDoubleComplex *n2_FFT_out_device, cufftDoubleComplex *n3_FFT_out_device, 
	cufftDoubleComplex *nv1x_FFT_out_device, cufftDoubleComplex *nv1y_FFT_out_device, cufftDoubleComplex *nv1z_FFT_out_device, 
	cufftDoubleComplex *nv2x_FFT_out_device, cufftDoubleComplex *nv2y_FFT_out_device, cufftDoubleComplex *nv2z_FFT_out_device,
	cufftDoubleComplex *n0_FFT_in_device,
	cufftDoubleComplex *n1_FFT_in_device, cufftDoubleComplex *n2_FFT_in_device, cufftDoubleComplex *n3_FFT_in_device, 
	cufftDoubleComplex *nv1x_FFT_in_device, cufftDoubleComplex *nv1y_FFT_in_device, cufftDoubleComplex *nv1z_FFT_in_device, 
	cufftDoubleComplex *nv2x_FFT_in_device, cufftDoubleComplex *nv2y_FFT_in_device, cufftDoubleComplex *nv2z_FFT_in_device)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	double tn0, tn1, tn2, tn3;
	// double tnv1x, tnv1y, tnv1z;
	double tnv2x, tnv2y, tnv2z;
	double tn30, tn31, tn32, tn33, tnv21, tnv22, t2n3, t3n3;

	for (i=index; i<(na_device[0]*nb_device[0]*nc_device[0]); i+=stride)
	{
		tn0 = n0_FFT_out_device[i].x;
		tn1 = n1_FFT_out_device[i].x;
		tn2 = n2_FFT_out_device[i].x;
		tn3 = n3_FFT_out_device[i].x;

		if (tn0<0)
		{
			tn0 = 0;
			// if (tn0<-0.0001)
			// {
			// 	// printf("Fatal Error!!!!!!! tn0!!!!!!!\n");
			// }
		}

		if (tn3<0)
		{
			tn3 = 0;
			// if (tn3<-0.0001)
			// {
			// 	// printf("Fatal Error!!!!!!! tn3!!!!!!!\n");
			// }
		}
		else if (tn3>0.99)
		{
			tn3 = 0.99;
			// printf("Fatal Error!!!!!!! tn3!!!!!!!\n");
		}

		tnv2x = nv2x_FFT_out_device[i].x;
		tnv2y = nv2y_FFT_out_device[i].x;
		tnv2z = nv2z_FFT_out_device[i].x;

		if (tn3<1e-5)
		{
			n0_FFT_in_device[i].x = 0;
			n0_FFT_in_device[i].y = 0;
			n1_FFT_in_device[i].x = 0;
			n1_FFT_in_device[i].y = 0;
			n2_FFT_in_device[i].x = 0;
			n2_FFT_in_device[i].y = 0;
			n3_FFT_in_device[i].x = 0;
			n3_FFT_in_device[i].y = 0;

			nv1x_FFT_in_device[i].x = 0;
			nv1x_FFT_in_device[i].y = 0;
			nv1y_FFT_in_device[i].x = 0;
			nv1y_FFT_in_device[i].y = 0;
			nv1z_FFT_in_device[i].x = 0;
			nv1z_FFT_in_device[i].y = 0;

			nv2x_FFT_in_device[i].x = 0;
			nv2x_FFT_in_device[i].y = 0;
			nv2y_FFT_in_device[i].x = 0;
			nv2y_FFT_in_device[i].y = 0;
			nv2z_FFT_in_device[i].x = 0;
			nv2z_FFT_in_device[i].y = 0;
		}
		else
		{
			tn31 = 1-tn3;
            tn30 = log(tn31);
            tn31 = 1.0/tn31;
            tn32 = tn31*tn31;
            tn33 = tn31*tn32;
            tnv22 = tnv2x*tnv2x + tnv2y*tnv2y + tnv2z*tnv2z;
            tnv21 = tnv2x*nv1x_FFT_out_device[i].x + tnv2y*nv1y_FFT_out_device[i].x + tnv2z*nv1z_FFT_out_device[i].x;
            t2n3 = 1.0/tn3/tn3;
            t3n3 = t2n3/tn3;
            // printf("%.5e\t%.5e\n", nv1y_FFT_out_device[i].x, (1.0/(1-tn3)));
            
            // nv2y_FFT_in_device[i].x = -1.0*(tn30/tn3+tn32)*tn2*tnv2y/(6.0*PI*tn3) - tnv1y*tn31;


            n0_FFT_in_device[i].x = -tn30;
			n0_FFT_in_device[i].y = 0;
			n1_FFT_in_device[i].x = tn2*tn31;
			n1_FFT_in_device[i].y = 0;
			n2_FFT_in_device[i].x = tn1*tn31 + 1.0*(1.0*tn30/tn3+tn32)*(tn2*tn2-tnv22)/(12.0*PI*tn3);
			n2_FFT_in_device[i].y = 0;
			n3_FFT_in_device[i].x = ((1.0*tn30*t3n3/(18.0*PI))+(1.0-3*tn3+1.0/tn32)*t2n3*tn33/(36.0*PI))*(3*tn2*tnv22-tn2*tn2*tn2)
	            					+ tn0*tn31 + (tn1*tn2-tnv21)*tn32 + 80.0*exp((tn3-0.7)*80.0);
			n3_FFT_in_device[i].y = 0;

			nv1x_FFT_in_device[i].x = -tnv2x*tn31;
			nv1x_FFT_in_device[i].y = 0;
			nv1y_FFT_in_device[i].x = -tnv2y*tn31;
			nv1y_FFT_in_device[i].y = 0;
			nv1z_FFT_in_device[i].x = -tnv2z*tn31;
			nv1z_FFT_in_device[i].y = 0;

			nv2x_FFT_in_device[i].x = -1.0*(tn30/tn3+tn32)*tn2*tnv2x/(6.0*PI*tn3) - nv1x_FFT_out_device[i].x*tn31;
			nv2x_FFT_in_device[i].y = 0;
			nv2y_FFT_in_device[i].x = -1.0*(log(1-tn3)/tn3+tn32)*tn2*tnv2y/(6.0*PI*tn3) - nv1y_FFT_out_device[i].x*(1.0/(1-n3_FFT_out_device[i].x));
			// nv2y_FFT_in_device[i].x = -1.0*(tn30/tn3+tn32)*tn2*tnv2y/(6.0*PI*tn3) - tnv1y*tn31;
			nv2y_FFT_in_device[i].y = 0;
			nv2z_FFT_in_device[i].x = -1.0*(log(1-tn3)/tn3+tn32)*tn2*tnv2z/(6.0*PI*tn3) - nv1z_FFT_out_device[i].x*(1.0/(1-n3_FFT_out_device[i].x));
			// nv2z_FFT_in_device[i].x = -1.0*(tn30/tn3+tn32)*tn2*tnv2z/(6.0*PI*tn3) - tnv1z*tn31;
			nv2z_FFT_in_device[i].y = 0;
		}
	}
}





__global__
void sum_F_deri(int *na_device, int *nb_device, int *nc_device, double *diameter_HS_device, int *k_gas_device, 
	cufftDoubleComplex *FFT_in_device, cufftDoubleComplex *n0_FFT_out_device, 
	cufftDoubleComplex *n1_FFT_out_device, cufftDoubleComplex *n2_FFT_out_device, cufftDoubleComplex *n3_FFT_out_device, 
	cufftDoubleComplex *nv1x_FFT_out_device, cufftDoubleComplex *nv1y_FFT_out_device, cufftDoubleComplex *nv1z_FFT_out_device, 
	cufftDoubleComplex *nv2x_FFT_out_device, cufftDoubleComplex *nv2y_FFT_out_device, cufftDoubleComplex *nv2z_FFT_out_device, 
	double *w_n0_FFT_device, double *w_n3_FFT_device, double *w_nv2x_FFT_device, double *w_nv2y_FFT_device, double *w_nv2z_FFT_device, 
	int j)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	int batch_i;
	for (i=index; i<(na_device[0]*nb_device[0]*nc_device[0]); i+=stride)
	{
		batch_i = i + j*na_device[0]*nb_device[0]*nc_device[0];

		FFT_in_device[i].x = n0_FFT_out_device[i].x*w_n0_FFT_device[batch_i];
		FFT_in_device[i].y = n0_FFT_out_device[i].y*w_n0_FFT_device[batch_i];

		FFT_in_device[i].x += n1_FFT_out_device[i].x*w_n0_FFT_device[batch_i]*diameter_HS_device[j]*0.5;
		FFT_in_device[i].y += n1_FFT_out_device[i].y*w_n0_FFT_device[batch_i]*diameter_HS_device[j]*0.5;

		FFT_in_device[i].x += n2_FFT_out_device[i].x*w_n0_FFT_device[batch_i]*PI*diameter_HS_device[j]*diameter_HS_device[j];
		FFT_in_device[i].y += n2_FFT_out_device[i].y*w_n0_FFT_device[batch_i]*PI*diameter_HS_device[j]*diameter_HS_device[j];

		FFT_in_device[i].x += n3_FFT_out_device[i].x*w_n3_FFT_device[batch_i];
		FFT_in_device[i].y += n3_FFT_out_device[i].y*w_n3_FFT_device[batch_i];



		FFT_in_device[i].x += nv1x_FFT_out_device[i].y*w_nv2x_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);
		FFT_in_device[i].y += -nv1x_FFT_out_device[i].x*w_nv2x_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);
		FFT_in_device[i].x += nv1y_FFT_out_device[i].y*w_nv2y_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);
		FFT_in_device[i].y += -nv1y_FFT_out_device[i].x*w_nv2y_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);
		FFT_in_device[i].x += nv1z_FFT_out_device[i].y*w_nv2z_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);
		FFT_in_device[i].y += -nv1z_FFT_out_device[i].x*w_nv2z_FFT_device[batch_i]*0.5/(PI*diameter_HS_device[j]);



		FFT_in_device[i].x += nv2x_FFT_out_device[i].y*w_nv2x_FFT_device[batch_i];
		FFT_in_device[i].y += -nv2x_FFT_out_device[i].x*w_nv2x_FFT_device[batch_i];
		FFT_in_device[i].x += nv2y_FFT_out_device[i].y*w_nv2y_FFT_device[batch_i];
		FFT_in_device[i].y += -nv2y_FFT_out_device[i].x*w_nv2y_FFT_device[batch_i];
		FFT_in_device[i].x += nv2z_FFT_out_device[i].y*w_nv2z_FFT_device[batch_i];
		FFT_in_device[i].y += -nv2z_FFT_out_device[i].x*w_nv2z_FFT_device[batch_i];
	}
}





void new_MBWR(int *k_gas, double *epsilon, double *sigma, double *temperature, double *rho_bulk, double *chem, double *bulk_pressure)
{
    int i, ii;

    double *x, *a, *b, *c, *d, *G, *gamma, *F;
    x = (double *) malloc(sizeof(double)*32);
    a = (double *) malloc(sizeof(double)*8);
    b = (double *) malloc(sizeof(double)*6);
    c = (double *) malloc(sizeof(double)*8);
    d = (double *) malloc(sizeof(double)*6);
    G = (double *) malloc(sizeof(double)*6);
    gamma = (double *) malloc(sizeof(double));
    F = (double *) malloc(sizeof(double));
    double *sum_a, *sum_b;
    sum_a = (double *) malloc(sizeof(double));
    sum_b = (double *) malloc(sizeof(double));

    double *rho_r, *T_r;
    rho_r = (double *) malloc(sizeof(double));
    T_r = (double *) malloc(sizeof(double));

    double *F_r, *P_r, *U_r;
    F_r = (double *) malloc(sizeof(double));
    P_r = (double *) malloc(sizeof(double));
    U_r = (double *) malloc(sizeof(double));


    // mixing LJ parameter for mixture
    double *epsilon_mix, *sigma_mix;
    double *tot_rho;
    double *sigma_x3;
    double *epsilon_x;
    epsilon_mix = (double *) malloc(sizeof(double));
    sigma_mix = (double *) malloc(sizeof(double));
    tot_rho = (double *) malloc(sizeof(double));
    sigma_x3 = (double *) malloc(sizeof(double));
    epsilon_x = (double *) malloc(sizeof(double));
    sigma_x3[0] = 0;
    epsilon_x[0] = 0;
    tot_rho[0] = 0;
    for (i=0; i<k_gas[0]; i++)
    {
        tot_rho[0] = tot_rho[0] + rho_bulk[i];
        for (ii=0; ii<k_gas[0]; ii++)
        {
            epsilon_mix[0] = sqrt(epsilon[i]*epsilon[ii]);
            sigma_mix[0] = 0.5*(sigma[i]+sigma[ii]);
            sigma_x3[0] = sigma_x3[0] + rho_bulk[i]*rho_bulk[ii]*pow(sigma_mix[0],3);
            epsilon_x[0] = epsilon_x[0] + rho_bulk[i]*rho_bulk[ii]*epsilon_mix[0]*pow(sigma_mix[0],3);
        }
    }
    sigma_x3[0] = 1.0*sigma_x3[0]/tot_rho[0]/tot_rho[0];
    epsilon_x[0] = 1.0*epsilon_x[0]/tot_rho[0]/tot_rho[0]/sigma_x3[0];

    rho_r[0] = sigma_x3[0]*tot_rho[0];
    T_r[0] = 1.0*temperature[0]/epsilon_x[0];

    // MBWR parameter from Keith E. Gubbins
    x[0] = 0.8623085097507421;
    x[1] = 2.976218765822098;
    x[2] = -8.402230115796038;
    x[3] = 0.1054136629203555;
    x[4] = -0.8564583828174598;
    x[5] = 1.582759470107601;
    x[6] = 0.7639421948305453;
    x[7] = 1.753173414312048;
    x[8] = 2.798291772190376e+03;
    x[9] = -4.8394220260857657e-02;
    x[10] = 0.9963265197721935;
    x[11] = -3.698000291272493e+01;
    x[12] = 2.084012299434647e+01;
    x[13] = 8.305402124717285e+01;
    x[14] = -9.574799715203068e+02;
    x[15] = -1.477746229234994e+02;
    x[16] = 6.398607852471505e+01;
    x[17] = 1.603993673294834e+01;
    x[18] = 6.805916615864377e+01;
    x[19] = -2.791293578795945e+03;
    x[20] = -6.245128304568454;
    x[21] = -8.116836104958410e+03;
    x[22] = 1.488735559561229e+01;
    x[23] = -1.059346754655084e+04;
    x[24] = -1.131607632802822e+02;
    x[25] = -8.867771540418822e+03;
    x[26] = -3.986982844450543e+01;
    x[27] = -4.689270299917261e+03;
    x[28] = 2.593535277438717e+02;
    x[29] = -2.694523589434903e+03;
    x[30] = -7.218487631550215e+02;
    x[31] = 1.721802063863269e+02;
    // MBWR parameters of a
    a[0] = x[0]*T_r[0] + x[1]*pow(T_r[0], 0.5) + x[2] + x[3]/T_r[0] + x[4]/pow(T_r[0], 2);
    a[1] = x[5]*T_r[0] + x[6] + x[7]/T_r[0] + x[8]/pow(T_r[0], 2);
    a[2] = x[9]*T_r[0] + x[10] + x[11]/T_r[0];
    a[3] = x[12];
    a[4] = x[13]/T_r[0] + x[14]/pow(T_r[0], 2);
    a[5] = x[15]/T_r[0];
    a[6] = x[16]/T_r[0] + x[17]/pow(T_r[0], 2);
    a[7] = x[18]/pow(T_r[0], 2);
    // MBWR parameters of b
    b[0] = x[19]/pow(T_r[0], 2) + x[20]/pow(T_r[0], 3);
    b[1] = x[21]/pow(T_r[0], 2) + x[22]/pow(T_r[0], 4);
    b[2] = x[23]/pow(T_r[0], 2) + x[24]/pow(T_r[0], 3);
    b[3] = x[25]/pow(T_r[0], 2) + x[26]/pow(T_r[0], 4);
    b[4] = x[27]/pow(T_r[0], 2) + x[28]/pow(T_r[0], 3);
    b[5] = x[29]/pow(T_r[0], 2) + x[30]/pow(T_r[0], 3) + x[31]/pow(T_r[0], 4);
    // MBWR parameters of c
    c[0] = x[1]*sqrt(T_r[0])*0.5 + x[2] + 2*x[3]/T_r[0] + 3*x[4]/pow(T_r[0],2);
    c[1] = x[6] + 2*x[7]/T_r[0] + 3*x[8]/pow(T_r[0],2);
    c[2] = x[10] + 2*x[11]/T_r[0];
    c[3] = x[12];
    c[4] = 2*x[13]/T_r[0] + 3*x[14]/pow(T_r[0],2);
    c[5] = 2*x[15]/T_r[0];
    c[6] = 2*x[16]/T_r[0] + 3*x[17]/pow(T_r[0],2);
    c[7] = 3*x[18]/pow(T_r[0],2);
    // MBWE parameters of d
    d[0] = 3*x[19]/pow(T_r[0],2) + 4*x[20]/pow(T_r[0],3);
    d[1] = 3*x[21]/pow(T_r[0],2) + 5*x[22]/pow(T_r[0],4);
    d[2] = 3*x[23]/pow(T_r[0],2) + 4*x[24]/pow(T_r[0],3);
    d[3] = 3*x[25]/pow(T_r[0],2) + 5*x[26]/pow(T_r[0],4);
    d[4] = 3*x[27]/pow(T_r[0],2) + 4*x[28]/pow(T_r[0],3);
    d[5] = 3*x[29]/pow(T_r[0],2) + 4*x[30]/pow(T_r[0],3) + 5*x[31]/pow(T_r[0],4);
    // MBWR parameter of G
    gamma[0] = 3;
    F[0] = exp(-gamma[0]*pow(rho_r[0],2));
    G[0] = (1-F[0])/(2*gamma[0]);
    for (i=1; i<6; i++)
    {
        G[i] = i/gamma[0]*G[i-1] - F[0]/(2*gamma[0])*pow(rho_r[0],2*i);
    }



    
    // calculate reduced residual Helmholtz Free energy
    sum_a[0] = 0;
    sum_b[0] = 0;
    for (i=0; i<=7; i++)
    {
        sum_a[0] = sum_a[0] + a[i]*pow(rho_r[0], (i+1))/(i+1);
    }
    for (i=0; i<=5; i++)
    {
        sum_b[0] = sum_b[0] + b[i]*G[i];
    }
    F_r[0] = sum_a[0] + sum_b[0];
    // calculate reduced reduced bulk pressure
    sum_a[0] = 0;
    sum_b[0] = 0;
    for (i=0; i<=7; i++)
    {
        sum_a[0] = sum_a[0] + a[i]*pow(rho_r[0], (i+2));
    }
    for (i=0; i<=5; i++)
    {
        sum_b[0] = sum_b[0] + F[0]*b[i]*pow(rho_r[0], (2*i+3));
    }
    P_r[0] = rho_r[0]*T_r[0] + sum_a[0] + sum_b[0];
    bulk_pressure[0] = P_r[0]*epsilon_x[0]*1.38064852*pow(10,-23)/(sigma_x3[0]*pow(10,-30))*pow(10,-5);
    // calculate reduced internal energy
    sum_a[0] = 0;
    sum_b[0] = 0;
    for (i=0; i<=7; i++)
    {
        sum_a[0] = sum_a[0] + c[i]*pow(rho_r[0], (i+1))/(i+1);
    }
    for (i=0; i<=5; i++)
    {
        sum_b[0] = sum_b[0] + d[i]*G[i];
    }
    U_r[0] = sum_a[0] + sum_b[0];


    // calculate the chemical potential for each species
    double *deriv1, *deriv2, *deriv3;
    deriv1 = (double *) malloc(sizeof(double));
    deriv2 = (double *) malloc(sizeof(double));
    deriv3 = (double *) malloc(sizeof(double));
    for (i=0; i<k_gas[0]; i++)
    {
        deriv3[0] = 0;
        for (ii=0; ii<k_gas[0]; ii++)
        {
            sigma_mix[0] = 0.5*(sigma[i]+sigma[ii]);
            deriv3[0] = deriv3[0] + rho_bulk[ii]*pow(sigma_mix[0],3);
        }
        deriv3[0] = deriv3[0]/tot_rho[0];
        deriv3[0] = deriv3[0] - sigma_x3[0];
        deriv3[0] = 2.0*deriv3[0]/tot_rho[0];

        deriv1[0] = 0;
        for (ii=0; ii<k_gas[0]; ii++)
        {
            epsilon_mix[0] = sqrt(epsilon[i]*epsilon[ii]);
            sigma_mix[0] = 0.5*(sigma[i]+sigma[ii]);
            deriv1[0] = deriv1[0] + rho_bulk[ii]*epsilon_mix[0]*pow(sigma_mix[0],3);
        }
        deriv1[0] = deriv1[0]/tot_rho[0]/sigma_x3[0] - epsilon_x[0];
        deriv1[0] = 2.0*deriv1[0]/tot_rho[0];
        deriv1[0] = -epsilon_x[0]*deriv3[0]/sigma_x3[0] + deriv1[0];

        deriv2[0] = (1.0*P_r[0]/pow(rho_r[0],2) - T_r[0]/rho_r[0]) * (sigma_x3[0]+tot_rho[0]*deriv3[0]);
        deriv2[0] = deriv2[0] - (F_r[0]-U_r[0])*deriv1[0]/epsilon_x[0];

        chem[i] = (F_r[0]*epsilon_x[0] + F_r[0]*tot_rho[0]*deriv1[0] + tot_rho[0]*epsilon_x[0]*deriv2[0])/temperature[0] + log(rho_bulk[i]);
        // chem[i] = (F_r[0]*epsilon_x[0] + F_r[0]*tot_rho[0]*deriv1[0] + tot_rho[0]*epsilon_x[0]*deriv2[0]) + log(rho_bulk[i]);
    }
}
























__global__
void cal_lj_pre_freq_batch(int *na_device, int *nb_device, int *nc_device, double *la_expanded_device, double *lb_expanded_device, 
	double *lc_expanded_device, double *sigma_device, double *epsilon_star_device, double *cutoff_device, int j, int jj, 
	cufftDoubleComplex *FFT_in_device)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	int index_x, index_y, index_z;
	double dx, dy, dz;
	double x, y, z;
	double sigma_mix, epsilon_mix;
	double dis;
	double pot_cutoff, pot;

	for (i=index; i<(na_device[0]*nb_device[0]*nc_device[0]); i+=stride)
	{
		dx = 1.0*la_expanded_device[0]/na_device[0];
		dy = 1.0*lb_expanded_device[0]/nb_device[0];
		dz = 1.0*lc_expanded_device[0]/nc_device[0];

		// decompose the row, column information from the overall index
		index_x = (int) ( (i)/(nb_device[0]*nc_device[0]) );
		index_y = (int) ( (i-index_x*nb_device[0]*nc_device[0])/(nc_device[0]) );
		index_z = (int) ( i-index_x*nb_device[0]*nc_device[0]-index_y*nc_device[0] );

		sigma_mix = 0.5*(sigma_device[j]+sigma_device[jj]);
		epsilon_mix = sqrt((epsilon_star_device[j]*epsilon_star_device[jj]));

		pot_cutoff = 4.0*epsilon_mix* ( pow((1.0*sigma_mix/cutoff_device[0]),12) - pow((1.0*sigma_mix/cutoff_device[0]),6) );

		if (index_x<0.5*na_device[0])
		{
			x = index_x*dx;
		}
		else
		{
			x = (index_x-na_device[0])*dx;
		}
		if (index_y<0.5*nb_device[0])
		{
			y = index_y*dy;
		}
		else
		{
			y = (index_y-nb_device[0])*dy;
		}
		if (index_z<0.5*nc_device[0])
		{
			z = index_z*dz;
		}
		else
		{
			z = (index_z-nc_device[0])*dz;
		}

		dis = sqrt(pow(x,2)+pow(y,2)+pow(z,2));

		if (dis < sigma_mix)
		{
			pot = 0;
		}
		else if  (dis>cutoff_device[0])
		{
			pot = 0;
		}
		else
		{
			pot = 4.0*epsilon_mix*(pow((1.0*sigma_mix/dis),12)-pow((1.0*sigma_mix/dis),6)) - pot_cutoff;
		}

		FFT_in_device[i].x = pot;
		FFT_in_device[i].y = 0;
	}
}








__global__
void cal_lj_freq_store_batch(int *na_device, int *nb_device, int *nc_device, double *la_expanded_device, 
	double *lb_expanded_device, double *lc_expanded_device, int j, int jj, cufftDoubleComplex *FFT_out_device, 
	double *uatt_device)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	int nijk_gas;
	double dx, dy, dz;

	for (i=index; i<(na_device[0]*nb_device[0]*nc_device[0]); i+=stride)
	{
		dx = 1.0*la_expanded_device[0]/na_device[0];
		dy = 1.0*lb_expanded_device[0]/nb_device[0];
		dz = 1.0*lc_expanded_device[0]/nc_device[0];

		nijk_gas = i + jj*na_device[0]*nb_device[0]*nc_device[0] + floor((j+0)*(j+1)*0.5)*na_device[0]*nb_device[0]*nc_device[0];
		uatt_device[nijk_gas] = 1.0*FFT_out_device[i].x*dx*dy*dz/(na_device[0]*nb_device[0]*nc_device[0]);
	}
}





__global__
void cal_grand_potential(int *na_device, int *nb_device, int *nc_device, int *k_gas_device, double *Vext_device, 
	double *chem_device,  double *sqrt_dens_device, double *rho_bulk_device, cufftDoubleComplex *n0_FFT_out_device, 
	cufftDoubleComplex *n1_FFT_out_device, cufftDoubleComplex *n2_FFT_out_device, cufftDoubleComplex *n3_FFT_out_device, 
	cufftDoubleComplex *nv1x_FFT_out_device, cufftDoubleComplex *nv1y_FFT_out_device, cufftDoubleComplex *nv1z_FFT_out_device, 
	cufftDoubleComplex *nv2x_FFT_out_device, cufftDoubleComplex *nv2y_FFT_out_device, cufftDoubleComplex *nv2z_FFT_out_device,
	double *denskr_device, double *grand_pontential_device)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i, j;
	int nijk_gas;
	double tn0, tn1, tn2, tn3;
	double tnv1x, tnv1y, tnv1z;
	double tnv2x, tnv2y, tnv2z;
	double tn30, tn31, tn32, tn33, tnv21, tnv22, t2n3, t3n3;
	for (i=index; i<(na_device[0]*nb_device[0]*nc_device[0]); i+=stride)
	{
		grand_pontential_device[i] = 0;
		for (j=0; j<k_gas_device[0]; j++)
		{
			nijk_gas = i + j*na_device[0]*nb_device[0]*nc_device[0];
			grand_pontential_device[i] += (Vext_device[nijk_gas]-chem_device[j])*sqrt_dens_device[nijk_gas]*sqrt_dens_device[nijk_gas];

			if ((sqrt_dens_device[nijk_gas]*sqrt_dens_device[nijk_gas])>1.0e-5*rho_bulk_device[j])
			{
				grand_pontential_device[i] += sqrt_dens_device[nijk_gas]*sqrt_dens_device[nijk_gas]*(log(sqrt_dens_device[nijk_gas]*sqrt_dens_device[nijk_gas])-1);
			}
		}

		tn0 = n0_FFT_out_device[i].x;
		tn1 = n1_FFT_out_device[i].x;
		tn2 = n2_FFT_out_device[i].x;
		tn3 = n3_FFT_out_device[i].x;

		if (tn0<0)
		{
			tn0 = 0;
			// if (tn0<-0.0001)
			// {
			// 	// printf("Fatal Error!!!!!!! tn0!!!!!!!\n");
			// }
		}

		if (tn3<0)
		{
			tn3 = 0;
			// if (tn3<-0.0001)
			// {
			// 	// printf("Fatal Error!!!!!!! tn3!!!!!!!\n");
			// }
		}
		else if (tn3>0.99)
		{
			tn3 = 0.99;
			// printf("Fatal Error!!!!!!! tn3!!!!!!!\n");
		}
		tnv1x = nv1x_FFT_out_device[i].x;
		tnv1y = nv1y_FFT_out_device[i].x;
		tnv1z = nv1z_FFT_out_device[i].x;
		tnv2x = nv2x_FFT_out_device[i].x;
		tnv2y = nv2y_FFT_out_device[i].x;
		tnv2z = nv2z_FFT_out_device[i].x;

		if (tn3>=1e-5)
		{
			tn31 = 1-tn3;
            tn30 = log(tn31);
            tn31 = 1.0/tn31;
            tn32 = tn31*tn31;
            tn33 = tn31*tn32;
            tnv22 = tnv2x*tnv2x + tnv2y*tnv2y + tnv2z*tnv2z;
            tnv21 = tnv2x*tnv1x + tnv2y*tnv1y + tnv2z*tnv1z;
            t2n3 = 1.0/tn3/tn3;
            t3n3 = t2n3/tn3;

            grand_pontential_device[i] += -tn0*tn30;

            grand_pontential_device[i] += (tn1*tn2-tnv21)*tn31;

            grand_pontential_device[i] += (pow(tn2,3)-3*tn2*tnv22)*(tn3+tn30/tn32)/(36.0*PI*tn3*tn3)*tn32;

            grand_pontential_device[i] += exp(80*(tn3-0.7));

            for (j=0; j<k_gas_device[0]; j++)
            {
                nijk_gas = i + j*na_device[0]*nb_device[0]*nc_device[0];
                grand_pontential_device[i] += 0.5*sqrt_dens_device[nijk_gas]*sqrt_dens_device[nijk_gas]*denskr_device[nijk_gas];
            }
		}
	}
}





__global__
void cal_grad(int *na_device, int *nb_device, int *nc_device, double *chem_device, cufftDoubleComplex *FFT_out_device, 
	double *denskr_device, double *Vext_device, double *sqrt_dens_device, double *rho_bulk_device, double *gradient_device, int j)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i;
	int batch_i;
	for (i=index; i<(na_device[0]*nb_device[0]*nc_device[0]); i+=stride)
	{
		batch_i = i + j*na_device[0]*nb_device[0]*nc_device[0];

		if ((sqrt_dens_device[batch_i]*sqrt_dens_device[batch_i])>(1.0e-5*rho_bulk_device[j]))
		{
			gradient_device[batch_i] = log(sqrt_dens_device[batch_i]*sqrt_dens_device[batch_i]) + Vext_device[batch_i] - chem_device[j] + FFT_out_device[i].x + denskr_device[batch_i];
		}
		else
		{
			gradient_device[batch_i] = 0;
		}
		gradient_device[batch_i] = gradient_device[batch_i]*2*sqrt_dens_device[batch_i];
	}
}





double myvalue(double *sqrt_dens, long int system_size)
{
	cudaMemcpy(sqrt_dens_device, sqrt_dens, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0], cudaMemcpyHostToDevice);
	// 
	int i, j;



	// initialize the weighted density for hard sphere term
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n0_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n1_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n2_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n3_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv1x_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv1y_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv1z_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv2x_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv2y_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv2z_FFT_in_device);

	ini_device_value_double<<<(int)((na[0]*nb[0]*nc[0]*k_gas[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0]*k_gas[0], denskr_device);
	ini_device_value_double<<<(int)((na[0]*nb[0]*nc[0]*k_gas[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0]*k_gas[0], denski_device);

	// calculate the weighted density
    for (j=0; j<k_gas[0]; j++)
    {
    	real_split_into_complex_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    	(na[0]*nb[0]*nc[0], sqrt_dens_device, FFT_in_device, j);

    	if (cufftExecZ2Z(fft_plan_device, FFT_in_device, FFT_out_device, CUFFT_FORWARD) != CUFFT_SUCCESS)
		{
			printf("CUFFT error: Forward FFT of density FAILED at specie:\n");
		}

		add_weighted_density<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
		(na_device, nb_device, nc_device, diameter_HS_device, k_gas_device, FFT_out_device, n0_FFT_in_device, 
		n1_FFT_in_device, n2_FFT_in_device, n3_FFT_in_device, nv1x_FFT_in_device, nv1y_FFT_in_device, nv1z_FFT_in_device, 
		nv2x_FFT_in_device, nv2y_FFT_in_device, nv2z_FFT_in_device, w_n0_FFT_device, w_n3_FFT_device, w_nv2x_FFT_device, 
		w_nv2y_FFT_device, w_nv2z_FFT_device, uatt_device, denskr_device, denski_device, j);
    }

    cufftExecZ2Z(fft_plan_device, n0_FFT_in_device, n0_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, n1_FFT_in_device, n1_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, n2_FFT_in_device, n2_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, n3_FFT_in_device, n3_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv1x_FFT_in_device, nv1x_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv1y_FFT_in_device, nv1y_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv1z_FFT_in_device, nv1z_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv2x_FFT_in_device, nv2x_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv2y_FFT_in_device, nv2y_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv2z_FFT_in_device, nv2z_FFT_out_device, CUFFT_INVERSE);

	for (j=0; j<k_gas[0]; j++)
    {
    	real_split_into_complex_separate_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    	(na[0]*nb[0]*nc[0], denskr_device, denski_device, FFT_in_device, j);

    	if (cufftExecZ2Z(fft_plan_device, FFT_in_device, FFT_out_device, CUFFT_INVERSE) != CUFFT_SUCCESS)
		{
			printf("CUFFT error: Backward FFT of density FAILED at specie:\n");
		}

		complext_split_into_real_separate_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
		(na[0]*nb[0]*nc[0], denskr_device, denski_device, FFT_out_device, j);
    }

    // calculate grand potential
    cal_grand_potential<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    (na_device, nb_device, nc_device, k_gas_device, Vext_device, chem_device, sqrt_dens_device, rho_bulk_device, n0_FFT_out_device, 
	n1_FFT_out_device, n2_FFT_out_device, n3_FFT_out_device, nv1x_FFT_out_device, nv1y_FFT_out_device, nv1z_FFT_out_device, 
	nv2x_FFT_out_device, nv2y_FFT_out_device, nv2z_FFT_out_device, denskr_device, grand_pontential_device);

    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
	cudaMemcpy(grand_pontential, grand_pontential_device, sizeof(double)*na[0]*nb[0]*nc[0], cudaMemcpyDeviceToHost);
	double result = 0;
	for (i=0; i<na[0]*nb[0]*nc[0]; i++)
	{
		result += grand_pontential[i];
	}
	result = result * (la_expanded[0]/(na[0])) * (lb_expanded[0]/(nb[0])) * (lc_expanded[0]/(nc[0]));

	// printf("free energy: %lf\n", result);

	return result;
}





void mygrad(double *gradient, double *sqrt_dens, long int system_size)
{
	cudaMemcpy(sqrt_dens_device, sqrt_dens, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0], cudaMemcpyHostToDevice);
	// 
	int i, j;



	// initialize the weighted density for hard sphere term
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n0_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n1_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n2_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n3_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv1x_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv1y_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv1z_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv2x_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv2y_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv2z_FFT_in_device);

	ini_device_value_double<<<(int)((na[0]*nb[0]*nc[0]*k_gas[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0]*k_gas[0], denskr_device);
	ini_device_value_double<<<(int)((na[0]*nb[0]*nc[0]*k_gas[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0]*k_gas[0], denski_device);

	// calculate the weighted density
    for (j=0; j<k_gas[0]; j++)
    {
    	real_split_into_complex_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    	(na[0]*nb[0]*nc[0], sqrt_dens_device, FFT_in_device, j);

    	if (cufftExecZ2Z(fft_plan_device, FFT_in_device, FFT_out_device, CUFFT_FORWARD) != CUFFT_SUCCESS)
		{
			printf("CUFFT error: Forward FFT of density FAILED at specie:\n");
		}

		add_weighted_density<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
		(na_device, nb_device, nc_device, diameter_HS_device, k_gas_device, FFT_out_device, n0_FFT_in_device, 
		n1_FFT_in_device, n2_FFT_in_device, n3_FFT_in_device, nv1x_FFT_in_device, nv1y_FFT_in_device, nv1z_FFT_in_device, 
		nv2x_FFT_in_device, nv2y_FFT_in_device, nv2z_FFT_in_device, w_n0_FFT_device, w_n3_FFT_device, w_nv2x_FFT_device, 
		w_nv2y_FFT_device, w_nv2z_FFT_device, uatt_device, denskr_device, denski_device, j);
    }

    cufftExecZ2Z(fft_plan_device, n0_FFT_in_device, n0_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, n1_FFT_in_device, n1_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, n2_FFT_in_device, n2_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, n3_FFT_in_device, n3_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv1x_FFT_in_device, nv1x_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv1y_FFT_in_device, nv1y_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv1z_FFT_in_device, nv1z_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv2x_FFT_in_device, nv2x_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv2y_FFT_in_device, nv2y_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv2z_FFT_in_device, nv2z_FFT_out_device, CUFFT_INVERSE);

	for (j=0; j<k_gas[0]; j++)
    {
    	real_split_into_complex_separate_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    	(na[0]*nb[0]*nc[0], denskr_device, denski_device, FFT_in_device, j);

    	if (cufftExecZ2Z(fft_plan_device, FFT_in_device, FFT_out_device, CUFFT_INVERSE) != CUFFT_SUCCESS)
		{
			printf("CUFFT error: Backward FFT of density FAILED at specie:\n");
		}

		complext_split_into_real_separate_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
		(na[0]*nb[0]*nc[0], denskr_device, denski_device, FFT_out_device, j);
    }

    //calculate the derivative term of weighted density
    cal_F_deri<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    (na_device, nb_device, nc_device, n0_FFT_out_device, n1_FFT_out_device, n2_FFT_out_device, n3_FFT_out_device, 
	nv1x_FFT_out_device, nv1y_FFT_out_device, nv1z_FFT_out_device, nv2x_FFT_out_device, nv2y_FFT_out_device, nv2z_FFT_out_device,
	n0_FFT_in_device, n1_FFT_in_device, n2_FFT_in_device, n3_FFT_in_device, nv1x_FFT_in_device, nv1y_FFT_in_device, nv1z_FFT_in_device, 
	nv2x_FFT_in_device, nv2y_FFT_in_device, nv2z_FFT_in_device);



	cufftExecZ2Z(fft_plan_device, n0_FFT_in_device, n0_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, n1_FFT_in_device, n1_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, n2_FFT_in_device, n2_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, n3_FFT_in_device, n3_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv1x_FFT_in_device, nv1x_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv1y_FFT_in_device, nv1y_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv1z_FFT_in_device, nv1z_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv2x_FFT_in_device, nv2x_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv2y_FFT_in_device, nv2y_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv2z_FFT_in_device, nv2z_FFT_out_device, CUFFT_FORWARD);



	for (j=0; j<k_gas[0]; j++)
    {
    	sum_F_deri<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    	(na_device, nb_device, nc_device, diameter_HS_device, k_gas_device, FFT_in_device, n0_FFT_out_device, 
		n1_FFT_out_device, n2_FFT_out_device, n3_FFT_out_device, nv1x_FFT_out_device, nv1y_FFT_out_device, nv1z_FFT_out_device, 
		nv2x_FFT_out_device, nv2y_FFT_out_device, nv2z_FFT_out_device, w_n0_FFT_device, w_n3_FFT_device, w_nv2x_FFT_device, 
		w_nv2y_FFT_device, w_nv2z_FFT_device, j);

		if (cufftExecZ2Z(fft_plan_device, FFT_in_device, FFT_out_device, CUFFT_INVERSE) != CUFFT_SUCCESS)
		{
			printf("CUFFT error: Forward FFT of density FAILED at specie:\t%d\n", i);
		}

		cal_grad<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
		(na_device, nb_device, nc_device, chem_device, FFT_out_device, denskr_device, Vext_device, sqrt_dens_device, rho_bulk_device, gradient_device, j);
    }

    cudaMemcpy(gradient, gradient_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0], cudaMemcpyDeviceToHost);
}





double myvalgrad(double *gradient, double *sqrt_dens, long int system_size)
{
	cudaMemcpy(sqrt_dens_device, sqrt_dens, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0], cudaMemcpyHostToDevice);
	// 
	int i, j;



	// initialize the weighted density for hard sphere term
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n0_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n1_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n2_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], n3_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv1x_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv1y_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv1z_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv2x_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv2y_FFT_in_device);
	ini_device_value_complex<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0], nv2z_FFT_in_device);

	ini_device_value_double<<<(int)((na[0]*nb[0]*nc[0]*k_gas[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0]*k_gas[0], denskr_device);
	ini_device_value_double<<<(int)((na[0]*nb[0]*nc[0]*k_gas[0]-1)/running_block_size+1),running_block_size>>>
	(na[0]*nb[0]*nc[0]*k_gas[0], denski_device);

	// calculate the weighted density
    for (j=0; j<k_gas[0]; j++)
    {
    	real_split_into_complex_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    	(na[0]*nb[0]*nc[0], sqrt_dens_device, FFT_in_device, j);

    	if (cufftExecZ2Z(fft_plan_device, FFT_in_device, FFT_out_device, CUFFT_FORWARD) != CUFFT_SUCCESS)
		{
			printf("CUFFT error: Forward FFT of density FAILED at specie:\n");
		}

		add_weighted_density<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
		(na_device, nb_device, nc_device, diameter_HS_device, k_gas_device, FFT_out_device, n0_FFT_in_device, 
		n1_FFT_in_device, n2_FFT_in_device, n3_FFT_in_device, nv1x_FFT_in_device, nv1y_FFT_in_device, nv1z_FFT_in_device, 
		nv2x_FFT_in_device, nv2y_FFT_in_device, nv2z_FFT_in_device, w_n0_FFT_device, w_n3_FFT_device, w_nv2x_FFT_device, 
		w_nv2y_FFT_device, w_nv2z_FFT_device, uatt_device, denskr_device, denski_device, j);

    }

    cufftExecZ2Z(fft_plan_device, n0_FFT_in_device, n0_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, n1_FFT_in_device, n1_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, n2_FFT_in_device, n2_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, n3_FFT_in_device, n3_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv1x_FFT_in_device, nv1x_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv1y_FFT_in_device, nv1y_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv1z_FFT_in_device, nv1z_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv2x_FFT_in_device, nv2x_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv2y_FFT_in_device, nv2y_FFT_out_device, CUFFT_INVERSE);
	cufftExecZ2Z(fft_plan_device, nv2z_FFT_in_device, nv2z_FFT_out_device, CUFFT_INVERSE);

	for (j=0; j<k_gas[0]; j++)
    {
    	real_split_into_complex_separate_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    	(na[0]*nb[0]*nc[0], denskr_device, denski_device, FFT_in_device, j);

    	if (cufftExecZ2Z(fft_plan_device, FFT_in_device, FFT_out_device, CUFFT_INVERSE) != CUFFT_SUCCESS)
		{
			printf("CUFFT error: Backward FFT of density FAILED at specie:\n");
		}

		complext_split_into_real_separate_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
		(na[0]*nb[0]*nc[0], denskr_device, denski_device, FFT_out_device, j);
    }

    // calculate grand potential
    cal_grand_potential<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    (na_device, nb_device, nc_device, k_gas_device, Vext_device, chem_device, sqrt_dens_device, rho_bulk_device, n0_FFT_out_device, 
	n1_FFT_out_device, n2_FFT_out_device, n3_FFT_out_device, nv1x_FFT_out_device, nv1y_FFT_out_device, nv1z_FFT_out_device, 
	nv2x_FFT_out_device, nv2y_FFT_out_device, nv2z_FFT_out_device, denskr_device, grand_pontential_device);

    // calculate the derivative term of weighted density
    cal_F_deri<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    (na_device, nb_device, nc_device, n0_FFT_out_device, n1_FFT_out_device, n2_FFT_out_device, n3_FFT_out_device, 
	nv1x_FFT_out_device, nv1y_FFT_out_device, nv1z_FFT_out_device, nv2x_FFT_out_device, nv2y_FFT_out_device, nv2z_FFT_out_device,
	n0_FFT_in_device, n1_FFT_in_device, n2_FFT_in_device, n3_FFT_in_device, nv1x_FFT_in_device, nv1y_FFT_in_device, nv1z_FFT_in_device, 
	nv2x_FFT_in_device, nv2y_FFT_in_device, nv2z_FFT_in_device);

	//Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
    //Revision needed for the speed!!!!!!!!!
	cudaMemcpy(grand_pontential, grand_pontential_device, sizeof(double)*na[0]*nb[0]*nc[0], cudaMemcpyDeviceToHost);
	double result = 0;
	for (i=0; i<na[0]*nb[0]*nc[0]; i++)
	{
		result += grand_pontential[i];
	}
	result = result * (la_expanded[0]/(na[0])) * (lb_expanded[0]/(nb[0])) * (lc_expanded[0]/(nc[0]));

	cufftExecZ2Z(fft_plan_device, n0_FFT_in_device, n0_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, n1_FFT_in_device, n1_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, n2_FFT_in_device, n2_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, n3_FFT_in_device, n3_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv1x_FFT_in_device, nv1x_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv1y_FFT_in_device, nv1y_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv1z_FFT_in_device, nv1z_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv2x_FFT_in_device, nv2x_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv2y_FFT_in_device, nv2y_FFT_out_device, CUFFT_FORWARD);
	cufftExecZ2Z(fft_plan_device, nv2z_FFT_in_device, nv2z_FFT_out_device, CUFFT_FORWARD);
	
	for (j=0; j<k_gas[0]; j++)
    {
    	sum_F_deri<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
    	(na_device, nb_device, nc_device, diameter_HS_device, k_gas_device, FFT_in_device, n0_FFT_out_device, 
		n1_FFT_out_device, n2_FFT_out_device, n3_FFT_out_device, nv1x_FFT_out_device, nv1y_FFT_out_device, nv1z_FFT_out_device, 
		nv2x_FFT_out_device, nv2y_FFT_out_device, nv2z_FFT_out_device, w_n0_FFT_device, w_n3_FFT_device, w_nv2x_FFT_device, 
		w_nv2y_FFT_device, w_nv2z_FFT_device, j);

		if (cufftExecZ2Z(fft_plan_device, FFT_in_device, FFT_out_device, CUFFT_INVERSE) != CUFFT_SUCCESS)
		{
			printf("CUFFT error: Forward FFT of density FAILED at specie:\t%d\n", i);
		}

		cal_grad<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
		(na_device, nb_device, nc_device, chem_device, FFT_out_device, denskr_device, Vext_device, sqrt_dens_device, rho_bulk_device, gradient_device, j);
    }

    cudaMemcpy(gradient, gradient_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0], cudaMemcpyDeviceToHost);
    return result;
}























int main(int argc, char *argv[])
{
	clock_t t;
	t = clock();
	clock_t temp_t;
	double t_ini, t_w_d=0, t_f_deri=0;
	// clock_t t_ini, t_w_d, t_f_deri;

	// define variables used for file I/O
	FILE *fp1;
    int buffersize = 256;
    char str[buffersize];
    int int_a;
    // define variables assisst the main fxn
    int i, ii, iii, iiii;
    int j, jj;
    
    
    
    // define key variables read from the input file and independent of kind of gas
    int *Nmax_a, *Nmax_b, *Nmax_c;
    double *la, *lb, *lc, *dl;
    double *alpha, *beta, *gamma;
    double *temperature;
    int *signal_read_Vext, *signal_EOS;
    double *cutoff, *cutoff_device;
    int *N_atoms;
    // allocate memory for those variables
    Nmax_a = (int *) malloc(sizeof(int));
    Nmax_b = (int *) malloc(sizeof(int));
    Nmax_c = (int *) malloc(sizeof(int));
    la = (double *) malloc(sizeof(double));
    lb = (double *) malloc(sizeof(double));
    lc = (double *) malloc(sizeof(double));
    dl = (double *) malloc(sizeof(double));
    alpha = (double *) malloc(sizeof(double));
    beta = (double *) malloc(sizeof(double));
    gamma = (double *) malloc(sizeof(double));
    temperature = (double *) malloc(sizeof(double));
    signal_read_Vext = (int *) malloc(sizeof(int));
    signal_EOS = (int *) malloc(sizeof(int));
    N_atoms = (int *) malloc(sizeof(int));
    // allocate pinned memory
    cudaMallocHost(&k_gas, sizeof(int));
    cudaMallocHost(&cutoff, sizeof(double));
    // allocate GPU memory
    cudaMalloc((void **)&k_gas_device, sizeof(int)*1);
    cudaMalloc((void **)&cutoff_device, sizeof(double)*1);
    // read varaiables from input
    fp1 = fopen(argv[1], "r");
    fgets(str, buffersize, fp1);
    fscanf(fp1, "%d %d %d\n", &Nmax_a[0], &Nmax_b[0], &Nmax_c[0]);
    fgets(str, buffersize, fp1);
    fscanf(fp1, "%lf %lf %lf %lf\n", &la[0], &lb[0], &lc[0], &dl[0]);
    fgets(str, buffersize, fp1);
    fscanf(fp1, "%lf %lf %lf\n", &alpha[0], &beta[0], &gamma[0]);
    fgets(str, buffersize, fp1);
    fscanf(fp1, "%lf ", &temperature[0]);
    fgets(str, buffersize, fp1);
    fgets(str, buffersize, fp1);
    fscanf(fp1, "%d\n", &k_gas[0]);
    cudaMemcpy(k_gas_device, k_gas, sizeof(int), cudaMemcpyHostToDevice);
    // cudaMemcpy(cutoff_device, cutoff, sizeof(double), cudaMemcpyHostToDevice);
    // check_int<<<1,32>>>(1, k_gas_device);
    // check_double<<<1,32>>>(1, cutoff_device);
    

    
    // define variables dependent on the kind of gas
    int *index_gas;
    double *epsilon;
    double *epsilon_star, *epsilon_star_device;
    double *sigma, *sigma_device;
    double *rho_bulk;
    // allocate memory for the variables
    index_gas = (int *) malloc(sizeof(int)*k_gas[0]);
    epsilon = (double *) malloc(sizeof(double)*k_gas[0]);
    epsilon_star = (double *) malloc(sizeof(double)*k_gas[0]);
    sigma = (double *) malloc(sizeof(double)*k_gas[0]);
    rho_bulk = (double *) malloc(sizeof(double)*k_gas[0]);
    // allocate pinned memory
    cudaMallocHost(&epsilon_star, sizeof(double)*k_gas[0]);
    cudaMallocHost(&sigma, sizeof(double)*k_gas[0]);
    cudaMallocHost(&rho_bulk, sizeof(double)*k_gas[0]);
    // allocate GPU memory
    cudaMalloc((void **)&epsilon_star_device, sizeof(double)*k_gas[0]);
    cudaMalloc((void **)&sigma_device, sizeof(double)*k_gas[0]);
    cudaMalloc((void **)&rho_bulk_device, sizeof(double)*k_gas[0]);

    // read variables from input
    fgets(str, buffersize, fp1);
    for (i=0; i<k_gas[0]; i++)
    {
        fscanf(fp1, "%d %lf %lf %lf\n", &index_gas[i], &epsilon[i], &sigma[i], &rho_bulk[i]);
        epsilon_star[i] = 1.0*epsilon[i]/temperature[0];
    }
    cudaMemcpy(epsilon_star_device, epsilon_star, sizeof(double)*k_gas[0], cudaMemcpyHostToDevice);
    cudaMemcpy(sigma_device, sigma, sizeof(double)*k_gas[0], cudaMemcpyHostToDevice);
    cudaMemcpy(rho_bulk_device, rho_bulk, sizeof(double)*k_gas[0], cudaMemcpyHostToDevice);


    // check_double<<<1,32>>>(k_gas[0], epsilon_star_device);
    // check_double<<<1,32>>>(k_gas[0], sigma_device);
    // check_double<<<1,32>>>(k_gas[0], rho_bulk_device);
    
    
    
    // skip special case!!!!!!!!!!!!!
    // skip special case!!!!!!!!!!!!!
    // skip special case!!!!!!!!!!!!!
    fgets(str, buffersize, fp1);
    fgets(str, buffersize, fp1);
    fgets(str, buffersize, fp1);
    
    
    
    // read other variables
    fgets(str, buffersize, fp1);
    fscanf(fp1, "%d %d %lf", &signal_read_Vext[0], &signal_EOS[0], &cutoff[0]);
    cudaMemcpy(cutoff_device, cutoff, sizeof(double)*1, cudaMemcpyHostToDevice);
    // check_double<<<1,32>>>(1, cutoff_device);
    fgets(str, buffersize, fp1);
    
    
    
    // skip excess entropy scalling parameters!!!!!!!!!!!!!!!!!
    // skip excess entropy scalling parameters!!!!!!!!!!!!!!!!!
    // skip excess entropy scalling parameters!!!!!!!!!!!!!!!!!
    fgets(str, buffersize, fp1);
    fgets(str, buffersize, fp1);
    
    
    
    // read host materials
    fgets(str, buffersize, fp1);
    fscanf(fp1, "%d\n", &N_atoms[0]);
    
    // creat variables based on input files
    double *sigma_host, *epsilon_host;
    double *x_host, *y_host, *z_host;
    // allocate memory
    sigma_host = (double *) malloc(N_atoms[0]*sizeof(double));
    epsilon_host = (double *) malloc(N_atoms[0]*sizeof(double));
    x_host = (double *) malloc(N_atoms[0]*sizeof(double));
    y_host = (double *) malloc(N_atoms[0]*sizeof(double));
    z_host = (double *) malloc(N_atoms[0]*sizeof(double));
    
    // continue read info from input file
    fgets(str, buffersize, fp1);
    for (i=0; i<N_atoms[0]; i++)
    {
        fscanf(fp1,"%d %lf %lf\n", &int_a, &sigma_host[i], &epsilon_host[i]);
    }
    fgets(str, buffersize, fp1);
    fgets(str, buffersize, fp1);
    fgets(str, buffersize, fp1);
    for (i=0; i<N_atoms[0]; i++)
    {
        fscanf(fp1,"%d %lf %lf %lf\n", &int_a, &x_host[i], &y_host[i], &z_host[i]);
    }
    fclose(fp1);
    
    
    
    // check whether input parameters are correct
    // printf("1:\t%d\t%d\t%d\n", Nmax_a[0], Nmax_b[0], Nmax_c[0]);
    // printf("2:\t%lf\t%lf\t%lf\t%lf\n", la[0], lb[0], lc[0], dl[0]);
    // printf("3:\t%lf\t%lf\t%lf\n", alpha[0], beta[0], gamma[0]);
    // printf("4:\t%lf\n", temperature[0]);
    // printf("5:\t%d\n", k_gas[0]);
    // for (i=0; i<k_gas[0]; i++)
    // {
    //     printf("index of gas:\t%d\t%lf\t%lf\t%lf\n", index_gas[i], epsilon[i], sigma[i], rho_bulk[i]);
    // }
    // printf("6:\t%d\t%d\t%lf\n", signal_read_Vext[0], signal_EOS[0], cutoff[0]);
    // printf("7:\t%d\n", N_atoms[0]);
    // for (i=0; i<N_atoms[0]; i++)
    // {
    //     printf("%d %lf %lf\n", i+1, sigma_host[i], epsilon_host[i]);
    // }
    // for (i=0; i<N_atoms[0]; i++)
    // {
    //     printf("%d %lf %lf %lf\n", i+1, x_host[i], y_host[i], z_host[i]);
    // }


    


    // The following section is correct assume the lattice cell is cubic.
    // Since Yu's program is not incorrect for triclinic cell, we are only
    // making sure the corresponding C code work for the cubic case

    // define variables used for expanding original to satisfy the periodic boundary condition
    int *expand_times_a, *expand_times_b, *expand_times_c, *expand_times_total;
    double *la_expanded_device;
    double *lb_expanded_device;
    double *lc_expanded_device;
    
    int *N_atoms_expanded, *N_atoms_expanded_device;
    // allocate memory
    expand_times_a = (int *) malloc(sizeof(int));
    expand_times_b = (int *) malloc(sizeof(int));
    expand_times_c = (int *) malloc(sizeof(int));
    expand_times_total = (int *) malloc(sizeof(int));
    // allocate pinned memory
    cudaMallocHost(&la_expanded, sizeof(double));
    cudaMallocHost(&lb_expanded, sizeof(double));
    cudaMallocHost(&lc_expanded, sizeof(double));
    cudaMallocHost(&na, sizeof(int));
    cudaMallocHost(&nb, sizeof(int));
    cudaMallocHost(&nc, sizeof(int));
    cudaMallocHost(&N_atoms_expanded, sizeof(int));
    // allocate GPU memory
    cudaMalloc((void **)&la_expanded_device, sizeof(double)*1);
    cudaMalloc((void **)&lb_expanded_device, sizeof(double)*1);
    cudaMalloc((void **)&lc_expanded_device, sizeof(double)*1);
    cudaMalloc((void **)&na_device, sizeof(int)*1);
    cudaMalloc((void **)&nb_device, sizeof(int)*1);
    cudaMalloc((void **)&nc_device, sizeof(int)*1);
    cudaMalloc((void **)&N_atoms_expanded_device, sizeof(int)*1);

    // calculate expanding time in each axis
    expand_times_a[0] = (int) ((2*cutoff[0]/la[0]) + 1);
    expand_times_b[0] = (int) ((2*cutoff[0]/lb[0]) + 1);
    expand_times_c[0] = (int) ((2*cutoff[0]/lc[0]) + 1);
    expand_times_total[0] = (int) (expand_times_a[0]*expand_times_b[0]*expand_times_c[0]);
    N_atoms_expanded[0] = N_atoms[0] * expand_times_total[0];
    la_expanded[0] = la[0] * expand_times_a[0];
    lb_expanded[0] = lb[0] * expand_times_b[0];
    lc_expanded[0] = lc[0] * expand_times_c[0];
    na[0] = (int) (la_expanded[0]/dl[0] + 1);
    if (na[0] > Nmax_a[0])
    {
        na[0] = Nmax_a[0];
    }
    nb[0] = (int) (lb_expanded[0]/dl[0] + 1);
    if (nb[0] > Nmax_b[0])
    {
        nb[0] = Nmax_b[0];
    }
    nc[0] = (int) (lc_expanded[0]/dl[0] + 1);
    if (nc[0] > Nmax_c[0])
    {
        nc[0] = Nmax_c[0];
    }
    cudaMemcpy(na_device, na, sizeof(int)*1, cudaMemcpyHostToDevice);
    cudaMemcpy(nb_device, nb, sizeof(int)*1, cudaMemcpyHostToDevice);
    cudaMemcpy(nc_device, nc, sizeof(int)*1, cudaMemcpyHostToDevice);
    cudaMemcpy(la_expanded_device, la_expanded, sizeof(double)*1, cudaMemcpyHostToDevice);
    cudaMemcpy(lb_expanded_device, lb_expanded, sizeof(double)*1, cudaMemcpyHostToDevice);
    cudaMemcpy(lc_expanded_device, lc_expanded, sizeof(double)*1, cudaMemcpyHostToDevice);
    cudaMemcpy(N_atoms_expanded_device, N_atoms_expanded, sizeof(int)*1, cudaMemcpyHostToDevice);
    // check_int<<<1,32>>>(1, na_device);
    // check_int<<<1,32>>>(1, nb_device);
    // check_int<<<1,32>>>(1, nc_device);
    // check_double<<<1,32>>>(1, la_expanded_device);
    // check_double<<<1,32>>>(1, lb_expanded_device);
    // check_double<<<1,32>>>(1, lc_expanded_device);
    // check_int<<<1,32>>>(1, N_atoms_expanded_device);

    // define variables dependent on the expand times in each axis
    double *sigma_host_expanded, *sigma_host_expanded_device;
    double *epsilon_host_star_expanded, *epsilon_host_star_expanded_device;
    double *x_host_expanded, *x_host_expanded_device;
    double *y_host_expanded, *y_host_expanded_device;
    double *z_host_expanded, *z_host_expanded_device;
    // allocate pinned memroy
    cudaMallocHost(&sigma_host_expanded, sizeof(double)*expand_times_total[0]*N_atoms[0]);
    cudaMallocHost(&epsilon_host_star_expanded, sizeof(double)*expand_times_total[0]*N_atoms[0]);
    cudaMallocHost(&x_host_expanded, sizeof(double)*expand_times_total[0]*N_atoms[0]);
    cudaMallocHost(&y_host_expanded, sizeof(double)*expand_times_total[0]*N_atoms[0]);
    cudaMallocHost(&z_host_expanded, sizeof(double)*expand_times_total[0]*N_atoms[0]);
    // allocate GPU memory
    cudaMalloc((void **)&sigma_host_expanded_device, sizeof(double)*expand_times_total[0]*N_atoms[0]);
    cudaMalloc((void **)&epsilon_host_star_expanded_device, sizeof(double)*expand_times_total[0]*N_atoms[0]);
    cudaMalloc((void **)&x_host_expanded_device, sizeof(double)*expand_times_total[0]*N_atoms[0]);
    cudaMalloc((void **)&y_host_expanded_device, sizeof(double)*expand_times_total[0]*N_atoms[0]);
    cudaMalloc((void **)&z_host_expanded_device, sizeof(double)*expand_times_total[0]*N_atoms[0]);


    
    j = 0;
    for (i=0; i<N_atoms[0]; i++)
    {
        for (ii=0; ii<expand_times_a[0]; ii++)
        {
            for (iii=0; iii<expand_times_b[0]; iii++)
            {
                for (iiii=0; iiii<expand_times_c[0]; iiii++)
                {
                    epsilon_host_star_expanded[j] = 1.0*epsilon_host[i]/temperature[0];
                    sigma_host_expanded[j] = sigma_host[i];
                    x_host_expanded[j] = x_host[i] + ii*la[0];
                    y_host_expanded[j] = y_host[i] + iii*lb[0];
                    z_host_expanded[j] = z_host[i] + iiii*lc[0];
                    j++;
                }
            }
        }
    }
    cudaMemcpy(sigma_host_expanded_device, sigma_host_expanded, sizeof(double)*expand_times_total[0]*N_atoms[0], cudaMemcpyHostToDevice);
    cudaMemcpy(epsilon_host_star_expanded_device, epsilon_host_star_expanded, sizeof(double)*expand_times_total[0]*N_atoms[0], cudaMemcpyHostToDevice);
    cudaMemcpy(x_host_expanded_device, x_host_expanded, sizeof(double)*expand_times_total[0]*N_atoms[0], cudaMemcpyHostToDevice);
    cudaMemcpy(y_host_expanded_device, y_host_expanded, sizeof(double)*expand_times_total[0]*N_atoms[0], cudaMemcpyHostToDevice);
    cudaMemcpy(z_host_expanded_device, z_host_expanded, sizeof(double)*expand_times_total[0]*N_atoms[0], cudaMemcpyHostToDevice);
    // check_double<<<1,32>>>(expand_times_total[0]*N_atoms[0], sigma_host_expanded_device);
    // check_double<<<1,32>>>(expand_times_total[0]*N_atoms[0], epsilon_host_star_expanded);
    // check_double<<<1,32>>>(expand_times_total[0]*N_atoms[0], x_host_expanded_device);
    // check_double<<<1,32>>>(expand_times_total[0]*N_atoms[0], y_host_expanded_device);
    // check_double<<<1,32>>>(expand_times_total[0]*N_atoms[0], z_host_expanded_device);
    if (j!=N_atoms_expanded[0])
    {
        printf("error in expanding the original cell!!!\n");
    }
    // printf("expand_time:\t%d\t%d\t%d\n", expand_times_a[0], expand_times_b[0], expand_times_c[0]);
    // printf("li_expanded:\t%lf\t%lf\t%lf\n", la_expanded[0], lb_expanded[0], lc_expanded[0]);
    // printf("li_expanded:\t%d\t%d\t%d\n", na[0], nb[0], nc[0]);

    // free the memory occupied by the original cell
    free(epsilon_host);
    free(sigma_host);
    free(x_host);
    free(y_host);
    free(z_host);
    free(la);
    free(lb);
    free(lc);
    free(N_atoms);





    // calculate corresponding hard sphere diameter according to LJ parameters
    double *Tstar;
    double *diameter_HS;
    Tstar = (double *) malloc(sizeof(double)*k_gas[0]);
    cudaMallocHost(&diameter_HS, sizeof(double)*k_gas[0]);
    cudaMalloc((void **)&diameter_HS_device, sizeof(double)*k_gas[0]);

    for (i=0; i<k_gas[0]; i++)
    {
        Tstar[i] = temperature[0]/epsilon[i];
        if ((Tstar[i]>0) && (Tstar[i]<15))
        {
            // printf("Reduced temperature is within good range.\n");
        }
        else
        {
            // printf("WARNING: reduced temperature is out of reasonable range!!!!!!\n");
            // printf("Index of gas: %d\tReduced Temperature: %lf\n", i, Tstar[i]);
        }
        diameter_HS[i] = (1+0.2977*Tstar[i]) / (1+0.33163*Tstar[i]+0.0010477*pow(Tstar[i],2))*sigma[i];
        // printf("%lf\t%lf\n", Tstar[i], diameter_HS[i]);
    }
    cudaMemcpy(diameter_HS_device, diameter_HS, sizeof(double)*k_gas[0], cudaMemcpyHostToDevice);
    // check_double<<<1,32>>>(k_gas[0], diameter_HS_device);





    // variables for multifunction used in the initiliazatiion
    int *nijk_gas;
    nijk_gas = (int *) malloc(sizeof(int));


    
    // allocate GPU memory
    cudaMalloc((void **)&w_n0_FFT_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);
    cudaMalloc((void **)&w_n3_FFT_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);
    cudaMalloc((void **)&w_nv2x_FFT_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);
    cudaMalloc((void **)&w_nv2y_FFT_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);
    cudaMalloc((void **)&w_nv2z_FFT_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);
    cudaMalloc((void **)&uatt_device, sizeof(double)*na[0]*nb[0]*nc[0]*(k_gas[0] + (int) (0.5*(k_gas[0]*(k_gas[0]-1)))));



    // This sectin works the same as the interpolation way Yu Liu has for the weight density but with GPU
    // However, this following section follows a more standard routine
    cal_weight_density<<<(int)((na[0]*nb[0]*nc[0]*k_gas[0]-1)/running_block_size+1),running_block_size>>>
    (na_device, nb_device, nc_device, la_expanded_device, lb_expanded_device, lc_expanded_device, k_gas_device, 
    diameter_HS_device, w_n0_FFT_device, w_n3_FFT_device, w_nv2x_FFT_device, w_nv2y_FFT_device, w_nv2z_FFT_device);


	
    cudaMalloc((void **)&FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    
    if (cufftPlan3d(&fft_plan_device, na[0], nb[0], nc[0], CUFFT_Z2Z) != CUFFT_SUCCESS)
	{
		printf("CUFFT error: Plan creation failed\n");
	}




    
    // putting attraction terms on GPU in order to speed up GPU program
    // and also to avoid the call of fftw3
    for (j=0; j<k_gas[0]; j++)
    {
    	for (jj=0; jj<=j; jj++)
    	{
    		cal_lj_pre_freq_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
		    (na_device, nb_device, nc_device, la_expanded_device, lb_expanded_device, lc_expanded_device, sigma_device, 
		    epsilon_star_device, cutoff_device, j, jj, FFT_in_device);

		    cufftExecZ2Z(fft_plan_device, FFT_in_device, FFT_out_device, CUFFT_FORWARD);

		    cal_lj_freq_store_batch<<<(int)((na[0]*nb[0]*nc[0]-1)/running_block_size+1),running_block_size>>>
		    (na_device, nb_device, nc_device, la_expanded_device, lb_expanded_device, lc_expanded_device, j, jj, 
		    FFT_out_device, uatt_device);


    	}
    }





    // calculate the residual chemical potential and full chemical potential
    double *chem;
    double *bulk_pressure;
    bulk_pressure = (double *) malloc(sizeof(double));
    // allocate pinned memory
    cudaMallocHost(&chem, sizeof(double)*k_gas[0]);
    // allocate GPU memory
    cudaMalloc((void **)&chem_device, sizeof(double)*k_gas[0]);

    // calculate the system pressure and reduced chemical potential
    new_MBWR(k_gas, epsilon, sigma, temperature, rho_bulk, chem, bulk_pressure);
    
    cudaMemcpy(chem_device, chem, sizeof(double)*k_gas[0], cudaMemcpyHostToDevice);
    // printf("chemical potential\n");
    // check_double<<<1,32>>>(k_gas[0], chem_device);





    // define variables
    double *sqrt_dens;
    // allocate pinned memroy
    cudaMallocHost(&sqrt_dens, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);
    // allocate GPU memory
    cudaMalloc((void **)&sqrt_dens_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);


    
    cudaMalloc((void **)&Vext_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);



    cal_Vext<<<(int)((na[0]*nb[0]*nc[0]*k_gas[0]-1)/running_block_size+1),running_block_size>>>
    (na_device, nb_device, nc_device, la_expanded_device, lb_expanded_device, lc_expanded_device, 
	N_atoms_expanded_device, epsilon_host_star_expanded_device, sigma_host_expanded_device,
	x_host_expanded_device, y_host_expanded_device, z_host_expanded_device,
	k_gas_device, epsilon_star_device, sigma_device, rho_bulk_device, 
	diameter_HS_device, sqrt_dens_device, cutoff_device, Vext_device);



    cudaMemcpy(sqrt_dens, sqrt_dens_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0], cudaMemcpyDeviceToHost);


	
    cudaMalloc((void **)&n0_FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&n1_FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&n2_FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&n3_FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv1x_FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv1y_FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv1z_FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv2x_FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv2y_FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv2z_FFT_in_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);

    
    cudaMalloc((void **)&n0_FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&n1_FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&n2_FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&n3_FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv1x_FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv1y_FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv1z_FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv2x_FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv2y_FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);
    cudaMalloc((void **)&nv2z_FFT_out_device, sizeof(cufftDoubleComplex)*na[0]*nb[0]*nc[0]);



    
    cudaMalloc((void **)&denskr_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);
    cudaMalloc((void **)&denski_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);



    double *tot_density;
    tot_density = (double *) malloc(sizeof(double)*k_gas[0]);



    // allocate pinned memory
    cudaMallocHost(&grand_pontential, sizeof(double)*na[0]*nb[0]*nc[0]);
    // allocate GPU memory
    cudaMalloc((void **)&grand_pontential_device, sizeof(double)*na[0]*nb[0]*nc[0]);

    int system_size = na[0]*nb[0]*nc[0]*k_gas[0];
    
    double *gradient;
    // allocate pinned memory
    cudaMallocHost(&gradient, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);
    // allocate GPU memory
    cudaMalloc((void **)&gradient_device, sizeof(double)*na[0]*nb[0]*nc[0]*k_gas[0]);




    // myvalue(sqrt_dens, system_size);
    // check_double<<<1,32>>>(k_gas[0], rho_bulk_device);
    // mygrad(gradient, sqrt_dens, system_size);
    // myvalgrad(gradient, sqrt_dens, system_size);
    
    // cudaDeviceSynchronize();
    // return 0;


    int output_signal;
    output_signal = cg_descent(sqrt_dens, system_size, NULL, NULL, 1.e-3, myvalue, mygrad, myvalgrad, NULL) ;
    cudaDeviceSynchronize();
    // printf("%s\t", argv[1]);
    printf("%d\t", output_signal);
    // return 0;

    t = clock() - t;
    double time_used;
    time_used = ((double)t)/CLOCKS_PER_SEC;
    printf("%lf\t", bulk_pressure[0]);
    for (j=0; j<k_gas[0]; j++)
	{
		tot_density[j] = 0;
		for (i=0; i<na[0]*nb[0]*nc[0]; i++)
		{
			nijk_gas[0] = i + j*na[0]*nb[0]*nc[0];
			// dens[nijk_gas[0]] = dens_new[nijk_gas[0]]*mix_F[0] + (1-mix_F[0])*dens[nijk_gas[0]];
			tot_density[j] = tot_density[j] + sqrt_dens[nijk_gas[0]]*sqrt_dens[nijk_gas[0]];
		}
		tot_density[j] = 1.0*tot_density[j]*1.0e4/6.023/na[0]/nb[0]/nc[0];
		// tot_density[j] = 1.0*tot_density[j]/na[0]/nb[0]/nc[0];
		// printf("adsorption amount: %lf (mol/L)\n", tot_density[j]);
		printf("%.5e\t", tot_density[j]);
	}
	printf("\n");
}
