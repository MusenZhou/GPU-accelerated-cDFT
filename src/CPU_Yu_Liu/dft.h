#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include "cg_user.h"

#define PI 3.1415926535898
#define MIN 0.000000000001
#define PRE 0.0000000001
#define STEP 10000
#define ZOOM 3.
#define DE 0.001//数值微分比例
class LJ
{
public:
	int len,kind;
	double gama;
	double **dia,**eff;//直径，e/k
	double tem;//温度
	double x[32],a[8],b[6];
	double aa[5][5];
	bool initial();
	//double balj(double);
	//double debalj(double);
	double out();
	void debaex(double*,double*);
	double baex(double*);
	double buex(double*);
	void r_chem(double*,double*);
	void chem_r(double*,double*);
	double press(double*);
	double s_per_n(double*);//熵（对粒子）
	double f_per_n(double*);//自由能（对粒子）
	double sex_nk(double*);//Sex/(NkB)
	//double chempot(double*);
	LJ();
	LJ(double,double,double);
	LJ(int,double*,double*,double);
	~LJ();
};

class FLJ
{
public:
  FLJ(double,double);//initialize (beta*epsilon,sigma,rou)
  double run(double);//calculate c function cattr
private:
  double be,sig;//beta*epsilon, sigma,rou
  double d,y,besd1,besd2,bey[2],lam[2];
};

class Vector
{
public:
  double x;
  double y;
  double z;
  void initial(double,double,double);
  double rabs();
};

class Weight_K
{
public:
  void initial(double);
  double w0(double);
  double w1(double);
  double w2(double);
  double w3(double);
  Vector wv1(Vector);
  Vector wv2(Vector);  

  double dia;
  double dw;
  int nw;
  double fw0[5000],fw3[5000],fwcor[5000];
};

class DFT
{
public:
  DFT();
  int run();  
  ~DFT();
  
  void fdf_dd(double*,double*,long);
  double fdf(double*,long); 
  double fdf_t(double*,double*,long);
  long tot_grids();
  long tot_sites();
  double time_cost();
  
  void cal_ck();
  void final();
  void get_sq_dens(double*);
  
  double *m_x;
  
  double sex();
private:
  //void cal_ck();
  void cal_dens();
  //void final();
  double dis(double,double,double,double,double,double);
  int norx(int);
  int nory(int);
  int norz(int);
  Vector tr_cu_k(Vector);//translate the corrdinates in tr axis to cubic axis
  Vector tr_cu_k(double,double,double);
  Vector tr_cu_r(Vector);//translate the corrdinates in tr axis to cubic axis
  Vector tr_cu_r(double,double,double);
  double dis_tri_r(double,double,double);
  double dis_tri_k(double,double,double);
  long n_mx(int,int,int,int);
  //void fixpoint();
  //void fixpoint2();
  //void fixpoint3();
  
  int nx,ny,nz,k_gas;
  long nxyz;
  double dx,dy,dz,lx,ly,lz,dkx,dky,dkz,lxu,lyu,lzu,dv;
  double angle_a,angle_b,angle_c;//three angles of the unit cell
  double cosa2,cosb2,cosc2;
  double *be,*sig,*dbh,*rb,*bmhs,tem;//epsilon,sigma,BH diameter,bulk density,beta*miu_hs
  double **bemix,**sigmix;
  int *name_gas;
  double ****vext;
  int *id,*aid;
  double *ax,*ay,*az;
  double *pot,*dia;
  //double ***dens;
  //double ***rou;
  double kapa,detarou;
  //double ***lamd;
  double *****uattk;
  double *fatt,*rdfatt,*ds_rou;
  double ****tw0,****tw3,****twcor,****tw2x,****tw2y,****tw2z;
  double ****denskr,****denski;
  //double ****tw1,****tw2,****tw1x,****tw1y,****tw1z;
  double dfa;
  double rcut;
  double torr;
  int eos;//0 for FMSA, 1 for MWBR
  int deos;//0 for Rosenfeld, 1 for Silva (Chemical Engineering Science, 79, 153-162)
  int nfa;
  int step;
  int cut_l;
  int pdx,pdy,pdz,pdxyz;
  double cutoff;
  double srav,srmax,*rouav,rou_tot,*max_dens;//,srmax_old,srav_old,rouav_old;
  bool conv;
  double pbulk;
  double mmof,matom;
  double fex;
  double knudsen;//knudsen diffusion coeffcient
  double p_mfv;//mixing parameter in entropy scaling
  double mass;//mass of gas
  Weight_K *wfun;
  Vector a1,a2,a3,b1,b2,b3;
  double a123;
  int ercode;
  clock_t start,finish;
  ////////////test////////////
  //double ***newrou;
  ////////////////////////////
  ////////// fix point ///////////
  //int nfix,*xfix,*yfix,*zfix;
  double fix_cut;
  ////////////////////////////////

  fftw_complex *n0,*n3,*nv2x,*nv2y,*nv2z,*cor,*n1,*n2,*nv1x,*nv1y,*nv1z;
  fftw_complex *fft_in,*fft_out;
  fftw_plan fft_for,fft_bak,fft_n0_for,fft_n0_bak,fft_n3_for,fft_n3_bak,fft_cor_for,fft_cor_bak;
  fftw_plan fft_n1_for,fft_n2_for,fft_n1_bak,fft_n2_bak,fft_nv1x_for,fft_nv1x_bak,fft_nv1y_for,fft_nv1y_bak,fft_nv1z_for,fft_nv1z_bak;
  fftw_plan fft_nv2x_for,fft_nv2x_bak,fft_nv2y_for,fft_nv2y_bak,fft_nv2z_for,fft_nv2z_bak;
};

extern DFT *MINIMIZATION;

double myvalue(double*,long);

void mygrad(double*, double*, long);

double myvalgrad(double*, double*, long num);