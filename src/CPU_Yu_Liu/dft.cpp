#include "dft.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int min(int i,int j)
{
	if(i>j)return j;
	else return i;
}

int max(int i,int j)
{
	if(i>j)return i;
	else return j;
}

double min(double i,double j)
{
	if(i<j)return i;
	else return j;
}

double max(double i,double j)
{
	if(i>j)return i;
	else return j;
}

double juedui(double x)
{
	if(x<0)return -x;
	return x;
}

double dmod(double x,double y) //x%y
{
  int i;
  i=floor(x/y);
  return x-i*y;
}




double lj_to_hs(double tstar)
{
  double ds;
  ds=tstar-0.05536*tstar*tstar+0.0007278*pow(tstar,4);
  ds=ds/1.1287;
  ds=pow(sqrt(1.0+ds)+1.0,-1.0/6.0);
  ds*=pow(2.0,1.0/6.0);
  return ds;
}



LJ::LJ()
{
	int i,j;
	//dia=1.;
	//be=0.7;
	gama=3.;
	len=1;
	ifstream ip("input.dat",ios::in);
	ip>>kind>>tem;
	dia=new double*[kind];
	eff=new double*[kind];
	for(i=0;i<kind;i++)
	{
		dia[i]=new double[kind];
		eff[i]=new double[kind];
	}
	
	for(i=0;i<kind;i++)
	{
		ip>>dia[i][i]>>eff[i][i];
		eff[i][i]*=1./tem;
	}
	for(i=0;i<kind;i++)
	{
		for(j=0;j<i;j++)
		{
			dia[i][j]=(dia[i][i]+dia[j][j])*0.5;
			eff[i][j]=sqrt(eff[i][i]*eff[j][j]);
			dia[j][i]=dia[i][j];
			eff[j][i]=eff[i][j];
		}
	}
	ip.close();
	/////////////////x////////////////////
	x[0]=0.8623085097507421;
	x[1]=2.976218765822098;
	x[2]=-8.402230115796038;
	x[3]=0.1054136629203555;
	x[4]=-0.8564583828174598;
	x[5]=1.582759470107601;
	x[6]=0.7639421948305453;
	x[7]=1.753173414312048;
	x[8]=2.798291772190376e+03;
	x[9]=-4.8394220260857657e-02;
	x[10]=0.9963265197721935;
	x[11]=-3.698000291272493e+01;
	x[12]=2.084012299434647e+01;
	x[13]=8.305402124717285e+01;
	x[14]=-9.574799715203068e+02;
	x[15]=-1.477746229234994e+02;
	x[16]=6.398607852471505e+01;
	x[17]=1.603993673294834e+01;
	x[18]=6.805916615864377e+01;
	x[19]=-2.791293578795945e+03;
	x[20]=-6.245128304568454;
	x[21]=-8.116836104958410e+03;
	x[22]=1.488735559561229e+01;
	x[23]=-1.059346754655084e+04;
	x[24]=-1.131607632802822e+02;
	x[25]=-8.867771540418822e+03;
	x[26]=-3.986982844450543e+01;
	x[27]=-4.689270299917261e+03;
	x[28]=2.593535277438717e+02;
	x[29]=-2.694523589434903e+03;
	x[30]=-7.218487631550215e+02;
	x[31]=1.721802063863269e+02;
    //////////////////////////////////////
	////////////////a/////////////////////
	/*
	a[0]=x[0]/be+x[1]*sqrt(1./be)+x[2]+x[3]*be+x[4]*be*be;
	a[1]=x[5]/be+x[6]+x[7]*be+x[8]*be*be;
	a[2]=x[9]/be+x[10]+x[11]*be;
	a[3]=x[12];
	a[4]=x[13]*be+x[14]*be*be;
	a[5]=x[15]*be;
	a[6]=x[16]*be+x[17]*be*be;
	a[7]=x[18]*be*be;
	//////////////////////////////////////
	/////////////////b////////////////////
	b[0]=x[19]*be*be+x[20]*be*be*be;
	b[1]=x[21]*be*be+x[22]*be*be*be*be;
	b[2]=x[23]*be*be+x[24]*be*be*be;
	b[3]=x[25]*be*be+x[26]*be*be*be*be;
	b[4]=x[27]*be*be+x[28]*be*be*be;
	b[5]=x[29]*be*be+x[30]*be*be*be+x[31]*be*be*be*be;
	*//////////////////////////////////////
	////////////////aa////////////////////
	aa[0][0]=0.49304346593882;
	aa[0][1]=2.1528349894745;
	aa[0][2]=-15.955682329017;
	aa[0][3]=24.035999666294;
	aa[0][4]=-8.6437958513990;
	aa[1][0]=-0.47031983115362;
	aa[1][1]=1.1471647487376;
	aa[1][2]=37.889828024211;
	aa[1][3]=-84.667121491179;
	aa[1][4]=39.643914108411;
	aa[2][0]=5.0325486243620;
	aa[2][1]=-25.915399226419;
	aa[2][2]=-18.862251310090; 
	aa[2][3]=107.63707381726;
	aa[2][4]=-66.602649735720;
	aa[3][0]=-7.3633150434385;
	aa[3][1]=51.553565337453;
	aa[3][2]=-40.519369256098;
	aa[3][3]=-38.796692647218;
	aa[3][4]=44.605139198318;
	aa[4][0]=2.9043607296043;
	aa[4][1]=-24.418812869291;
	aa[4][2]=31.500186765040;
	aa[4][3]=-5.3368920371407;
	aa[4][4]=-9.5183440180133;
	//////////////////////////////////////
	cout<<"initial complete"<<endl;
	return ;
}

LJ::LJ(double sig, double pot,double tp)//pure component
{
	int i,j;
	//dia=1.;
	//be=0.7;
	gama=3.;
	len=1;
	kind=1;
	tem=tp;
	dia=new double*[kind];
	eff=new double*[kind];
	for(i=0;i<kind;i++)
	{
		dia[i]=new double[kind];
		eff[i]=new double[kind];
	}
	
	for(i=0;i<kind;i++)
	{
		dia[i][i]=sig;
		eff[i][i]=pot;
	}
	for(i=0;i<kind;i++)
	{
		for(j=0;j<i;j++)
		{
			dia[i][j]=(dia[i][i]+dia[j][j])*0.5;
			eff[i][j]=sqrt(eff[i][i]*eff[j][j]);
			dia[j][i]=dia[i][j];
			eff[j][i]=eff[i][j];
		}
	}
	/////////////////x////////////////////
	x[0]=0.8623085097507421;
	x[1]=2.976218765822098;
	x[2]=-8.402230115796038;
	x[3]=0.1054136629203555;
	x[4]=-0.8564583828174598;
	x[5]=1.582759470107601;
	x[6]=0.7639421948305453;
	x[7]=1.753173414312048;
	x[8]=2.798291772190376e+03;
	x[9]=-4.8394220260857657e-02;
	x[10]=0.9963265197721935;
	x[11]=-3.698000291272493e+01;
	x[12]=2.084012299434647e+01;
	x[13]=8.305402124717285e+01;
	x[14]=-9.574799715203068e+02;
	x[15]=-1.477746229234994e+02;
	x[16]=6.398607852471505e+01;
	x[17]=1.603993673294834e+01;
	x[18]=6.805916615864377e+01;
	x[19]=-2.791293578795945e+03;
	x[20]=-6.245128304568454;
	x[21]=-8.116836104958410e+03;
	x[22]=1.488735559561229e+01;
	x[23]=-1.059346754655084e+04;
	x[24]=-1.131607632802822e+02;
	x[25]=-8.867771540418822e+03;
	x[26]=-3.986982844450543e+01;
	x[27]=-4.689270299917261e+03;
	x[28]=2.593535277438717e+02;
	x[29]=-2.694523589434903e+03;
	x[30]=-7.218487631550215e+02;
	x[31]=1.721802063863269e+02;
    //////////////////////////////////////
	////////////////a/////////////////////
	/*
	a[0]=x[0]/be+x[1]*sqrt(1./be)+x[2]+x[3]*be+x[4]*be*be;
	a[1]=x[5]/be+x[6]+x[7]*be+x[8]*be*be;
	a[2]=x[9]/be+x[10]+x[11]*be;
	a[3]=x[12];
	a[4]=x[13]*be+x[14]*be*be;
	a[5]=x[15]*be;
	a[6]=x[16]*be+x[17]*be*be;
	a[7]=x[18]*be*be;
	//////////////////////////////////////
	/////////////////b////////////////////
	b[0]=x[19]*be*be+x[20]*be*be*be;
	b[1]=x[21]*be*be+x[22]*be*be*be*be;
	b[2]=x[23]*be*be+x[24]*be*be*be;
	b[3]=x[25]*be*be+x[26]*be*be*be*be;
	b[4]=x[27]*be*be+x[28]*be*be*be;
	b[5]=x[29]*be*be+x[30]*be*be*be+x[31]*be*be*be*be;
	*//////////////////////////////////////
	////////////////aa////////////////////
	aa[0][0]=0.49304346593882;
	aa[0][1]=2.1528349894745;
	aa[0][2]=-15.955682329017;
	aa[0][3]=24.035999666294;
	aa[0][4]=-8.6437958513990;
	aa[1][0]=-0.47031983115362;
	aa[1][1]=1.1471647487376;
	aa[1][2]=37.889828024211;
	aa[1][3]=-84.667121491179;
	aa[1][4]=39.643914108411;
	aa[2][0]=5.0325486243620;
	aa[2][1]=-25.915399226419;
	aa[2][2]=-18.862251310090; 
	aa[2][3]=107.63707381726;
	aa[2][4]=-66.602649735720;
	aa[3][0]=-7.3633150434385;
	aa[3][1]=51.553565337453;
	aa[3][2]=-40.519369256098;
	aa[3][3]=-38.796692647218;
	aa[3][4]=44.605139198318;
	aa[4][0]=2.9043607296043;
	aa[4][1]=-24.418812869291;
	aa[4][2]=31.500186765040;
	aa[4][3]=-5.3368920371407;
	aa[4][4]=-9.5183440180133;
	//////////////////////////////////////

	return ;
}

LJ::LJ(int k_in,double *dia_in,double *eff_in,double tp)
{
	int i,j;
	//dia=1.;
	//be=0.7;
	gama=3.;
	len=1;
	//ifstream ip("input.dat",ios::in);
	//ip>>kind>>tem;
	kind=k_in;
	tem=tp;
	dia=new double*[kind];
	eff=new double*[kind];
	for(i=0;i<kind;i++)
	{
		dia[i]=new double[kind];
		eff[i]=new double[kind];
	}
	
	for(i=0;i<kind;i++)
	{
		dia[i][i]=dia_in[i];
		eff[i][i]=eff_in[i];
		//eff[i][i]*=1./tem;
	}
	for(i=0;i<kind;i++)
	{
		for(j=0;j<i;j++)
		{
			dia[i][j]=(dia[i][i]+dia[j][j])*0.5;
			eff[i][j]=sqrt(eff[i][i]*eff[j][j]);
			dia[j][i]=dia[i][j];
			eff[j][i]=eff[i][j];
		}
	}
	//ip.close();
	/////////////////x////////////////////
	x[0]=0.8623085097507421;
	x[1]=2.976218765822098;
	x[2]=-8.402230115796038;
	x[3]=0.1054136629203555;
	x[4]=-0.8564583828174598;
	x[5]=1.582759470107601;
	x[6]=0.7639421948305453;
	x[7]=1.753173414312048;
	x[8]=2.798291772190376e+03;
	x[9]=-4.8394220260857657e-02;
	x[10]=0.9963265197721935;
	x[11]=-3.698000291272493e+01;
	x[12]=2.084012299434647e+01;
	x[13]=8.305402124717285e+01;
	x[14]=-9.574799715203068e+02;
	x[15]=-1.477746229234994e+02;
	x[16]=6.398607852471505e+01;
	x[17]=1.603993673294834e+01;
	x[18]=6.805916615864377e+01;
	x[19]=-2.791293578795945e+03;
	x[20]=-6.245128304568454;
	x[21]=-8.116836104958410e+03;
	x[22]=1.488735559561229e+01;
	x[23]=-1.059346754655084e+04;
	x[24]=-1.131607632802822e+02;
	x[25]=-8.867771540418822e+03;
	x[26]=-3.986982844450543e+01;
	x[27]=-4.689270299917261e+03;
	x[28]=2.593535277438717e+02;
	x[29]=-2.694523589434903e+03;
	x[30]=-7.218487631550215e+02;
	x[31]=1.721802063863269e+02;
    //////////////////////////////////////
	////////////////a/////////////////////
	/*
	a[0]=x[0]/be+x[1]*sqrt(1./be)+x[2]+x[3]*be+x[4]*be*be;
	a[1]=x[5]/be+x[6]+x[7]*be+x[8]*be*be;
	a[2]=x[9]/be+x[10]+x[11]*be;
	a[3]=x[12];
	a[4]=x[13]*be+x[14]*be*be;
	a[5]=x[15]*be;
	a[6]=x[16]*be+x[17]*be*be;
	a[7]=x[18]*be*be;
	//////////////////////////////////////
	/////////////////b////////////////////
	b[0]=x[19]*be*be+x[20]*be*be*be;
	b[1]=x[21]*be*be+x[22]*be*be*be*be;
	b[2]=x[23]*be*be+x[24]*be*be*be;
	b[3]=x[25]*be*be+x[26]*be*be*be*be;
	b[4]=x[27]*be*be+x[28]*be*be*be;
	b[5]=x[29]*be*be+x[30]*be*be*be+x[31]*be*be*be*be;
	*//////////////////////////////////////
	////////////////aa////////////////////
	aa[0][0]=0.49304346593882;
	aa[0][1]=2.1528349894745;
	aa[0][2]=-15.955682329017;
	aa[0][3]=24.035999666294;
	aa[0][4]=-8.6437958513990;
	aa[1][0]=-0.47031983115362;
	aa[1][1]=1.1471647487376;
	aa[1][2]=37.889828024211;
	aa[1][3]=-84.667121491179;
	aa[1][4]=39.643914108411;
	aa[2][0]=5.0325486243620;
	aa[2][1]=-25.915399226419;
	aa[2][2]=-18.862251310090; 
	aa[2][3]=107.63707381726;
	aa[2][4]=-66.602649735720;
	aa[3][0]=-7.3633150434385;
	aa[3][1]=51.553565337453;
	aa[3][2]=-40.519369256098;
	aa[3][3]=-38.796692647218;
	aa[3][4]=44.605139198318;
	aa[4][0]=2.9043607296043;
	aa[4][1]=-24.418812869291;
	aa[4][2]=31.500186765040;
	aa[4][3]=-5.3368920371407;
	aa[4][4]=-9.5183440180133;
	//////////////////////////////////////
	cout<<"initial complete"<<endl;
	return ;
}


LJ::~LJ()
{
	int i;
	for(i=0;i<kind;i++)
	{
		delete []dia[i];
		delete []eff[i];
	}
	delete []dia;
	delete []eff;
}
/*
double LJ::balj(double rou)
{
	int i,j;
	double roux,ba,f,ylj;
	double g[6];
	/////////////////g////////////////////
	roux=rou*dia*dia*dia;
	f=exp(-gama*roux*roux);
	g[0]=(1.-f)/(2*gama);
	for(i=1;i<6;i++)
	{
		g[i]=i/gama*g[i-1]-f/(2*gama)*pow(roux,2*i);
	}
	//////////////////////////////////////
	for(i=0,ba=0.;i<8;i++)
	{
		ba+=a[i]/(i+1)*pow(roux,i+1);
	}
	for(i=0;i<6;i++)
	{
		ba+=b[i]*g[i];
	}
	///////////////ylj////////////////////
	for(i=0,ylj=1.;i<5;i++)
	{
		for(j=0;j<5;j++)
		{
			ylj+=aa[i][j]*pow(roux,i+1)*pow(be,j);
		}
	}
	//////////////////////////////////////
	ba+=(1-len)*log(ylj)/(len*1.);
	return ba;
}

double LJ::debalj(double rou)
{
	int i,j;
	double roux,deba,f,ylj,deylj;
	double deg[6];
    /////////////////g////////////////////
	roux=rou*dia*dia*dia;
	f=exp(-gama*roux*roux);
	deg[0]=roux*f;
	for(i=1;i<6;i++)
	{
		deg[i]=i/gama*deg[i-1]-f*pow(roux,2*i-1)*(i-gama*roux*roux)/gama;
	}
	//////////////////////////////////////
	for(i=0,deba=0.;i<8;i++)
	{
		deba+=a[i]*pow(roux,i);
	}
	for(i=0;i<6;i++)
	{
		deba+=b[i]*deg[i];
	}
	///////////////ylj////////////////////
	for(i=0,ylj=1.;i<5;i++)
	{
		for(j=0;j<5;j++)
		{
			ylj+=aa[i][j]*pow(roux,i+1)*pow(be,j);
		}
	}
	//////////////////////////////////////
	///////////////deylj//////////////////
	for(i=0,deylj=0.;i<5;i++)
	{
		for(j=0;j<5;j++)
		{
			deylj+=(i+1)*aa[i][j]*pow(roux,i)*pow(be,j);
		}
	}
	//////////////////////////////////////
	deba+=(1-len)*deylj/(len*1.)/ylj;
	return deba*dia*dia*dia;
}*/

double LJ::baex(double *r)
{
	int i,j,k;
	double roux,bex,rous;
	double ba,f,ylj,bach;
	double a[8],b[6],g[6];
	roux=0.;
	rous=0.;
	bex=0.;
	for(i=0;i<kind;i++)
	{
		rous+=r[i];
		for(j=0;j<kind;j++)
		{
			roux+=r[i]*r[j]*dia[i][j]*dia[i][j]*dia[i][j];
			bex+=r[i]*r[j]*dia[i][j]*dia[i][j]*dia[i][j]*eff[i][j];
		}
	}
	roux=roux/rous;
	bex=bex/rous/roux;
	////////////////a/////////////////////
	a[0]=x[0]/bex+x[1]*sqrt(1./bex)+x[2]+x[3]*bex+x[4]*bex*bex;
	a[1]=x[5]/bex+x[6]+x[7]*bex+x[8]*bex*bex;
	a[2]=x[9]/bex+x[10]+x[11]*bex;
	a[3]=x[12];
	a[4]=x[13]*bex+x[14]*bex*bex;
	a[5]=x[15]*bex;
	a[6]=x[16]*bex+x[17]*bex*bex;
	a[7]=x[18]*bex*bex;
	//////////////////////////////////////
	/////////////////b////////////////////
	b[0]=x[19]*bex*bex+x[20]*bex*bex*bex;
	b[1]=x[21]*bex*bex+x[22]*bex*bex*bex*bex;
	b[2]=x[23]*bex*bex+x[24]*bex*bex*bex;
	b[3]=x[25]*bex*bex+x[26]*bex*bex*bex*bex;
	b[4]=x[27]*bex*bex+x[28]*bex*bex*bex;
	b[5]=x[29]*bex*bex+x[30]*bex*bex*bex+x[31]*bex*bex*bex*bex;//cout<<bex<<endl;
	//////////////////////////////////////
	f=exp(-gama*roux*roux);
	g[0]=(1.-f)/(2*gama);
	for(i=1;i<6;i++)
	{
		g[i]=i/gama*g[i-1]-f/(2*gama)*pow(roux,2*i);
		
	}
	//////////////////////////////////////
	for(i=0,ba=0.;i<8;i++)
	{
		ba+=a[i]/(i+1)*pow(roux,i+1);
	}
	for(i=0;i<6;i++)
	{
		ba+=b[i]*g[i];
		
	}
	ba*=bex;
	bach=0.;
	//状态方程中，混合LJ成链部分有问题
	/*for(k=0;k<kind;k++)
	{
		for(i=0,ylj=1.;i<5;i++)
		{
		    for(j=0;j<5;j++)
			{
				ylj+=aa[i][j]*pow(roux,i+1)*pow(be,j);
			}
		}
		bach+=(1-len[k])*log(ylj)/(len[k]*1.)*r[k];
	}
	bach=bach/rous;
	*/
	ba+=bach;
	return ba;
}


double LJ::buex(double *r)
{
	int i,j,k;
	double roux,bex,rous;
	double ba,f,ylj,bach;
	double a[8],b[6],g[6];
	roux=0.;
	rous=0.;
	bex=0.;
	for(i=0;i<kind;i++)
	{
		rous+=r[i];
		for(j=0;j<kind;j++)
		{
			roux+=r[i]*r[j]*dia[i][j]*dia[i][j]*dia[i][j];
			bex+=r[i]*r[j]*dia[i][j]*dia[i][j]*dia[i][j]*eff[i][j];
		}
	}
	roux=roux/rous;
	bex=bex/rous/roux;
	////////////////a/////////////////////
	a[0]=0.5*x[1]*sqrt(1./bex)+x[2]+2*x[3]*bex+3*x[4]*bex*bex;
	a[1]=x[6]+2*x[7]*bex+3*x[8]*bex*bex;
	a[2]=x[10]+2*x[11]*bex;
	a[3]=x[12];
	a[4]=2*x[13]*bex+3*x[14]*bex*bex;
	a[5]=2*x[15]*bex;
	a[6]=2*x[16]*bex+3*x[17]*bex*bex;
	a[7]=3*x[18]*bex*bex;
	//////////////////////////////////////
	/////////////////b////////////////////
	b[0]=3*x[19]*bex*bex+4*x[20]*bex*bex*bex;
	b[1]=3*x[21]*bex*bex+5*x[22]*bex*bex*bex*bex;
	b[2]=3*x[23]*bex*bex+4*x[24]*bex*bex*bex;
	b[3]=3*x[25]*bex*bex+5*x[26]*bex*bex*bex*bex;
	b[4]=3*x[27]*bex*bex+4*x[28]*bex*bex*bex;
	b[5]=3*x[29]*bex*bex+4*x[30]*bex*bex*bex+5*x[31]*bex*bex*bex*bex;//cout<<bex<<endl;
	//////////////////////////////////////
	f=exp(-gama*roux*roux);
	g[0]=(1.-f)/(2*gama);
	for(i=1;i<6;i++)
	{
		g[i]=i/gama*g[i-1]-f/(2*gama)*pow(roux,2*i);
		
	}
	//////////////////////////////////////
	for(i=0,ba=0.;i<8;i++)
	{
		ba+=a[i]/(i+1)*pow(roux,i+1);
	}
	for(i=0;i<6;i++)
	{
		ba+=b[i]*g[i];
		
	}
	ba*=bex;
	bach=0.;
	//状态方程中，混合LJ成链部分有问题
	/*for(k=0;k<kind;k++)
	{
		for(i=0,ylj=1.;i<5;i++)
		{
		    for(j=0;j<5;j++)
			{
				ylj+=aa[i][j]*pow(roux,i+1)*pow(be,j);
			}
		}
		bach+=(1-len[k])*log(ylj)/(len[k]*1.)*r[k];
	}
	bach=bach/rous;
	*/
	ba+=bach;
	return ba;
}

void LJ::debaex(double *r,double *da)//r为各粒子的密度，da为返回的相应粒子的导数
{
	int i	;
	//double roux,bex,rous,sx3;
	//double ba,f,ylj,bach;
	//double a[8],b[6],g[6],c[8],d[6];
	//double px,ux;
	double *dexdr,*dsxdr;
	double *ru,*rd;
	double dr;
	dr=0.0001;
	ru=new double[kind];
	rd=new double[kind];
	for(i=0;i<kind;i++)
	{
		ru[i]=r[i];
		rd[i]=r[i];
	}
	for(i=0;i<kind;i++)
	{
		if(r[i]<MIN)
		{
			da[i]=0.;
		}
		else
		{
			ru[i]=r[i]*(1.+dr);
		    rd[i]=r[i]*(1.-dr);
		    da[i]=(baex(ru)-baex(rd))/(2.*r[i]*dr);
		    ru[i]=r[i];
		    rd[i]=r[i];
		}
	}
	delete []ru;
	delete []rd;
	return;
	
/*int i,j,k;
	double roux,bex,rous,sx3;
	double ba,f,ylj,bach;
	double a[8],b[6],g[6],c[8],d[6];
	double px,ux;
	double *dexdr,*dsxdr;
	dexdr=new double[kind];
	dsxdr=new double[kind];
	sx3=0.;
	rous=0.;
	bex=0.;
	for(i=0;i<kind;i++)
	{
		rous+=r[i];
		for(j=0;j<kind;j++)
		{
			sx3+=r[i]*r[j]*dia[i][j]*dia[i][j]*dia[i][j];
			bex+=r[i]*r[j]*dia[i][j]*dia[i][j]*dia[i][j]*eff[i][j];
		}
	}
	roux=sx3/rous;
	bex=bex/sx3;
	////////////////a/////////////////////
	a[0]=x[0]/bex+x[1]*sqrt(1./bex)+x[2]+x[3]*bex+x[4]*bex*bex;
	a[1]=x[5]/bex+x[6]+x[7]*bex+x[8]*bex*bex;
	a[2]=x[9]/bex+x[10]+x[11]*bex;
	a[3]=x[12];
	a[4]=x[13]*bex+x[14]*bex*bex;
	a[5]=x[15]*bex;
	a[6]=x[16]*bex+x[17]*bex*bex;
	a[7]=x[18]*bex*bex;
	//////////////////////////////////////
	/////////////////b////////////////////
	b[0]=x[19]*bex*bex+x[20]*bex*bex*bex;
	b[1]=x[21]*bex*bex+x[22]*bex*bex*bex*bex;
	b[2]=x[23]*bex*bex+x[24]*bex*bex*bex;
	b[3]=x[25]*bex*bex+x[26]*bex*bex*bex*bex;
	b[4]=x[27]*bex*bex+x[28]*bex*bex*bex;
	b[5]=x[29]*bex*bex+x[30]*bex*bex*bex+x[31]*bex*bex*bex*bex;
	//////////////////////////////////////
	////////////////c/////////////////////
	c[0]=0.5*x[2]*sqrt(1./bex)+x[3]+2*x[4]*bex+3*x[5]*bex*bex;
	c[1]=x[7]+2*x[8]*bex+3*x[9]*bex*bex;
	c[2]=x[11]+2*x[12]*bex;
	c[3]=x[13];
	c[4]=2*x[14]*bex+3*x[15]*bex*bex;
	c[5]=2*x[16]*bex;
	c[6]=2*x[17]*bex+3*x[18]*bex*bex;
	c[7]=3*x[19]*bex*bex;
	//////////////////////////////////////
	////////////////d/////////////////////
	d[0]=3*x[20]*bex*bex+4*x[21]*bex*bex*bex;
	d[1]=3*x[22]*bex*bex+5*x[23]*bex*bex*bex*bex;
	d[2]=3*x[24]*bex*bex+4*x[25]*bex*bex*bex;
	d[3]=3*x[26]*bex*bex+5*x[27]*bex*bex*bex*bex;
	d[4]=3*x[28]*bex*bex+4*x[29]*bex*bex*bex;
	d[5]=3*x[30]*bex*bex+4*x[31]*bex*bex*bex+5*x[32]*bex*bex*bex*bex;
	//////////////////////////////////////
	f=exp(-gama*roux*roux);
	g[0]=(1.-f)/(2*gama);
	for(i=1;i<6;i++)
	{
		g[i]=i/gama*g[i-1]-f/(2*gama)*pow(roux,2*i);
	}
	//////////////////////////////////////
	ba=0.;
	px=roux/bex;
	ux=0.;
	for(i=0;i<8;i++)
	{
		ba+=a[i]/(i+1)*pow(roux,i+1);
		px+=a[i]*pow(roux,i+2);
		ux+=c[i]/(i+1)*pow(roux,i+1);
	}
	for(i=0;i<6;i++)
	{
		ba+=b[i]*g[i];
		px+=f*b[i]*pow(roux,2*i+3);
		ux+=d[i]*g[i];
	}
	////////////////计算dsxdr,dexdr/////////////////
	for(i=0;i<kind;i++)
	{
		dsxdr[i]=0.;
		dexdr[i]=0.;
		for(k=0;k<kind;k++)
		{
			dsxdr[i]+=r[k]*dia[i][k]*dia[i][k]*dia[i][k];
			dexdr[i]+=r[k]*eff[i][k]*dia[i][k]*dia[i][k]*dia[i][k];
		}
		dsxdr[i]=dsxdr[i]/rous;zzzzzzzzzzzzzzzzzzzzzzz
		dsxdr[i]+=-sx3;
		dsxdr[i]*=2./rous;
		///
		dexdr[i]=dexdr[i]/roux;
		dexdr[i]+=-bex;
		dexdr[i]*=2./rous;
		dexdr[i]+=-bex/sx3*dsxdr[i];
	}
	/////////////////////////////////////////////////
	///////////////////单体项贡献/////////////////////
	for(i=0;i<kind;i++)
	{
		da[i]=px/roux/roux-1./bex/roux;
		da[i]*=sx3+rous*dsxdr[i];
		da[i]+=(ux-ba)/bex*dexdr[i];
		da[i]*=bex;
	}
	//////////////////////////////////////////////////
	//////////////////成链项贡献//////////////////////
	//////////////////////////////////////////////////
	delete []dexdr;
	delete []dsxdr;*/
}

void LJ::r_chem(double *r,double *chem)
{
	int i;
	double *da, totr;
	da=new double[kind];
	debaex(r,da);
	totr=0.;
	for(i=0;i<kind;i++)
	{
		totr+=r[i];
	}
	for(i=0;i<kind;i++)
	{
		chem[i]=baex(r)+totr*da[i]+log(r[i]);
	}
	delete []da;
}

void LJ::chem_r(double *chem,double *r)
{
	int i;
	double dd;
	double *rd,*ru,*mu,*md,*rn,*mn;
	rd=new double[kind];
	ru=new double[kind];
	mu=new double[kind];
	md=new double[kind];
	rn=new double[kind];
	mn=new double[kind];
	if(kind==1)
	{
		rn[0]=exp(chem[0]);
		r_chem(rn,mn);
		rd[0]=rn[0]/ZOOM;
		ru[0]=rn[0]*ZOOM;
		r_chem(rd,md);
		r_chem(ru,mu);
		for(i=0;;i++)
		{
			//rn[0]=pow(pow(rd[0],mu[0]-chem[0])*pow(ru[0],chem[0]-md[0]),1./(mu[0]-md[0]));
			rn[0]=(rd[0]+ru[0])*0.5;
			r_chem(rn,mn);
			if(mn[0]-chem[0]>0.)
			{
				ru[0]=rn[0];
				mu[0]=mn[0];
			}
			else
			{
				rd[0]=rn[0];
				md[0]=mn[0];
			}
			if(juedui(ru[0]-rd[0])<PRE)
			{
				cout<<"in precision"<<endl;
				break;
			}
		}
		cout<<ru[0]-rd[0]<<endl;
		r[0]=0.5*(ru[0]+rd[0]);
	}
	delete []rd;
	delete []ru;
	delete []mu;
	delete []md;
	delete []rn;
	delete []mn;
}

double LJ::press(double *r)
{
	double p,*da;
	int i,j,k;
	double roux,bex,rous;
	double ba,f,ylj,bach;
	double a[8],b[6],g[6];
	if(kind==1)
	{
		
	    da=new double[kind];
	    debaex(r,da);
	    p=r[0]+r[0]*r[0]*da[0];
	    p*=1.3626E2*tem;
	    return p;
	}
	////////////////////////////////////////
	roux=0.;
	rous=0.;
	bex=0.;
	for(i=0;i<kind;i++)
	{
		rous+=r[i];
		for(j=0;j<kind;j++)
		{
			roux+=r[i]*r[j]*dia[i][j]*dia[i][j]*dia[i][j];
			bex+=r[i]*r[j]*dia[i][j]*dia[i][j]*dia[i][j]*eff[i][j];
		}
	}
	roux=roux/rous;
	bex=bex/rous/roux;
	////////////////a/////////////////////
	a[0]=x[0]/bex+x[1]*sqrt(1./bex)+x[2]+x[3]*bex+x[4]*bex*bex;
	a[1]=x[5]/bex+x[6]+x[7]*bex+x[8]*bex*bex;
	a[2]=x[9]/bex+x[10]+x[11]*bex;
	a[3]=x[12];
	a[4]=x[13]*bex+x[14]*bex*bex;
	a[5]=x[15]*bex;
	a[6]=x[16]*bex+x[17]*bex*bex;
	a[7]=x[18]*bex*bex;
	//////////////////////////////////////
	/////////////////b////////////////////
	b[0]=x[19]*bex*bex+x[20]*bex*bex*bex;
	b[1]=x[21]*bex*bex+x[22]*bex*bex*bex*bex;
	b[2]=x[23]*bex*bex+x[24]*bex*bex*bex;
	b[3]=x[25]*bex*bex+x[26]*bex*bex*bex*bex;
	b[4]=x[27]*bex*bex+x[28]*bex*bex*bex;
	b[5]=x[29]*bex*bex+x[30]*bex*bex*bex+x[31]*bex*bex*bex*bex;//cout<<bex<<endl;
	//////////////////////////////////////
	f=exp(-gama*roux*roux);
	g[0]=(1.-f)/(2*gama);
	for(i=1;i<6;i++)
	{
		g[i]=i/gama*g[i-1]-f/(2*gama)*pow(roux,2*i);
		
	}
	//////////////////////////////////////
	for(i=0,ba=0.;i<8;i++)
	{
		ba+=a[i]*pow(roux,i+2);
	}
	for(i=0;i<6;i++)
	{
		ba+=b[i]*f*pow(roux,2*i+3);
		
	}
	ba+=roux;
	ba*=rous/roux;
	ba*=1.3626E2*tem;
	return ba;
}

double LJ::s_per_n(double *r)
{
	int i;
	double tbak,bf1,bf2,totr,re;
	tbak=tem;
	///////////////
	bf1=0.;
	totr=0.;
	for(i=0;i<kind;i++)
	{
		bf1+=r[i]*(log(r[i])-1.);
		totr+=r[i];
	}
	bf1=bf1/totr;
	bf2=bf1;//理想项
	tem=tbak*(1.+DE);
	bf1+=baex(r);
	bf1*=tem;
	tem=tbak*(1.-DE);
	bf2+=baex(r);
	bf2*=tem;
	re=(bf2-bf1)/(2.*DE*tbak)*8.31;
	///////////////
	tem=tbak;
	return re;
}

double LJ::sex_nk(double *r)
{
  int i;
  double tbak,bf1,bf2,re;
  tbak=tem;
  tem=tbak*(1.+DE);
  bf1=baex(r)*tem;
  tem=tbak*(1.-DE);
  bf2=baex(r)*tem;
  re=(bf2-bf1)/(2.*DE*tbak);
  tem=tbak;
  return re;
}


double LJ::f_per_n(double *r)
{
	int i;
	double totr,re;
	re=0.;
	totr=0.;
	for(i=0;i<kind;i++)
	{
		re+=r[i]*(log(r[i])-1.);
		totr+=r[i];
	}
	re=re/totr;
	re+=baex(r);
    re*=tem*8.31;
	return re;
}

double slj_nk(double sig,double pot,double tp,double* r)
{
  double efact,t1,t2,f1,f2,pot1,pot2,re;
  efact=pot*tp;
  t1=tp*(1.0-DE);
  pot1=efact/t1;
  LJ ref1(sig,pot1,t1);
  f1=ref1.baex(r)*t1;
  
  t2=tp*(1.0+DE);
  pot2=efact/t2;
  LJ ref2(sig,pot2,t2);
  f2=ref2.baex(r)*t2;
  
  re=(f1-f2)/(2*DE*tp);
  return re;
}

////////////////////////////////////////////////////////////////////////////////

//JCP 131,024704 (2009)

FLJ::FLJ(double pot, double dia)
{
  double ds,k0,tstar;
  be=pot;
  sig=dia;
  
  tstar=1.0/be;
  ds=(1.0+0.2977*tstar)/(1.0+0.33163*tstar+1.0477E-3*tstar*tstar);//d/sigma
  /////////
  //ds=lj_to_hs(tstar);
  /////////
  //ds=1.0;
  d=sig*ds;//d
  
  lam[0]=2.9637*ds;
  lam[1]=14.0167*ds;
  k0=2.1717/ds;
  bey[0]=k0*exp(lam[0]*(sig/d-1.0));
  bey[1]=-k0*exp(lam[1]*(sig/d-1.0));
  
  //// besd1=48be*[1/9(sigma/d)^12-1/3(sigma/d)^6]
  besd1=pow(ds,-12)/9-pow(ds,-6)/3;
  besd1*=48*be;
  //// besd2=24be*[1/9(sigma/d)^12-1/3(sigma/d)^6+2/9(sigma/d)^3]
  besd2=pow(ds,-12)/9-pow(ds,-6)/3+2*pow(ds,-3)/9;
  besd2*=24*be;
}


double FLJ::run(double rou)
{
  int i,j;
  double f1,f2,fmf; 
  double y,y1,y2,y5,ghs;
  double q[2],l[2],s[2];
  double tp;//
  
  y=PI*rou*d*d*d/6;
  y1=1.0-y;
  y2=1.0+2*y;
  y5=1.0+0.5*y;

  //ghs=(1.0-0.5*y)/y1/y1/y1;
  ghs=(1.0+0.5*y)/y1/y1;
  ///////////////s,l,q//////////////////  
  for(i=0;i<2;i++)
  {
    s[i]=y1*y1*pow(lam[i],3)+6*y*y1*lam[i]*lam[i]+18*y*y*lam[i]-12*y*y2;
    l[i]=y5*lam[i]+y2;
    q[i]=s[i]+12*y*l[i]*exp(-lam[i]);
    q[i]=q[i]/y1/y1/pow(lam[i],3);
  }
  //////////////////////////////////////
  
  /////////////f1//////////////
  f1=0.0;
  for(i=0;i<2;i++)
  {
    tp=l[i]/y1/y1/q[i]-1.0-lam[i];
    tp*=bey[i]/lam[i]/lam[i];
    f1+=tp;
  }
  f1*=-12*y;
  f1+=y*besd1-2*besd2*y*ghs;
  /////////////////////////////
  
  ////////////f2///////////////
  f2=0.0;
  tp=0.0;
  for(i=0;i<2;i++)
  {
    for(j=0;j<2;j++)
    {
      f2+=bey[i]*bey[j]/(lam[i]+lam[j])/pow(q[i]*q[j],2);
    }
    tp+=bey[i]/q[i]/q[i];
  }
  f2=-6*y*f2-tp*besd2*y;
  /////////////////////////////
  
  fmf=-(16.0/9.0)*PI*be*rou*sig*sig*sig;
  
  return f1+f2-fmf;
}



void Vector::initial(double a, double b,double c)
{
  x=a;
  y=b;
  z=c;
}

double Vector::rabs()
{
  return sqrt(x*x+y*y+z*z);
}

Vector operator+(Vector a, Vector b)
{
  Vector c;
  c.x=a.x+b.x;
  c.y=a.y+b.y;
  c.z=a.z+b.z;
  return c;
}

double operator^(Vector a, Vector b)
{
  return a.x*b.x+a.y*b.y+a.z*b.z;
}

Vector operator*(Vector a, double b)
{
  Vector c;
  c.x=a.x*b;
  c.y=a.y*b;
  c.z=a.z*b;
  return c;
}

Vector operator*(double b,Vector a)
{
  Vector c;
  c.x=a.x*b;
  c.y=a.y*b;
  c.z=a.z*b;
  return c;
}

Vector operator/(Vector a, double b)
{
  Vector c;
  c.x=a.x/b;
  c.y=a.y/b;
  c.z=a.z/b;
  return c;
}

Vector operator*(Vector a, Vector b)
{
  Vector c;
  c.x=a.y*b.z-a.z*b.y;
  c.y=a.z*b.x-a.x*b.z;
  c.z=a.x*b.y-a.y*b.x;
  return c;
}



void Weight_K::initial(double x)
{
  int i;
  double kr,ks2;
  
  dia=x;
  nw=5000;
  dw=0.01;
  
  for(i=1;i<nw;i++)
  {
    kr=i*dw;
    ks2=kr*dia*0.5;
    fw0[i]=sin(ks2)/ks2;
    fw3[i]=4*PI*(sin(ks2)-ks2*cos(ks2))/kr/kr/kr;
  }
  fw0[0]=1.0;
  fw3[0]=PI/6*dia*dia*dia;
    
  for(i=1;i<nw;i++)// wcor
  {
    kr=i*dw;
    ks2=kr*dia;
    fwcor[i]=3*(sin(ks2)-ks2*cos(ks2))/ks2/ks2/ks2;
  }
  fwcor[0]=1.0;
  
  return;
}

double Weight_K::w0(double r)
{
  int i;
  double dr;
  if(r>=(nw-1)*dw)
  {
    return 0.0;
  }
  i=int(r/dw);
  dr=r-i*dw;
  return fw0[i]*(1.0-dr)+fw0[i+1]*dr;
}

double Weight_K::w1(double r)
{
  return w0(r)*0.5*dia;
}

double Weight_K::w2(double r)
{
  return w0(r)*PI*dia*dia;
}

double Weight_K::w3(double r)
{
  int i;
  double dr;
  if(r>=(nw-1)*dw)
  {
    return 0.0;
  }
  i=int(r/dw);
  dr=r-i*dw;
  return fw3[i]*(1.0-dr)+fw3[i+1]*dr;
}

Vector Weight_K::wv2(Vector vr)
{
  double r,wv2_abs;
  Vector vn;
  
  r=vr.rabs();
  wv2_abs=-w3(r);
  
  if(r<=1.0e-10)
  {
    vn.x=0.0;
    vn.y=0.0;
    vn.z=0.0;
  }
  else
  {
    vn.x=vr.x*wv2_abs;
    vn.y=vr.y*wv2_abs;
    vn.z=vr.z*wv2_abs;
  }
  return vn;
}

Vector Weight_K::wv1(Vector vr)
{
  double r,wv2_abs;
  Vector vn;
  
  r=vr.rabs();
  wv2_abs=-w3(r)/(2*PI*dia);
  
  if(r<=1.0e-10)
  {
    vn.x=0.0;
    vn.y=0.0;
    vn.z=0.0;
  }
  else
  {
    vn.x=vr.x*wv2_abs;
    vn.y=vr.y*wv2_abs;
    vn.z=vr.z*wv2_abs;
  }
  return vn;
}



DFT::DFT()
{
  int i,j,k,n_atom,k_atom,n,m,n_atom_in,npd,ig;
  char tip[50];
  double tstar,yita; 
  double dmix,emix,x,y,z,r,x1,y1,z1,x2,y2,z2;
  bool find_atom;
  int readv;
  int cutex;
  int ix1,ix2,iy1,iy2,iz1,iz2,ix,iy,iz;
  double cutoff2,cuty,cutz,dgrid;
  int maxnx,maxny,maxnz;
  Vector a12,a23,a31;
  double atot,dtem,shift,mv;
  long nijk;
  double acdens;
  int nspecial;
  double e_in,d_in;
  
  start=clock();
  
  ifstream ip("input.dat",ios::in);
  ofstream op("output.dat",ios::out);
  
  ip>>tip>>tip>>tip;
  ip>>maxnx>>maxny>>maxnz;
  ip>>tip>>tip>>tip>>tip;
  ip>>lxu>>lyu>>lzu>>dgrid;
  ip>>tip>>tip>>tip;
  ip>>angle_a>>angle_b>>angle_c;  
  ip>>tip>>tip>>tip>>tip;
  ip>>tem>>kapa>>detarou>>torr;
  ip>>tip;
  ip>>k_gas;
  be=new double[k_gas];
  sig=new double[k_gas];
  rb=new double[k_gas];
  bmhs=new double[k_gas];
  dbh=new double[k_gas];
  name_gas=new int[k_gas];
  wfun=new Weight_K[k_gas];
  rouav=new double[k_gas];
  max_dens=new double[k_gas];
  
  bemix=new double*[k_gas];
  sigmix=new double*[k_gas];
  

  ip>>tip>>tip>>tip>>tip;
  for(i=0;i<k_gas;i++)
  {
	  ip>>name_gas[i]>>be[i]>>sig[i]>>rb[i];
	  bemix[i]=new double[k_gas];
	  sigmix[i]=new double[k_gas];
  }
  
  for(i=0;i<k_gas;i++)
  {
    for(j=0;j<k_gas;j++)
    {
      bemix[i][j]=sqrt(be[i]*be[j]);
      sigmix[i][j]=(sig[i]+sig[j])/2;
    }
  }
  
  ip>>tip;
  ip>>nspecial;
  ip>>tip>>tip>>tip>>tip;
  for(k=0;k<nspecial;k++)
  {
    ip>>i>>j>>e_in>>d_in;
    bemix[i][j]=e_in;
    sigmix[i][j]=d_in;
    bemix[j][i]=bemix[i][j];
    sigmix[j][i]=sigmix[j][i];
  }
  
  ip>>tip>>tip>>tip>>tip>>tip;
  ip>>readv>>eos>>cutoff>>deos>>fix_cut;
  ip>>tip>>tip>>tip;
  ip>>mass>>knudsen>>p_mfv;
  ip>>tip;
  ip>>k_atom;
  
  //cout<<mass<<" "<<knudsen<<" "<<p_mfv<<endl;
  
  angle_a*=PI/180.0;
  angle_b*=PI/180.0;
  angle_c*=PI/180.0;
  
  a1.initial(1.0,0.0,0.0);
  a2.initial(cos(angle_c),sin(angle_c),0.0);
  x=cos(angle_b);
  y=(cos(angle_a)-cos(angle_b)*cos(angle_c))/sin(angle_c);
  z=pow(sin(angle_c),2)-pow(cos(angle_a),2)-pow(cos(angle_b),2);
  z+=2*cos(angle_a)*cos(angle_b)*cos(angle_c);
  atot=sqrt(z);
  z=atot/sin(angle_c);
  a3.initial(x,y,z);
  
  cosa2=cos(angle_a)*2;
  cosb2=cos(angle_b)*2;
  cosc2=cos(angle_c)*2;
  
  a12=a1*a2;
  a23=a2*a3;
  a31=a3*a1;
  a123=a1^a23;
  b1=a23/a23.rabs();
  b2=a31/a31.rabs();
  b3=a12/a12.rabs();
  
  
  
  cout<<b1.x<<" "<<b1.y<<" "<<b1.z<<endl;
  
  pdx=int(cutoff*2/lxu)+1;
  pdy=int(cutoff*2/lyu)+1;
  pdz=int(cutoff*2/lzu)+1;
  //pdx=1.0;
  //pdy=1.0;
  //pdz=1.0;
  pdxyz=pdx*pdy*pdz;
  cout<<pdx<<" "<<pdy<<" "<<pdz<<endl;
  lx=lxu*pdx;
  ly=lyu*pdy;
  lz=lzu*pdz;
  cout<<lx<<" "<<ly<<" "<<lz<<endl;
  nx=int(lx/dgrid)+1;
  ny=int(ly/dgrid)+1;
  nz=int(lz/dgrid)+1;
  if(nx>maxnx)nx=maxnx;
  if(ny>maxny)ny=maxny;
  if(nz>maxnz)nz=maxnz;
  cout<<nx<<" "<<ny<<" "<<nz<<endl;
  
  pot=new double[k_atom];
  dia=new double[k_atom];
  id=new int[k_atom];
  
  ip>>tip>>tip>>tip;
  
  for(i=0;i<k_atom;i++)
  {
    ip>>id[i]>>dia[i]>>pot[i];
    pot[i]=pot[i]/tem;
  }
  
  ip>>tip>>tip>>tip;
  ip>>n_atom_in>>mmof>>matom;
  
  n_atom=n_atom_in*pdxyz;
  aid=new int[n_atom];
  ax=new double[n_atom];
  ay=new double[n_atom];
  az=new double[n_atom];
  
  ip>>tip>>tip>>tip>>tip;
  
  for(n=0;n<n_atom_in;n++)
  {
    ip>>aid[n]>>ax[n]>>ay[n]>>az[n];
    ax[n]=dmod(ax[n],lxu);
    ay[n]=dmod(ay[n],lyu);
    az[n]=dmod(az[n],lzu);
    //cout<<aid[n]<<" "<<ax[n]<<" "<<ay[n]<<" "<<az[n]<<endl;
    for(i=0;i<pdx;i++)
    {
      for(j=0;j<pdy;j++)
      {
	for(k=0;k<pdz;k++)
	{
	  npd=((i*pdy+j)*pdz+k)*n_atom_in+n;
	  if(i!=0||j!=0||k!=0)
	  {	    
	    aid[npd]=aid[n];
	    ax[npd]=ax[n]+i*lxu;
	    ay[npd]=ay[n]+j*lyu;
	    az[npd]=az[n]+k*lzu;
	  }
	  //if(n==0)cout<<'\n'<<((i*pdy+j)*pdz+k);
	}
      }
    }
  }
  
  nfa=10000;
  dfa=3.0e-6;
   
  dx=lx/nx;
  dy=ly/ny;
  dz=lz/nz;
  dv=dx*dy*dz;
  rcut=lx/2;
  if(rcut>ly/2)rcut=ly/2;
  if(rcut>lz/2)rcut=lz/2;
  dkx=2*PI/lx*sin(angle_a)/atot;
  dky=2*PI/ly*sin(angle_b)/atot;
  dkz=2*PI/lz*sin(angle_c)/atot;
  nxyz=nx*ny*nz;

  cout<<dkx<<" "<<dky<<" "<<dkz<<endl;
  
  for(i=0;i<k_gas;i++)
  {
	  be[i]=be[i]/tem;
	  tstar=1.0/be[i];
	  dbh[i]=(1.0+0.2977*tstar)/(1.0+0.33163*tstar+1.0477E-3*tstar*tstar)*sig[i];
	  wfun[i].initial(dbh[i]);
	  //yita=0.5235987755983*rb*dbh*dbh*dbh;//1/6*pi*rou
  //bmhs=(8.0-9*yita+3*yita*yita)*yita/pow(1-yita,3);
	  for(j=0;j<k_gas;j++)
	  {
	    bemix[i][j]=bemix[i][j]/tem;
	  }
  }
  
  fatt=new double[nfa];
  rdfatt=new double[nfa];
  ds_rou=new double[nfa];
  
  vext=new double***[nx];
  //dens=new double**[nx];
  denskr=new double***[nx];
  denski=new double***[nx];
  uattk=new double****[nx];
  
  tw0=new double***[nx];
  //tw1=new double***[nx];
  //tw2=new double***[nx];
  tw3=new double***[nx];
  twcor=new double***[nx];
  //tw1x=new double***[nx];
  //tw1y=new double***[nx];
  //tw1z=new double***[nx];
  tw2x=new double***[nx];
  tw2y=new double***[nx];
  tw2z=new double***[nx];
  
  m_x = new double [nx*ny*nz*k_gas];
  
  n0=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  n1=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  n2=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  n3=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  nv1x=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  nv1y=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  nv1z=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  nv2x=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  nv2y=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  nv2z=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  cor=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  fft_in=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  fft_out=(fftw_complex *)malloc(nx*ny*nz*sizeof(fftw_complex));
  for(i=0;i<nx;i++)
  {
    vext[i]=new double**[ny];
	denskr[i]=new double**[ny];
	denski[i]=new double**[ny];
    uattk[i]=new double***[ny];
    tw0[i]=new double**[ny];
	//tw1[i]=new double**[ny];
	//tw2[i]=new double**[ny];
    tw3[i]=new double**[ny];
    twcor[i]=new double**[ny];
	//tw1x[i]=new double**[ny];
    //tw1y[i]=new double**[ny];
    //tw1z[i]=new double**[ny];
    tw2x[i]=new double**[ny];
    tw2y[i]=new double**[ny];
    tw2z[i]=new double**[ny];
    for(j=0;j<ny;j++)
    {
      vext[i][j]=new double*[nz];
	  denskr[i][j]=new double*[nz];
	  denski[i][j]=new double*[nz];
      uattk[i][j]=new double**[nz];
      tw0[i][j]=new double*[nz];
	  //tw1[i][j]=new double*[nz];
	  //tw2[i][j]=new double*[nz];
      tw3[i][j]=new double*[nz];
      twcor[i][j]=new double*[nz];
	  //tw1x[i][j]=new double*[nz];
      //tw1y[i][j]=new double*[nz];
      //tw1z[i][j]=new double*[nz];
      tw2x[i][j]=new double*[nz];
      tw2y[i][j]=new double*[nz];
      tw2z[i][j]=new double*[nz];
      for(k=0;k<nz;k++)
      {
		  vext[i][j][k]=new double[k_gas];
		  denskr[i][j][k]=new double[k_gas];
		  denski[i][j][k]=new double[k_gas];
		  uattk[i][j][k]=new double*[k_gas];
		  tw0[i][j][k]=new double[k_gas];
		  //tw1[i][j][k]=new double[k_gas];
		  //tw2[i][j][k]=new double[k_gas];
		  tw3[i][j][k]=new double[k_gas];
		  twcor[i][j][k]=new double[k_gas];
		  //tw1x[i][j][k]=new double[k_gas];
		  //tw1y[i][j][k]=new double[k_gas];
		  //tw1z[i][j][k]=new double[k_gas];
		  tw2x[i][j][k]=new double[k_gas];
		  tw2y[i][j][k]=new double[k_gas];
		  tw2z[i][j][k]=new double[k_gas];
		  for(n=0;n<k_gas;n++)
		  {
			  vext[i][j][k][n]=0.0;
			  uattk[i][j][k][n]=new double[n+1];
		  }
      }
    }
  }
  fft_for=fftw_plan_dft_3d(nx,ny,nz,fft_in,fft_out,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_bak=fftw_plan_dft_3d(nx,ny,nz,fft_in,fft_out,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_n0_for=fftw_plan_dft_3d(nx,ny,nz,n0,n0,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_n1_for=fftw_plan_dft_3d(nx,ny,nz,n1,n1,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_n2_for=fftw_plan_dft_3d(nx,ny,nz,n2,n2,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_n3_for=fftw_plan_dft_3d(nx,ny,nz,n3,n3,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_nv1x_for=fftw_plan_dft_3d(nx,ny,nz,nv1x,nv1x,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_nv1y_for=fftw_plan_dft_3d(nx,ny,nz,nv1y,nv1y,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_nv1z_for=fftw_plan_dft_3d(nx,ny,nz,nv1z,nv1z,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_nv2x_for=fftw_plan_dft_3d(nx,ny,nz,nv2x,nv2x,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_nv2y_for=fftw_plan_dft_3d(nx,ny,nz,nv2y,nv2y,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_nv2z_for=fftw_plan_dft_3d(nx,ny,nz,nv2z,nv2z,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_cor_for=fftw_plan_dft_3d(nx,ny,nz,cor,cor,FFTW_FORWARD,FFTW_ESTIMATE);
  fft_n0_bak=fftw_plan_dft_3d(nx,ny,nz,n0,n0,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_n1_bak=fftw_plan_dft_3d(nx,ny,nz,n1,n1,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_n2_bak=fftw_plan_dft_3d(nx,ny,nz,n2,n2,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_n3_bak=fftw_plan_dft_3d(nx,ny,nz,n3,n3,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_nv1x_bak=fftw_plan_dft_3d(nx,ny,nz,nv1x,nv1x,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_nv1y_bak=fftw_plan_dft_3d(nx,ny,nz,nv1y,nv1y,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_nv1z_bak=fftw_plan_dft_3d(nx,ny,nz,nv1z,nv1z,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_nv2x_bak=fftw_plan_dft_3d(nx,ny,nz,nv2x,nv2x,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_nv2y_bak=fftw_plan_dft_3d(nx,ny,nz,nv2y,nv2y,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_nv2z_bak=fftw_plan_dft_3d(nx,ny,nz,nv2z,nv2z,FFTW_BACKWARD,FFTW_ESTIMATE);
  fft_cor_bak=fftw_plan_dft_3d(nx,ny,nz,cor,cor,FFTW_BACKWARD,FFTW_ESTIMATE);
  ip.close();
  op<<"Box: "<<lx<<"*"<<ly<<"*"<<lz<<"*"<<a123<<"="<<lx*ly*lz*a123<<"A^3"<<'\n';
  op<<"Grid: "<<nx<<"*"<<ny<<"*"<<nz<<"="<<nxyz<<'\n';
  for(i=0;i<k_gas;i++)
  {
	  op<<"Gas "<<i<<': ';
	  op<<"bulk density: "<<rb[i]<<" molec/A^3   ";
	  op<<"sigma: "<<sig[i]<<"A     beta*epsilon:"<<be[i];
  }
  op.close();
  
  
  conv=true;
  
  double dout,eout,mx_cut;
  mx_cut=1.0/dbh[0]/dbh[0]/dbh[0]/k_gas;

  if(readv==0)
  {
	  //////////////////////////////////////
	  cutoff2=cutoff*cutoff;
	  for(n=0;n<n_atom;n++)
	  {
		  find_atom=false;
		  for(m=0;m<k_atom;m++)
		  {
			  if(id[m]==aid[n])
			  {
				  dout=dia[m];//0.5*(sig+dia[m]);
				  eout=pot[m];//sqrt(be*pot[m]);
				  find_atom=true;
				  break;
			  }
		  }
		  if(!find_atom)
		  {
			  dout=0.0;
			  eout=0.0;
			  continue;
		  }
		  for(ig=0;ig<k_gas;ig++)
		  {
			  dmix=(dout+sig[ig])*0.5;
			  emix=sqrt(eout*be[ig]);
			  shift=pow(dmix/cutoff,6);
			  shift=4*emix*(shift*shift-shift);
			  for(i=0;i<nx;i++)
			  {
				  //x=i*dx;
				  x=(i+0.5)*dx-ax[n];
				  if(x>lx/2)
				  {
					  x=x-lx;
				  }
				  if(x<-lx/2)
				  {
					  x=x+lx;
				  }
				  for(j=0;j<ny;j++)
				  {
					  //y=j*dy;
					  y=(j+0.5)*dy-ay[n];
					  if(y>ly/2)
					  {
						  y=y-ly;
					  }
					  if(y<-ly/2)
					  {
						  y=y+ly;
					  }
					  for(k=0;k<nz;k++)
					  {
						  //z=k*dz;
						  z=(k+0.5)*dz-az[n];
						  if(z>lz/2)
						  {
							  z=z-lz;
						  }
						  if(z<-lz/2)
						  {
							  z=z+lz;
						  }
						  r=dis_tri_r(x,y,z);
						  if(r<cutoff)
						  {			  
							  if(r<dmix*0.1)
							  {
								  r=dmix*0.1;
							  }
							  r=pow(dmix/r,6);
							  vext[i][j][k][ig]+=4*emix*(r*r-r)-shift;
						  }
						  ///////////////////////////////////////
					  }
				  }
			  }
		  }
	  }
	  cout<<"writing Vext.dat"<<endl;
	  ofstream vop("Vext.dat",ios::out);
	  for(ig=0;ig<k_gas;ig++)
	  {
		  rouav[ig]=0.0;
	  }
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  vop<<i*dx<<" "<<j*dy<<" "<<k*dz;
				  for(ig=0;ig<k_gas;ig++)
				  {				    
					  nijk=n_mx(i,j,k,ig);
					  /*if(k<nz/4||k>nz*3/4)
					  {
					    vext[i][j][k][ig]=1E7;
					  }*/
					  vop<<" "<<vext[i][j][k][ig]*tem;
					  m_x[nijk]=rb[ig]*exp(-vext[i][j][k][ig]);
					  if(m_x[nijk]>mx_cut)
					  {
						  m_x[nijk]=mx_cut;
					  }
					  rouav[ig]+=m_x[nijk];
					  m_x[nijk]=sqrt(m_x[nijk]);
				  }
				  vop<<'\n';
			  }
		  }
	  }  
  vop.close();  
/////////////////////////////////////////////////////
  }
  else
  {
	  /////////////////////////////////////////////////////    
	  ifstream vop("Vext.dat",ios::in);
	  for(ig=0;ig<k_gas;ig++)
	  {
		  rouav[ig]=0.0;
	  }
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  vop>>x>>y>>z;
				  for(ig=0;ig<k_gas;ig++)
				  {
					  nijk=n_mx(i,j,k,ig);
					  vop>>vext[i][j][k][ig];
					  vext[i][j][k][ig]=vext[i][j][k][ig]/tem;
					  m_x[nijk]=rb[ig]*exp(-vext[i][j][k][ig]);
					  if(m_x[nijk]>mx_cut)
					  {
						  m_x[nijk]=mx_cut;
					  }
					  rouav[ig]+=m_x[nijk];
					  m_x[nijk]=sqrt(m_x[nijk]);
				  }
			  }
		  }
	  }
	  vop.close(); 
//////////////////////////////////////////////////////  
  }
  for(ig=0;ig<k_gas;ig++)
  {
	  rouav[ig]=rouav[ig]/nxyz;
  }
}

long DFT::n_mx(int i,int j,int k,int m)
{
	return k+nz*(j+ny*i)+m*nxyz;
}

DFT::~DFT()
{
  int i,j,k,n;

  fftw_destroy_plan(fft_for);
  fftw_destroy_plan(fft_bak);
  fftw_destroy_plan(fft_n0_for);
  fftw_destroy_plan(fft_n3_for);
  fftw_destroy_plan(fft_nv2x_for);
  fftw_destroy_plan(fft_nv2y_for);
  fftw_destroy_plan(fft_nv2z_for);
  fftw_destroy_plan(fft_cor_for);
  fftw_destroy_plan(fft_n0_bak);
  fftw_destroy_plan(fft_n3_bak);
  fftw_destroy_plan(fft_nv2x_bak);
  fftw_destroy_plan(fft_nv2y_bak);
  fftw_destroy_plan(fft_nv2z_bak);
  fftw_destroy_plan(fft_cor_bak);
  cout<<"plan free ok"<<endl;//test
  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
		for(k=0;k<nz;k++)
		{
			for(n=0;n<k_gas;n++)
			{
				delete []uattk[i][j][k][n];
			}
			delete []vext[i][j][k];
			delete []uattk[i][j][k];
			delete []tw0[i][j][k];
			delete []tw3[i][j][k];
			//delete []tw1[i][j][k];
			//delete []tw2[i][j][k];
			//delete []tw1x[i][j][k];
			//delete []tw1y[i][j][k];
			//delete []tw1z[i][j][k];
			delete []tw2x[i][j][k];
			delete []tw2y[i][j][k];
			delete []tw2z[i][j][k];
		}
      delete []vext[i][j];
      delete []uattk[i][j];
      delete []tw0[i][j];
      delete []tw3[i][j];
	  //delete []tw1[i][j];
      //delete []tw2[i][j];
      //delete []tw1x[i][j];
      //delete []tw1y[i][j];
      //delete []tw1z[i][j];
      delete []tw2x[i][j];
      delete []tw2y[i][j];
      delete []tw2z[i][j];
    }
    delete []vext[i];
    delete []uattk[i];
    delete []tw0[i];
    delete []tw3[i];
	//delete []tw1[i];
    //delete []tw2[i];
    //delete []tw1x[i];
    //delete []tw1y[i];
    //delete []tw1z[i];
    delete []tw2x[i];
    delete []tw2y[i];
    delete []tw2z[i];
  }
  delete []vext;
  delete []uattk;
  delete []fatt;
  delete []rdfatt;
  delete []ds_rou;
  delete []max_dens;
  delete []tw0;
  delete []tw3;
  //delete []tw1;
  //delete []tw2;
  //delete []tw1x;
  //delete []tw1y;
  //delete []tw1z;
  delete []tw2x;
  delete []tw2y;
  delete []tw2z;
  delete []be;
  delete []sig;
  delete []rb;
  delete []bmhs;
  delete []dbh;
  delete []name_gas;
  delete []wfun;
  delete []rouav;
  delete []pot;
  delete []dia;
  delete []id;
  delete []aid;
  delete []ax;
  delete []ay;
  delete []az;
  for(i=0;i<k_gas;i++)
  {
    delete []bemix[i];
    delete []sigmix[i];
  }
  delete []bemix;
  delete []sigmix;
  cout<<"array free ok"<<endl;//test
  if(n0!=NULL)fftw_free(n0);
  if(n3!=NULL)fftw_free(n3);
  if(n1!=NULL)fftw_free(n1);
  if(n2!=NULL)fftw_free(n2);
  if(nv1x!=NULL)fftw_free(nv1x);
  if(nv1y!=NULL)fftw_free(nv1y);
  if(nv1z!=NULL)fftw_free(nv1z);
  if(nv2x!=NULL)fftw_free(nv2x);
  if(nv2y!=NULL)fftw_free(nv2y);
  if(nv2z!=NULL)fftw_free(nv2z);
  if(fft_in!=NULL)fftw_free(fft_in);
  if(fft_out!=NULL)fftw_free(fft_out);
  if(cor!=NULL)fftw_free(cor);
  cout<<"fftw free ok"<<endl;//test
  delete []m_x;
  cout<<"mx free ok"<<endl;//test
}


int DFT::norx(int i)
{
  int ni;
  ni=i%nx;
  if(ni<0)ni+=nx;
  return ni;
}
int DFT::nory(int i)
{
  int ni;
  ni=i%ny;
  if(ni<0)ni+=ny;
  return ni;
}
int DFT::norz(int i)
{
  int ni;
  ni=i%nz;
  if(ni<0)ni+=nz;
  return ni;
}

Vector DFT::tr_cu_k(double xx,double yy,double zz)
{
  Vector cu;
  cu=xx*b1+yy*b2+zz*b3;
  return cu;
}

Vector DFT::tr_cu_k(Vector tr)
{
  return tr_cu_k(tr.x,tr.y,tr.z);
}

Vector DFT::tr_cu_r(double xx,double yy,double zz)
{
  Vector cu;
  cu=xx*a1+yy*a2+zz*a3;
  return cu;
}

Vector DFT::tr_cu_r(Vector tr)
{
  return tr_cu_r(tr.x,tr.y,tr.z);
}

double DFT::dis(double x1,double y1,double z1,double x2,double y2,double z2)
{
  double dex,dey,dez;
  dex=fabs(x1-x2);
  dey=fabs(y1-y2);
  dez=fabs(z1-z2);
  while(dex>=lx)
  {
    dex=dex-lx;
  }
  while(dey>=ly)
  {
    dey=dey-ly;
  }
  while(dez>=lz)
  {
    dez=dez-lz;
  }
  if(dex>lx/2)
  {
    dex=lx-dex;
  }
  if(dey>ly/2)
  {
    dey=ly-dey;
  }
  if(dez>lz/2)
  {
    dez=lz-dez;
  }
  return sqrt(dex*dex+dey*dey+dez*dez);
}

double DFT::dis_tri_r(double xx,double yy,double zz)
{
  Vector rcu;
  rcu=tr_cu_r(xx,yy,zz);
  return rcu.rabs();
}

double DFT::dis_tri_k(double xx,double yy,double zz)
{
  Vector rcu;
  rcu=tr_cu_k(xx,yy,zz);
  return rcu.rabs();
}

void DFT::cal_ck()
{
  int i,j,k,ig,ig2;
  long nijk;
  double hs,y,wa,x,z,r,rm,wm,wtot;
  double rho,*chem,*rin;
  double xk,yk,zk,rkc;
  Vector kc;
  int irkc;
  double rrkc;
  double urc;

  chem=new double[k_gas];
  rin=new double[k_gas];

  LJ ch(k_gas,sig,be,tem);

  for(i=0;i<k_gas;i++)
  {
	  rin[i]=rb[i];
  }
  ch.r_chem(rin,chem);
  pbulk=ch.press(rin);
  for(i=0;i<k_gas;i++)
  {
	  bmhs[i]=chem[i];
	  cout<<'\n'<<"chemical potential of "<<i<<" : "<<bmhs[i]<<endl;
  }
  
  for(ig=0;ig<k_gas;ig++)
  {
	  for(ig2=0;ig2<=ig;ig2++)
	  {
		  urc=(pow(cutoff,-12)-pow(cutoff,-6))*4*bemix[ig][ig2];
		  for(i=0;i<nx;i++)
		  {
			  if(i<nx/2)
			  {
				  x=i*dx;
				  xk=i*dkx;
			  }
			  else
			  {
				  x=(i-nx)*dx;
				  xk=(i-nx)*dkx;
			  }
			  for(j=0;j<ny;j++)
			  {
				  if(j<ny/2)
				  {				  
					  y=(j)*dy;
					  yk=j*dky;
				  }
				  else
				  {
					  y=(j-ny)*dy;
					  yk=(j-ny)*dky;
				  }
				  for(k=0;k<nz;k++)
				  {
					  nijk=k+nz*(j+ny*i);
					  if(k<nz/2)
					  {
						  z=(k)*dz;
						  zk=k*dkz;
					  }
					  else
					  {
						  z=(k-nz)*dz;
						  zk=(k-nz)*dkz;
					  }
					  if(ig2==0)
					  {
						  kc=tr_cu_k(xk,yk,zk);
						  rkc=kc.rabs();
						  irkc=int(rkc/wfun[ig].dw);
						  rrkc=rkc/wfun[ig].dw-irkc;
						  if(irkc>=4999)
						  {
							  tw0[i][j][k][ig]=0.0;
							  tw3[i][j][k][ig]=0.0;					  
					  //tw1[i][j][k][ig]=0.0;
					  //tw2[i][j][k][ig]=0.0;
					  //tw1x[i][j][k][ig]=0.0;
					  //tw1y[i][j][k][ig]=0.0;
					  //tw1z[i][j][k][ig]=0.0;
							  tw2x[i][j][k][ig]=0.0;
							  tw2y[i][j][k][ig]=0.0;
							  tw2z[i][j][k][ig]=0.0;
							  twcor[i][j][k][ig]=0.0;
						  }
						  else
						  {
							  tw0[i][j][k][ig]=((1.0-rrkc)*wfun[ig].fw0[irkc]+rrkc*wfun[ig].fw0[irkc+1])/nxyz;
					   //tw1[i][j][k][ig]=tw0[i][j][k][ig]*0.5*dbh[ig];
					  //tw2[i][j][k][ig]=tw0[i][j][k][ig]*PI*dbh[ig]*dbh[ig];
							  tw3[i][j][k][ig]=((1.0-rrkc)*wfun[ig].fw3[irkc]+rrkc*wfun[ig].fw3[irkc+1])/nxyz;
							  tw2x[i][j][k][ig]=-tw3[i][j][k][ig]*kc.x;
							  tw2y[i][j][k][ig]=-tw3[i][j][k][ig]*kc.y;
							  tw2z[i][j][k][ig]=-tw3[i][j][k][ig]*kc.z;
					  //tw1x[i][j][k][ig]=tw2x[i][j][k][ig]/(2*PI*dbh[ig]);
					  //tw1y[i][j][k][ig]=tw2y[i][j][k][ig]/(2*PI*dbh[ig]);
					  //tw1z[i][j][k][ig]=tw2z[i][j][k][ig]/(2*PI*dbh[ig]);
							  twcor[i][j][k][ig]=((1.0-rrkc)*wfun[ig].fwcor[irkc]+rrkc*wfun[ig].fwcor[irkc+1])/nxyz;
						  }
					  }
					  r=dis_tri_r(x,y,z)/sigmix[ig][ig2];//(sig[ig]+sig[ig2]);
					  if(r<1.0)  //r<sig*2^(1/6)
					  {
						  wa=0.0;
					  }
					  else if(r>cutoff)
					  {
						  wa=0.0;
					  }
					  else
					  {
						  wa=(pow(r,-12)-pow(r,-6))*4*bemix[ig][ig2]-urc;
					  }
					  fft_in[nijk][0]=wa;
					  fft_in[nijk][1]=0.0;
				  }
			  }
		  }
		  fftw_execute(fft_for);
		  for(i=0;i<nx;i++)
		  {
			  for(j=0;j<ny;j++)
			  {
				  for(k=0;k<nz;k++)
				  {
					  nijk=k+nz*(j+ny*i);
					  uattk[i][j][k][ig][ig2]=fft_out[nijk][0]*dx*dy*dz/nxyz;
				  }
			  }
		  }
		  ////////////////////////////////
      }
  }
  return;
}


double DFT::fdf(double *sqrou,long nnn)
{
  int i,j,k,ir,idmx,idmy,idmz,ig,ig2,ig2t,igt;
  long nijk,nijk2;
  double x,y,z,r;
  double w0,w1,w2,w3,rir,wv2x,wv2y,wv2z,wcor;
  double tn0,tn1,tn2,tn3,tnv2x,tnv2y,tnv2z,tnv1x,tnv1y,tnv1z;
  double t2n3,t3n3,tn30,tn31,tn32,tn33,tnv22,tnv21;
  double rcor,fa;
  double lamda,roun,drou,rouavn;
  int it;
  double rtot1,rtot2;
  double tnv2xr,tnv2xi,tnv2yr,tnv2yi,tnv2zr,tnv2zi;
  double tnv1xr,tnv1xi,tnv1yr,tnv1yi,tnv1zr,tnv1zi;
  double tn0i,tn1i,tn2i,tn3i;
  Vector txyz,wxyz,tixyz,wixyz;
  double n3add,acdens,fre,faicor;
  double dbh05,pidbh2,dbh2pai;
  

  for(nijk=0;nijk<nxyz;nijk++)
  {
	  n0[nijk][0]=0.0;
	  n0[nijk][1]=0.0;
	  n1[nijk][0]=0.0;
	  n1[nijk][1]=0.0;
	  n2[nijk][0]=0.0;
	  n2[nijk][1]=0.0;
	  n3[nijk][0]=0.0;
	  n3[nijk][1]=0.0;
	  nv1x[nijk][0]=0.0;
	  nv1x[nijk][1]=0.0;
	  nv1y[nijk][0]=0.0;
	  nv1y[nijk][1]=0.0;
	  nv1z[nijk][0]=0.0;
	  nv1z[nijk][1]=0.0;
	  nv2x[nijk][0]=0.0;
	  nv2x[nijk][1]=0.0;
	  nv2y[nijk][0]=0.0;
	  nv2y[nijk][1]=0.0;
	  nv2z[nijk][0]=0.0;
	  nv2z[nijk][1]=0.0;
	  cor[nijk][0]=0.0;
	  cor[nijk][1]=0.0;	  
  }
  for(ig=0;ig<k_gas;ig++)
  {
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  denskr[i][j][k][ig]=0.0;
				  denski[i][j][k][ig]=0.0;
			  }
		  }
	  }
  }
  for(ig=0;ig<k_gas;ig++)
  {
	  dbh05=dbh[ig]*0.5;
	  pidbh2=PI*dbh[ig]*dbh[ig];
	  dbh2pai=2*PI*dbh[ig];
	  for(nijk=0;nijk<nxyz;nijk++)
	  {
		  nijk2=nijk+ig*nxyz;
		  fft_in[nijk][0]=sqrou[nijk2]*sqrou[nijk2];
		  fft_in[nijk][1]=0.0;		  
      }
	  fftw_execute(fft_for);
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);				  
				  n0[nijk][0]+=fft_out[nijk][0]*tw0[i][j][k][ig];
				  n0[nijk][1]+=fft_out[nijk][1]*tw0[i][j][k][ig];
				  n1[nijk][0]+=fft_out[nijk][0]*tw0[i][j][k][ig]*dbh05;
				  n1[nijk][1]+=fft_out[nijk][1]*tw0[i][j][k][ig]*dbh05;
				  n2[nijk][0]+=fft_out[nijk][0]*tw0[i][j][k][ig]*pidbh2;
				  n2[nijk][1]+=fft_out[nijk][1]*tw0[i][j][k][ig]*pidbh2;
				  n3[nijk][0]+=fft_out[nijk][0]*tw3[i][j][k][ig];
				  n3[nijk][1]+=fft_out[nijk][1]*tw3[i][j][k][ig];
				  nv2x[nijk][0]+=-fft_out[nijk][1]*tw2x[i][j][k][ig];
				  nv2x[nijk][1]+=fft_out[nijk][0]*tw2x[i][j][k][ig];
				  nv2y[nijk][0]+=-fft_out[nijk][1]*tw2y[i][j][k][ig];
				  nv2y[nijk][1]+=fft_out[nijk][0]*tw2y[i][j][k][ig];
				  nv2z[nijk][0]+=-fft_out[nijk][1]*tw2z[i][j][k][ig];
				  nv2z[nijk][1]+=fft_out[nijk][0]*tw2z[i][j][k][ig];
				  nv1x[nijk][0]+=-fft_out[nijk][1]*tw2x[i][j][k][ig]/dbh2pai;
				  nv1x[nijk][1]+=fft_out[nijk][0]*tw2x[i][j][k][ig]/dbh2pai;
				  nv1y[nijk][0]+=-fft_out[nijk][1]*tw2y[i][j][k][ig]/dbh2pai;
				  nv1y[nijk][1]+=fft_out[nijk][0]*tw2y[i][j][k][ig]/dbh2pai;
				  nv1z[nijk][0]+=-fft_out[nijk][1]*tw2z[i][j][k][ig]/dbh2pai;
				  nv1z[nijk][1]+=fft_out[nijk][0]*tw2z[i][j][k][ig]/dbh2pai;
				  //cor[nijk][0]=fft_out[nijk][0]*twcor[i][j][k];
				  //cor[nijk][1]=fft_out[nijk][1]*twcor[i][j][k];
				  for(ig2=0;ig2<k_gas;ig2++)
				  {
					  if(ig2>ig)
					  {
						  igt=ig2;
						  ig2t=ig;
					  }
					  else
					  {
						  igt=ig;
						  ig2t=ig2;
					  }
					  denskr[i][j][k][ig2]+=fft_out[nijk][0]*uattk[i][j][k][igt][ig2t];
					  denski[i][j][k][ig2]+=fft_out[nijk][1]*uattk[i][j][k][igt][ig2t];
				  }
			  }
		  }
	  }
  }
  fftw_execute(fft_n0_bak);
  fftw_execute(fft_n3_bak);
  fftw_execute(fft_n1_bak);
  fftw_execute(fft_n2_bak);
  fftw_execute(fft_nv1x_bak);
  fftw_execute(fft_nv1y_bak);
  fftw_execute(fft_nv1z_bak);
  fftw_execute(fft_nv2x_bak);
  fftw_execute(fft_nv2y_bak);
  fftw_execute(fft_nv2z_bak);

  for(ig=0;ig<k_gas;ig++)
  {
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);
				  fft_in[nijk][0]=denskr[i][j][k][ig];
				  fft_in[nijk][1]=denski[i][j][k][ig];
			  }
		  }
	  }
	  fftw_execute(fft_bak);
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);
				  denskr[i][j][k][ig]=fft_out[nijk][0];
				  denski[i][j][k][ig]=fft_out[nijk][1];//rou convolution with uatt
			  }
		  }
	  }
  }


  fre=0.0;
  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
      for(k=0;k<nz;k++)
      {
	nijk=k+nz*(j+ny*i);

	for(ig=0;ig<k_gas;ig++)
	{
		nijk2=nijk+ig*nxyz;
		acdens=sqrou[nijk2]*sqrou[nijk2];
		if(acdens>1.0e-5*rb[ig])
		{
			fre+=acdens*(log(acdens)-1.0);//idea gas term
		}
		fre+=(vext[i][j][k][ig]-bmhs[ig])*acdens;//external term
		fre+=0.5*acdens*denskr[i][j][k][ig];//attractive term
	}

	tn0=n0[nijk][0];//nxyz;	
	tn1=n1[nijk][0];
	tn2=n2[nijk][0];
	tn3=n3[nijk][0];//nxyz;
	
	if(tn0<0.0)
	{
	  tn0=0.0;
	  if(tn0<-0.0001)ercode=1;
	}
	if(tn3<0.0)
	{
	  tn3=0.0;
	  if(tn3<-0.0001)ercode=1;
	}
	
	tnv2x=nv2x[nijk][0];//nxyz;
	tnv2y=nv2y[nijk][0];//nxyz;
	tnv2z=nv2z[nijk][0];//nxyz;
	tnv1x=nv1x[nijk][0];//nxyz;
	tnv1y=nv1y[nijk][0];//nxyz;
	tnv1z=nv1z[nijk][0];//nxyz;	

	n3add=tn3-0.7;
	n3add=80.0*(exp(n3add*80.0));
	fre+=exp(80.0*(tn3-0.7));

	if(tn3>fix_cut)
	{
	  tn3=fix_cut;
	  ercode=2;
	}
	else
	{
	  //n3add=0.0;
	}
	
	txyz=tr_cu_r(tnv2x,tnv2y,tnv2z);	
	
	
	if(tn3<1.0e-5)
	{
	  n0[nijk][0]=0.0;
	  n1[nijk][0]=0.0;
	  n2[nijk][0]=0.0;
	  n3[nijk][0]=0.0;
	  nv1x[nijk][0]=0.0;
	  nv1y[nijk][0]=0.0;
	  nv1z[nijk][0]=0.0;
	  nv2x[nijk][0]=0.0;
	  nv2y[nijk][0]=0.0;
	  nv2z[nijk][0]=0.0;
	  n0[nijk][1]=0.0;
	  n1[nijk][1]=0.0;
	  n2[nijk][1]=0.0;
	  n3[nijk][1]=0.0;
	  nv1x[nijk][1]=0.0;
	  nv1y[nijk][1]=0.0;
	  nv1z[nijk][1]=0.0;
	  nv2x[nijk][1]=0.0;
	  nv2y[nijk][1]=0.0;
	  nv2z[nijk][1]=0.0;
	  continue;
	}
	
	
	
	tn31=1.0-tn3;
	tn30=log(tn31);//log(1-n3)
	tn31=1.0/tn31;//1/(1-n3)
	tn32=tn31*tn31;//1/(1-n3)^2
	tn33=tn31*tn32;//1/(1-n3)^3
	tnv22=tnv2x*tnv2x+tnv2y*tnv2y+tnv2z*tnv2z;//nv2*nv2
	tnv21=tnv2x*tnv1x+tnv2y*tnv1y+tnv2z*tnv1z;//nv2*nv1
	t2n3=1.0/tn3/tn3;
	t3n3=t2n3/tn3;
	
	n0[nijk][0]=-tn30;
	n1[nijk][0]=tn2*tn31;
	n2[nijk][0]=(tn1*tn31+(tn30/tn3+tn32)*(tn2*tn2-tnv22)/(12*PI*tn3));
	n0[nijk][1]=0.0;
	n1[nijk][1]=0.0;
	n2[nijk][1]=0.0;
	
	n3[nijk][0]=tn30*t3n3/(18*PI);
	n3[nijk][0]+=(1.0-3*tn3+1.0/tn32)*t2n3*tn33/(36*PI);
	n3[nijk][0]*=3*tn2*tnv22-tn2*tn2*tn2;
	n3[nijk][0]+=tn0*tn31;
	n3[nijk][0]+=(tn1*tn2-tnv21)*tn32;
	n3[nijk][0]+=n3add;
	n3[nijk][1]=0.0;
	
	nv1x[nijk][0]=-tnv2x*tn31;
	nv1x[nijk][0]=0.0;
	nv2x[nijk][0]=-(tn30/tn3+tn32)*tn2*tnv2x/(6*PI*tn3)-tnv1x*tn31;
	nv2x[nijk][1]=0.0;
	
	nv1y[nijk][0]=-tnv2y*tn31;
	nv1y[nijk][1]=0.0;
	nv2y[nijk][0]=-(tn30/tn3+tn32)*tn2*tnv2y/(6*PI*tn3)-tnv1y*tn31;
	nv2y[nijk][1]=0.0;
	
	nv1z[nijk][0]=-tnv2z*tn31;
	nv1z[nijk][0]=0.0;
	nv2z[nijk][0]=-(tn30/tn3+tn32)*tn2*tnv2z/(6*PI*tn3)-tnv1z*tn31;
	nv2z[nijk][1]=0.0;

	fre+=-tn0*tn30;
	fre+=(tn1*tn2-tnv21)*tn31;
	fre+=(pow(tn2,3)-3*tn2*tnv22)*(tn3+tn30/tn32)/(36*PI*tn3*tn3)*tn32;
	

      }
    }
  }
  fre*=dv;

  return fre;
}





void DFT::fdf_dd(double *df,double *sqrou,long nnn)
{
  int i,j,k,ir,idmx,idmy,idmz,ig,ig2,ig2t,igt;
  long nijk,nijk2;
  double x,y,z,r;
  double w0,w1,w2,w3,rir,wv2x,wv2y,wv2z,wcor;
  double tn0,tn1,tn2,tn3,tnv2x,tnv2y,tnv2z,tnv1x,tnv1y,tnv1z;
  double t2n3,t3n3,tn30,tn31,tn32,tn33,tnv22,tnv21;
  double rcor,fa;
  double lamda,roun,drou,rouavn;
  int it;
  double rtot1,rtot2;
  double tnv2xr,tnv2xi,tnv2yr,tnv2yi,tnv2zr,tnv2zi;
  double tnv1xr,tnv1xi,tnv1yr,tnv1yi,tnv1zr,tnv1zi;
  double tn0i,tn1i,tn2i,tn3i;
  Vector txyz,wxyz,tixyz,wixyz;
  double n3add,acdens,fre,faicor;
  double dbh05,pidbh2,dbh2pai;
  double fres;
  

  for(nijk=0;nijk<nxyz;nijk++)
  {
	  n0[nijk][0]=0.0;
	  n0[nijk][1]=0.0;
	  n1[nijk][0]=0.0;
	  n1[nijk][1]=0.0;
	  n2[nijk][0]=0.0;
	  n2[nijk][1]=0.0;
	  n3[nijk][0]=0.0;
	  n3[nijk][1]=0.0;
	  nv1x[nijk][0]=0.0;
	  nv1x[nijk][1]=0.0;
	  nv1y[nijk][0]=0.0;
	  nv1y[nijk][1]=0.0;
	  nv1z[nijk][0]=0.0;
	  nv1z[nijk][1]=0.0;
	  nv2x[nijk][0]=0.0;
	  nv2x[nijk][1]=0.0;
	  nv2y[nijk][0]=0.0;
	  nv2y[nijk][1]=0.0;
	  nv2z[nijk][0]=0.0;
	  nv2z[nijk][1]=0.0;
	  cor[nijk][0]=0.0;
	  cor[nijk][1]=0.0;	  
  }
  for(ig=0;ig<k_gas;ig++)
  {
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  denskr[i][j][k][ig]=0.0;
				  denski[i][j][k][ig]=0.0;
			  }
		  }
	  }
  }
  for(ig=0;ig<k_gas;ig++)
  {
	  dbh05=dbh[ig]*0.5;
	  pidbh2=PI*dbh[ig]*dbh[ig];
	  dbh2pai=2*PI*dbh[ig];
	  for(nijk=0;nijk<nxyz;nijk++)
	  {
		  nijk2=nijk+ig*nxyz;
		  fft_in[nijk][0]=sqrou[nijk2]*sqrou[nijk2];
		  fft_in[nijk][1]=0.0;		  
      }
	  fftw_execute(fft_for);
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);				  
				  n0[nijk][0]+=fft_out[nijk][0]*tw0[i][j][k][ig];
				  n0[nijk][1]+=fft_out[nijk][1]*tw0[i][j][k][ig];
				  n1[nijk][0]+=fft_out[nijk][0]*tw0[i][j][k][ig]*dbh05;
				  n1[nijk][1]+=fft_out[nijk][1]*tw0[i][j][k][ig]*dbh05;
				  n2[nijk][0]+=fft_out[nijk][0]*tw0[i][j][k][ig]*pidbh2;
				  n2[nijk][1]+=fft_out[nijk][1]*tw0[i][j][k][ig]*pidbh2;
				  n3[nijk][0]+=fft_out[nijk][0]*tw3[i][j][k][ig];
				  n3[nijk][1]+=fft_out[nijk][1]*tw3[i][j][k][ig];
				  nv2x[nijk][0]+=-fft_out[nijk][1]*tw2x[i][j][k][ig];
				  nv2x[nijk][1]+=fft_out[nijk][0]*tw2x[i][j][k][ig];
				  nv2y[nijk][0]+=-fft_out[nijk][1]*tw2y[i][j][k][ig];
				  nv2y[nijk][1]+=fft_out[nijk][0]*tw2y[i][j][k][ig];
				  nv2z[nijk][0]+=-fft_out[nijk][1]*tw2z[i][j][k][ig];
				  nv2z[nijk][1]+=fft_out[nijk][0]*tw2z[i][j][k][ig];
				  nv1x[nijk][0]+=-fft_out[nijk][1]*tw2x[i][j][k][ig]/dbh2pai;
				  nv1x[nijk][1]+=fft_out[nijk][0]*tw2x[i][j][k][ig]/dbh2pai;
				  nv1y[nijk][0]+=-fft_out[nijk][1]*tw2y[i][j][k][ig]/dbh2pai;
				  nv1y[nijk][1]+=fft_out[nijk][0]*tw2y[i][j][k][ig]/dbh2pai;
				  nv1z[nijk][0]+=-fft_out[nijk][1]*tw2z[i][j][k][ig]/dbh2pai;
				  nv1z[nijk][1]+=fft_out[nijk][0]*tw2z[i][j][k][ig]/dbh2pai;
				  //cor[nijk][0]=fft_out[nijk][0]*twcor[i][j][k];
				  //cor[nijk][1]=fft_out[nijk][1]*twcor[i][j][k];
				  for(ig2=0;ig2<k_gas;ig2++)
				  {
					  if(ig2>ig)
					  {
						  igt=ig2;
						  ig2t=ig;
					  }
					  else
					  {
						  igt=ig;
						  ig2t=ig2;
					  }
					  denskr[i][j][k][ig2]+=fft_out[nijk][0]*uattk[i][j][k][igt][ig2t];
					  denski[i][j][k][ig2]+=fft_out[nijk][1]*uattk[i][j][k][igt][ig2t];
				  }
			  }
		  }
	  }
  }
  fftw_execute(fft_n0_bak);
  fftw_execute(fft_n3_bak);
  fftw_execute(fft_n1_bak);
  fftw_execute(fft_n2_bak);
  fftw_execute(fft_nv1x_bak);
  fftw_execute(fft_nv1y_bak);
  fftw_execute(fft_nv1z_bak);
  fftw_execute(fft_nv2x_bak);
  fftw_execute(fft_nv2y_bak);
  fftw_execute(fft_nv2z_bak);

  for(ig=0;ig<k_gas;ig++)
  {
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);
				  fft_in[nijk][0]=denskr[i][j][k][ig];
				  fft_in[nijk][1]=denski[i][j][k][ig];
			  }
		  }
	  }
	  fftw_execute(fft_bak);
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);
				  denskr[i][j][k][ig]=fft_out[nijk][0];
				  denski[i][j][k][ig]=fft_out[nijk][1];//rou convolution with uatt
			  }
		  }
	  }
  }


  fre=0.0;
  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
      for(k=0;k<nz;k++)
      {
	nijk=k+nz*(j+ny*i);
	fres=0.0;
	for(ig=0;ig<k_gas;ig++)
	{
		nijk2=nijk+ig*nxyz;
		acdens=sqrou[nijk2]*sqrou[nijk2];
		if(acdens>1.0e-5*rb[ig])
		{
			fres+=acdens*(log(acdens)-1.0);//idea gas term
		}
		fres+=(vext[i][j][k][ig]-bmhs[ig])*acdens;//external term
		fre+=0.5*acdens*denskr[i][j][k][ig];//attractive term
	}

	tn0=n0[nijk][0];//nxyz;	
	tn1=n1[nijk][0];
	tn2=n2[nijk][0];
	tn3=n3[nijk][0];//nxyz;
	
	if(tn0<0.0)
	{
	  tn0=0.0;
	  if(tn0<-0.0001)ercode=1;
	}
	if(tn3<0.0)
	{
	  tn3=0.0;
	  if(tn3<-0.0001)ercode=1;
	}
	
	tnv2x=nv2x[nijk][0];//nxyz;
	tnv2y=nv2y[nijk][0];//nxyz;
	tnv2z=nv2z[nijk][0];//nxyz;
	tnv1x=nv1x[nijk][0];//nxyz;
	tnv1y=nv1y[nijk][0];//nxyz;
	tnv1z=nv1z[nijk][0];//nxyz;	
	
	fres+=exp(80.0*(tn3-0.7));
	n3add=tn3-0.7;
	n3add=80.0*(exp(n3add*80.0));

	if(tn3>fix_cut)
	{
	  tn3=fix_cut;
	  ercode=2;
	}
	else
	{
	  //n3add=0.0;
	}
	
	txyz=tr_cu_r(tnv2x,tnv2y,tnv2z);	
	
	
	if(tn3<1.0e-5)
	{
	  n0[nijk][0]=0.0;
	  n1[nijk][0]=0.0;
	  n2[nijk][0]=0.0;
	  n3[nijk][0]=0.0;
	  nv1x[nijk][0]=0.0;
	  nv1y[nijk][0]=0.0;
	  nv1z[nijk][0]=0.0;
	  nv2x[nijk][0]=0.0;
	  nv2y[nijk][0]=0.0;
	  nv2z[nijk][0]=0.0;
	  n0[nijk][1]=0.0;
	  n1[nijk][1]=0.0;
	  n2[nijk][1]=0.0;
	  n3[nijk][1]=0.0;
	  nv1x[nijk][1]=0.0;
	  nv1y[nijk][1]=0.0;
	  nv1z[nijk][1]=0.0;
	  nv2x[nijk][1]=0.0;
	  nv2y[nijk][1]=0.0;
	  nv2z[nijk][1]=0.0;
	  continue;
	}
	
	
	
	tn31=1.0-tn3;
	tn30=log(tn31);//log(1-n3)
	tn31=1.0/tn31;//1/(1-n3)
	tn32=tn31*tn31;//1/(1-n3)^2
	tn33=tn31*tn32;//1/(1-n3)^3
	tnv22=tnv2x*tnv2x+tnv2y*tnv2y+tnv2z*tnv2z;//nv2*nv2
	tnv21=tnv2x*tnv1x+tnv2y*tnv1y+tnv2z*tnv1z;//nv2*nv1
	t2n3=1.0/tn3/tn3;
	t3n3=t2n3/tn3;
	

	fres+=-tn0*tn30;
	fres+=(tn1*tn2-tnv21)*tn31;
	fres+=(pow(tn2,3)-3*tn2*tnv22)*(tn3+tn30/tn32)/(36*PI*tn3*tn3)*tn32;

	

	if(fres<-100.0)
	{
		//cout<<"something may be wrong"<<endl;
	}


	n0[nijk][0]=-tn30;
	n1[nijk][0]=tn2*tn31;
	n2[nijk][0]=(tn1*tn31+(tn30/tn3+tn32)*(tn2*tn2-tnv22)/(12*PI*tn3));
	n0[nijk][1]=0.0;
	n1[nijk][1]=0.0;
	n2[nijk][1]=0.0;
	
	n3[nijk][0]=tn30*t3n3/(18*PI);
	n3[nijk][0]+=(1.0-3*tn3+1.0/tn32)*t2n3*tn33/(36*PI);
	n3[nijk][0]*=3*tn2*tnv22-tn2*tn2*tn2;
	n3[nijk][0]+=tn0*tn31;
	n3[nijk][0]+=(tn1*tn2-tnv21)*tn32;
	
	
	n3[nijk][0]+=n3add;
	n3[nijk][1]=0.0;
	
	nv1x[nijk][0]=-tnv2x*tn31;
	nv1x[nijk][0]=0.0;
	nv2x[nijk][0]=-(tn30/tn3+tn32)*tn2*tnv2x/(6*PI*tn3)-tnv1x*tn31;
	nv2x[nijk][1]=0.0;
	
	nv1y[nijk][0]=-tnv2y*tn31;
	nv1y[nijk][1]=0.0;
	nv2y[nijk][0]=-(tn30/tn3+tn32)*tn2*tnv2y/(6*PI*tn3)-tnv1y*tn31;
	nv2y[nijk][1]=0.0;
	
	nv1z[nijk][0]=-tnv2z*tn31;
	nv1z[nijk][0]=0.0;
	nv2z[nijk][0]=-(tn30/tn3+tn32)*tn2*tnv2z/(6*PI*tn3)-tnv1z*tn31;
	nv2z[nijk][1]=0.0;



	fre+=fres;
      }
    }
  }
  fre*=dv;

  fftw_execute(fft_n0_for);
  fftw_execute(fft_n3_for);
  fftw_execute(fft_n1_for);
  fftw_execute(fft_n2_for);
  fftw_execute(fft_nv1x_for);
  fftw_execute(fft_nv1y_for);
  fftw_execute(fft_nv1z_for);
  fftw_execute(fft_nv2x_for);
  fftw_execute(fft_nv2y_for);
  fftw_execute(fft_nv2z_for);

  srav=0.0;
  srmax=0.0;
  rou_tot=0.0;

  for(ig=0;ig<k_gas;ig++)
  {
	  dbh05=dbh[ig]*0.5;
	  pidbh2=PI*dbh[ig]*dbh[ig];
	  dbh2pai=2*PI*dbh[ig];
	  rouav[ig]=0.0;

	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);

	fft_in[nijk][0]=n0[nijk][0]*tw0[i][j][k][ig];
	fft_in[nijk][1]=n0[nijk][1]*tw0[i][j][k][ig];//n0 term

	fft_in[nijk][0]+=n1[nijk][0]*tw0[i][j][k][ig]*dbh05;
	fft_in[nijk][1]+=n1[nijk][1]*tw0[i][j][k][ig]*dbh05;//n1

	fft_in[nijk][0]+=n2[nijk][0]*tw0[i][j][k][ig]*pidbh2;
	fft_in[nijk][1]+=n2[nijk][1]*tw0[i][j][k][ig]*pidbh2;//n2

	fft_in[nijk][0]+=n3[nijk][0]*tw3[i][j][k][ig];
	fft_in[nijk][1]+=n3[nijk][1]*tw3[i][j][k][ig];//n3

	fft_in[nijk][0]+=nv1x[nijk][1]*tw2x[i][j][k][ig]/dbh2pai;
	fft_in[nijk][1]+=-nv1x[nijk][0]*tw2x[i][j][k][ig]/dbh2pai;//nv1x

	fft_in[nijk][0]+=nv1y[nijk][1]*tw2y[i][j][k][ig]/dbh2pai;
	fft_in[nijk][1]+=-nv1y[nijk][0]*tw2y[i][j][k][ig]/dbh2pai;//nv1y

	fft_in[nijk][0]+=nv1z[nijk][1]*tw2z[i][j][k][ig]/dbh2pai;
	fft_in[nijk][1]+=-nv1z[nijk][0]*tw2z[i][j][k][ig]/dbh2pai;//nv1z

	fft_in[nijk][0]+=nv2x[nijk][1]*tw2x[i][j][k][ig];
	fft_in[nijk][1]+=-nv2x[nijk][0]*tw2x[i][j][k][ig];//nv2x

	fft_in[nijk][0]+=nv2y[nijk][1]*tw2y[i][j][k][ig];
	fft_in[nijk][1]+=-nv2y[nijk][0]*tw2y[i][j][k][ig];//nv2y

	fft_in[nijk][0]+=nv2z[nijk][1]*tw2z[i][j][k][ig];
	fft_in[nijk][1]+=-nv2z[nijk][0]*tw2z[i][j][k][ig];//nv2z
	
			  }
		  }
	  }

	  fftw_execute(fft_bak);

	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);
				  nijk2=nijk+ig*nxyz;
				  acdens=sqrou[nijk2]*sqrou[nijk2];
				  if(acdens>1.0E-5*rb[ig])
				  {
					  df[nijk2]=log(acdens)+vext[i][j][k][ig]-bmhs[ig]+fft_out[nijk][0]+denskr[i][j][k][ig];	
				  }
				  else
				  {
					  df[nijk2]=0.0;
				  }
				  df[nijk2]*=2*sqrou[nijk2];

				  //converge test//
				  if(vext[i][j][k][ig]<200.0)
				  {
					  roun=exp(bmhs[ig]-fft_out[nijk][0]-denskr[i][j][k][ig]-vext[i][j][k][ig]);
					  //roun=rb*exp(-lamd[i][j][k]-vext[i][j][k]);
				  }
				  else
				  {
					  roun=0.0;
				  }
				  drou=fabs(roun-acdens);
				  srav+=drou*drou;
				  if(drou>srmax)
				  {
					  srmax=drou;
					  idmx=i;
					  idmy=j;
					  idmz=k;
				  }
				  rouav[ig]+=acdens;
			  }
		  }
	  }
	  rouav[ig]=rouav[ig]/nxyz;	  
	  rou_tot+=rouav[ig];
  }
  srav=sqrt(srav/nxyz)/rou_tot;
  srmax=srmax/rou_tot;

  /*cout<<fre<<" "<<srav<<" "<<srmax;
  for(ig=0;ig<k_gas;ig++)
  {
	  cout<<" "<<rouav[ig];
  }
  cout<<endl;*/

  return ;
}













double DFT::fdf_t(double *df,double *sqrou,long nnn)
{
  int i,j,k,ir,idmx,idmy,idmz,ig,ig2,ig2t,igt;
  long nijk,nijk2;
  double x,y,z,r;
  double w0,w1,w2,w3,rir,wv2x,wv2y,wv2z,wcor;
  double tn0,tn1,tn2,tn3,tnv2x,tnv2y,tnv2z,tnv1x,tnv1y,tnv1z;
  double t2n3,t3n3,tn30,tn31,tn32,tn33,tnv22,tnv21;
  double rcor,fa;
  double lamda,roun,drou,rouavn;
  int it;
  double rtot1,rtot2;
  double tnv2xr,tnv2xi,tnv2yr,tnv2yi,tnv2zr,tnv2zi;
  double tnv1xr,tnv1xi,tnv1yr,tnv1yi,tnv1zr,tnv1zi;
  double tn0i,tn1i,tn2i,tn3i;
  Vector txyz,wxyz,tixyz,wixyz;
  double n3add,acdens,fre,faicor;
  double dbh05,pidbh2,dbh2pai;
  

  for(nijk=0;nijk<nxyz;nijk++)
  {
	  n0[nijk][0]=0.0;
	  n0[nijk][1]=0.0;
	  n1[nijk][0]=0.0;
	  n1[nijk][1]=0.0;
	  n2[nijk][0]=0.0;
	  n2[nijk][1]=0.0;
	  n3[nijk][0]=0.0;
	  n3[nijk][1]=0.0;
	  nv1x[nijk][0]=0.0;
	  nv1x[nijk][1]=0.0;
	  nv1y[nijk][0]=0.0;
	  nv1y[nijk][1]=0.0;
	  nv1z[nijk][0]=0.0;
	  nv1z[nijk][1]=0.0;
	  nv2x[nijk][0]=0.0;
	  nv2x[nijk][1]=0.0;
	  nv2y[nijk][0]=0.0;
	  nv2y[nijk][1]=0.0;
	  nv2z[nijk][0]=0.0;
	  nv2z[nijk][1]=0.0;
	  cor[nijk][0]=0.0;
	  cor[nijk][1]=0.0;	  
  }
  for(ig=0;ig<k_gas;ig++)
  {
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  denskr[i][j][k][ig]=0.0;
				  denski[i][j][k][ig]=0.0;
			  }
		  }
	  }
  }
  for(ig=0;ig<k_gas;ig++)
  {
	  dbh05=dbh[ig]*0.5;
	  pidbh2=PI*dbh[ig]*dbh[ig];
	  dbh2pai=2*PI*dbh[ig];
	  for(nijk=0;nijk<nxyz;nijk++)
	  {
		  nijk2=nijk+ig*nxyz;
		  fft_in[nijk][0]=sqrou[nijk2]*sqrou[nijk2];
		  fft_in[nijk][1]=0.0;		  
      }
	  fftw_execute(fft_for);
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);				  
				  n0[nijk][0]+=fft_out[nijk][0]*tw0[i][j][k][ig];
				  n0[nijk][1]+=fft_out[nijk][1]*tw0[i][j][k][ig];
				  n1[nijk][0]+=fft_out[nijk][0]*tw0[i][j][k][ig]*dbh05;
				  n1[nijk][1]+=fft_out[nijk][1]*tw0[i][j][k][ig]*dbh05;
				  n2[nijk][0]+=fft_out[nijk][0]*tw0[i][j][k][ig]*pidbh2;
				  n2[nijk][1]+=fft_out[nijk][1]*tw0[i][j][k][ig]*pidbh2;
				  n3[nijk][0]+=fft_out[nijk][0]*tw3[i][j][k][ig];
				  n3[nijk][1]+=fft_out[nijk][1]*tw3[i][j][k][ig];
				  nv2x[nijk][0]+=-fft_out[nijk][1]*tw2x[i][j][k][ig];
				  nv2x[nijk][1]+=fft_out[nijk][0]*tw2x[i][j][k][ig];
				  nv2y[nijk][0]+=-fft_out[nijk][1]*tw2y[i][j][k][ig];
				  nv2y[nijk][1]+=fft_out[nijk][0]*tw2y[i][j][k][ig];
				  nv2z[nijk][0]+=-fft_out[nijk][1]*tw2z[i][j][k][ig];
				  nv2z[nijk][1]+=fft_out[nijk][0]*tw2z[i][j][k][ig];
				  nv1x[nijk][0]+=-fft_out[nijk][1]*tw2x[i][j][k][ig]/dbh2pai;
				  nv1x[nijk][1]+=fft_out[nijk][0]*tw2x[i][j][k][ig]/dbh2pai;
				  nv1y[nijk][0]+=-fft_out[nijk][1]*tw2y[i][j][k][ig]/dbh2pai;
				  nv1y[nijk][1]+=fft_out[nijk][0]*tw2y[i][j][k][ig]/dbh2pai;
				  nv1z[nijk][0]+=-fft_out[nijk][1]*tw2z[i][j][k][ig]/dbh2pai;
				  nv1z[nijk][1]+=fft_out[nijk][0]*tw2z[i][j][k][ig]/dbh2pai;
				  //cor[nijk][0]=fft_out[nijk][0]*twcor[i][j][k];
				  //cor[nijk][1]=fft_out[nijk][1]*twcor[i][j][k];
				  for(ig2=0;ig2<k_gas;ig2++)
				  {
					  if(ig2>ig)
					  {
						  igt=ig2;
						  ig2t=ig;
					  }
					  else
					  {
						  igt=ig;
						  ig2t=ig2;
					  }
					  denskr[i][j][k][ig2]+=fft_out[nijk][0]*uattk[i][j][k][igt][ig2t];
					  denski[i][j][k][ig2]+=fft_out[nijk][1]*uattk[i][j][k][igt][ig2t];
				  }
			  }
		  }
	  }
  }
  fftw_execute(fft_n0_bak);
  fftw_execute(fft_n3_bak);
  fftw_execute(fft_n1_bak);
  fftw_execute(fft_n2_bak);
  fftw_execute(fft_nv1x_bak);
  fftw_execute(fft_nv1y_bak);
  fftw_execute(fft_nv1z_bak);
  fftw_execute(fft_nv2x_bak);
  fftw_execute(fft_nv2y_bak);
  fftw_execute(fft_nv2z_bak);

  for(ig=0;ig<k_gas;ig++)
  {
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);
				  fft_in[nijk][0]=denskr[i][j][k][ig];
				  fft_in[nijk][1]=denski[i][j][k][ig];
			  }
		  }
	  }
	  fftw_execute(fft_bak);
	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);
				  denskr[i][j][k][ig]=fft_out[nijk][0];
				  denski[i][j][k][ig]=fft_out[nijk][1];//rou convolution with uatt
			  }
		  }
	  }
  }


  fre=0.0;
  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
      for(k=0;k<nz;k++)
      {
	nijk=k+nz*(j+ny*i);

	for(ig=0;ig<k_gas;ig++)
	{
		nijk2=nijk+ig*nxyz;
		acdens=sqrou[nijk2]*sqrou[nijk2];
		if(acdens>1.0e-5*rb[ig])
		{
			fre+=acdens*(log(acdens)-1.0);//idea gas term
		}
		fre+=(vext[i][j][k][ig]-bmhs[ig])*acdens;//external term
		fre+=0.5*acdens*denskr[i][j][k][ig];//attractive term
	}

	tn0=n0[nijk][0];//nxyz;	
	tn1=n1[nijk][0];
	tn2=n2[nijk][0];
	tn3=n3[nijk][0];//nxyz;
	
	if(tn0<0.0)
	{
	  tn0=0.0;
	  if(tn0<-0.0001)ercode=1;
	}
	if(tn3<0.0)
	{
	  tn3=0.0;
	  if(tn3<-0.0001)ercode=1;
	}
	
	tnv2x=nv2x[nijk][0];//nxyz;
	tnv2y=nv2y[nijk][0];//nxyz;
	tnv2z=nv2z[nijk][0];//nxyz;
	tnv1x=nv1x[nijk][0];//nxyz;
	tnv1y=nv1y[nijk][0];//nxyz;
	tnv1z=nv1z[nijk][0];//nxyz;	
	
	fre+=exp(80.0*(tn3-0.7));
	n3add=tn3-0.7;
	n3add=80.0*(exp(n3add*80.0));

	if(tn3>fix_cut)
	{
	  tn3=fix_cut;
	  ercode=2;
	}
	else
	{
	  //n3add=0.0;
	}
	
	txyz=tr_cu_r(tnv2x,tnv2y,tnv2z);	
	
	
	if(tn3<1.0e-5)
	{
	  n0[nijk][0]=0.0;
	  n1[nijk][0]=0.0;
	  n2[nijk][0]=0.0;
	  n3[nijk][0]=0.0;
	  nv1x[nijk][0]=0.0;
	  nv1y[nijk][0]=0.0;
	  nv1z[nijk][0]=0.0;
	  nv2x[nijk][0]=0.0;
	  nv2y[nijk][0]=0.0;
	  nv2z[nijk][0]=0.0;
	  n0[nijk][1]=0.0;
	  n1[nijk][1]=0.0;
	  n2[nijk][1]=0.0;
	  n3[nijk][1]=0.0;
	  nv1x[nijk][1]=0.0;
	  nv1y[nijk][1]=0.0;
	  nv1z[nijk][1]=0.0;
	  nv2x[nijk][1]=0.0;
	  nv2y[nijk][1]=0.0;
	  nv2z[nijk][1]=0.0;
	  continue;
	}
	
	
	
	tn31=1.0-tn3;
	tn30=log(tn31);//log(1-n3)
	tn31=1.0/tn31;//1/(1-n3)
	tn32=tn31*tn31;//1/(1-n3)^2
	tn33=tn31*tn32;//1/(1-n3)^3
	tnv22=tnv2x*tnv2x+tnv2y*tnv2y+tnv2z*tnv2z;//nv2*nv2
	tnv21=tnv2x*tnv1x+tnv2y*tnv1y+tnv2z*tnv1z;//nv2*nv1
	t2n3=1.0/tn3/tn3;
	t3n3=t2n3/tn3;
	
	n0[nijk][0]=-tn30;
	n1[nijk][0]=tn2*tn31;
	n2[nijk][0]=(tn1*tn31+(tn30/tn3+tn32)*(tn2*tn2-tnv22)/(12*PI*tn3));
	n0[nijk][1]=0.0;
	n1[nijk][1]=0.0;
	n2[nijk][1]=0.0;
	
	n3[nijk][0]=tn30*t3n3/(18*PI);
	n3[nijk][0]+=(1.0-3*tn3+1.0/tn32)*t2n3*tn33/(36*PI);
	n3[nijk][0]*=3*tn2*tnv22-tn2*tn2*tn2;
	n3[nijk][0]+=tn0*tn31;
	n3[nijk][0]+=(tn1*tn2-tnv21)*tn32;

	
	n3[nijk][0]+=n3add;
	n3[nijk][1]=0.0;
	
	nv1x[nijk][0]=-tnv2x*tn31;
	nv1x[nijk][0]=0.0;
	nv2x[nijk][0]=-(tn30/tn3+tn32)*tn2*tnv2x/(6*PI*tn3)-tnv1x*tn31;
	nv2x[nijk][1]=0.0;
	
	nv1y[nijk][0]=-tnv2y*tn31;
	nv1y[nijk][1]=0.0;
	nv2y[nijk][0]=-(tn30/tn3+tn32)*tn2*tnv2y/(6*PI*tn3)-tnv1y*tn31;
	nv2y[nijk][1]=0.0;
	
	nv1z[nijk][0]=-tnv2z*tn31;
	nv1z[nijk][0]=0.0;
	nv2z[nijk][0]=-(tn30/tn3+tn32)*tn2*tnv2z/(6*PI*tn3)-tnv1z*tn31;
	nv2z[nijk][1]=0.0;

	fre+=-tn0*tn30;
	fre+=(tn1*tn2-tnv21)*tn31;
	fre+=(pow(tn2,3)-3*tn2*tnv22)*(tn3+tn30/tn32)/(36*PI*tn3*tn3)*tn32;
	

      }
    }
  }
  fre*=dv;

  fftw_execute(fft_n0_for);
  fftw_execute(fft_n3_for);
  fftw_execute(fft_n1_for);
  fftw_execute(fft_n2_for);
  fftw_execute(fft_nv1x_for);
  fftw_execute(fft_nv1y_for);
  fftw_execute(fft_nv1z_for);
  fftw_execute(fft_nv2x_for);
  fftw_execute(fft_nv2y_for);
  fftw_execute(fft_nv2z_for);

  srav=0.0;
  srmax=0.0;
  rou_tot=0.0;

  for(ig=0;ig<k_gas;ig++)
  {
	  dbh05=dbh[ig]*0.5;
	  pidbh2=PI*dbh[ig]*dbh[ig];
	  dbh2pai=2*PI*dbh[ig];
	  rouav[ig]=0.0;

	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);

	fft_in[nijk][0]=n0[nijk][0]*tw0[i][j][k][ig];
	fft_in[nijk][1]=n0[nijk][1]*tw0[i][j][k][ig];//n0 term

	fft_in[nijk][0]+=n1[nijk][0]*tw0[i][j][k][ig]*dbh05;
	fft_in[nijk][1]+=n1[nijk][1]*tw0[i][j][k][ig]*dbh05;//n1

	fft_in[nijk][0]+=n2[nijk][0]*tw0[i][j][k][ig]*pidbh2;
	fft_in[nijk][1]+=n2[nijk][1]*tw0[i][j][k][ig]*pidbh2;//n2

	fft_in[nijk][0]+=n3[nijk][0]*tw3[i][j][k][ig];
	fft_in[nijk][1]+=n3[nijk][1]*tw3[i][j][k][ig];//n3

	fft_in[nijk][0]+=nv1x[nijk][1]*tw2x[i][j][k][ig]/dbh2pai;
	fft_in[nijk][1]+=-nv1x[nijk][0]*tw2x[i][j][k][ig]/dbh2pai;//nv1x

	fft_in[nijk][0]+=nv1y[nijk][1]*tw2y[i][j][k][ig]/dbh2pai;
	fft_in[nijk][1]+=-nv1y[nijk][0]*tw2y[i][j][k][ig]/dbh2pai;//nv1y

	fft_in[nijk][0]+=nv1z[nijk][1]*tw2z[i][j][k][ig]/dbh2pai;
	fft_in[nijk][1]+=-nv1z[nijk][0]*tw2z[i][j][k][ig]/dbh2pai;//nv1z

	fft_in[nijk][0]+=nv2x[nijk][1]*tw2x[i][j][k][ig];
	fft_in[nijk][1]+=-nv2x[nijk][0]*tw2x[i][j][k][ig];//nv2x

	fft_in[nijk][0]+=nv2y[nijk][1]*tw2y[i][j][k][ig];
	fft_in[nijk][1]+=-nv2y[nijk][0]*tw2y[i][j][k][ig];//nv2y

	fft_in[nijk][0]+=nv2z[nijk][1]*tw2z[i][j][k][ig];
	fft_in[nijk][1]+=-nv2z[nijk][0]*tw2z[i][j][k][ig];//nv2z
	
			  }
		  }
	  }

	  fftw_execute(fft_bak);

	  for(i=0;i<nx;i++)
	  {
		  for(j=0;j<ny;j++)
		  {
			  for(k=0;k<nz;k++)
			  {
				  nijk=k+nz*(j+ny*i);
				  nijk2=nijk+ig*nxyz;
				  acdens=sqrou[nijk2]*sqrou[nijk2];
				  if(acdens>1.0E-5*rb[ig])
				  {
					  df[nijk2]=log(acdens)+vext[i][j][k][ig]-bmhs[ig]+fft_out[nijk][0]+denskr[i][j][k][ig];	
				  }
				  else
				  {
					  df[nijk2]=0.0;
				  }
				  df[nijk2]*=2*sqrou[nijk2];

				  //converge test//
				  if(vext[i][j][k][ig]<200.0)
				  {
					  roun=exp(bmhs[ig]-fft_out[nijk][0]-denskr[i][j][k][ig]-vext[i][j][k][ig]);
					  //roun=rb*exp(-lamd[i][j][k]-vext[i][j][k]);
				  }
				  else
				  {
					  roun=0.0;
				  }
				  drou=fabs(roun-acdens);
				  srav+=drou*drou;
				  if(drou>srmax)
				  {
					  srmax=drou;
					  idmx=i;
					  idmy=j;
					  idmz=k;
				  }
				  rouav[ig]+=acdens;
			  }
		  }
	  }
	  rouav[ig]=rouav[ig]/nxyz;	  
	  rou_tot+=rouav[ig];
  }
  srav=sqrt(srav/nxyz)/rou_tot;
  srmax=srmax/rou_tot;

  cout<<fre<<" "<<srav<<" "<<srmax;
  for(ig=0;ig<k_gas;ig++)
  {
	  cout<<" "<<rouav[ig];
  }
  cout<<endl;

  return fre;
}












void DFT::final()
{
  int i,j,k,ig;
  long nijk,nijk2;
  double sum_ds,tds,acdens;
  int ir;
  double wcor,tcor,x,y,z,r,rir,sumcor;
  double tsex,fvol,mco,dbulk,dtot,dstar;
  
  ofstream fop1("density.dat",ios::out);
  ofstream fop2("output.dat",ios::app);
  fop1<<nx<<" "<<ny<<" "<<nz;

  
  for(ig=0;ig<k_gas;ig++)
  {
	  rouav[ig]=0.0;
	  max_dens[ig]=0.0;
	  fop1<<" "<<name_gas[ig];
  }
 

  for(i=0;i<nx;i++)
  {
	  for(j=0;j<ny;j++)
	  {
		  for(k=0;k<nz;k++)
		  {
			  nijk=k+nz*(j+ny*i);
			  fop1<<'\n'<<i*dx<<" "<<j*dy<<" "<<k*dz;
			  for(ig=0;ig<k_gas;ig++)
			  {
				  nijk2=nijk+ig*nxyz;
				  acdens=m_x[nijk2]*m_x[nijk2];
				  rouav[ig]+=acdens;
				  fop1<<" "<<acdens;
				  if(acdens>max_dens[ig])
				  {
					  max_dens[ig]=acdens;
				  }
			  }
		  }
	  }
  }
  fop1.close();

  rou_tot=0.0;
  for(ig=0;ig<k_gas;ig++)
  {
	  rouav[ig]=rouav[ig]/nxyz;
	  rou_tot+=rouav[ig];
  }

  fop2<<'\n'<<"Pressure of bulk: "<< pbulk<<" atm";
  fop2<<'\n'<<"Molecule    Average_density(mol/L)     Max_density/bulk";
  for(ig=0;ig<k_gas;ig++)
  {
	  fop2<<'\n'<<name_gas[ig]<<" "<<rouav[ig]*1.0e4/6.023<<" ("<<rouav[ig]*100/rou_tot<<"%)  "<<max_dens[ig]/rb[ig];
  }
  fop2<<'\n'<<" total density: "<<rou_tot;
  fop2.close();

  return;
}

int DFT::run()
{
  /*ercode=-1;
  cal_ck();
  //cout<<"Start iteration"<<endl;

  for(step=1;step < 1000;step++)
  {
    ercode=0;
    cal_dens();
    //if(step==10)break;

    cout<<step<<" "<<srav<<" "<<srmax<<" "<<rouav<<" "<<ercode<<endl;
    if(fabs(srmax)<torr)
    {
      break;
    }
  }
  final();*/
  return 0;
}


long DFT::tot_grids()
{
  return nx*ny*nz;
}

long DFT::tot_sites()
{
  return nx*ny*nz*k_gas;
}

double DFT::time_cost()
{
  return (double)(clock()-start)/CLOCKS_PER_SEC;;
}


//DFT mof;
DFT* MINIMIZATION;

double myvalue(double* xrho,long num)
{
	return MINIMIZATION->fdf(xrho, num);
}

void mygrad(double* grad, double* xrho, long num)
{
	return MINIMIZATION->fdf_dd(grad,xrho,num);
}

double myvalgrad(double* grad, double* xrho, long num)
{
	return MINIMIZATION->fdf_t(grad,xrho,num);
}




