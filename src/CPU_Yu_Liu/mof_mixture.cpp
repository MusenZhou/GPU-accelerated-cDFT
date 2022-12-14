#include "dft.h"
#include <iostream>
using namespace std;
int main()
{
	DFT mof;
	cout<<"initial OK"<<endl;
	MINIMIZATION=&mof;
	
	MINIMIZATION->cal_ck();
	cout<<"start minimization"<<endl;
	cg_descent (MINIMIZATION->m_x, MINIMIZATION->tot_sites(), NULL, NULL, 1.e-9, myvalue, mygrad, myvalgrad, NULL) ;
	MINIMIZATION->final();
	cout<<"time cost:  "<<MINIMIZATION->time_cost()<<endl;
	MINIMIZATION->~DFT();
	return 0;
}

