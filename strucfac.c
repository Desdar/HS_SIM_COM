#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include<sstream>
#include<fstream>
#include <time.h>
#include<md_sim.h>
//#include "nrutil.h"
//#include"ini.hpp "   
using namespace std;


double S_q(Particle *p, int NMOL, Cont *cont, double tag){

	int n;
	double twoPI_L_x,twoPI_L_y,twoPI_L_z, C_sum, S_sum, Q;
	double QMAX=15.;
	double qmax=0.;
	double XL2,YL2,ZL2;
	double rho=NMOL/(cont->x[0]*cont->x[1]*cont->x[2]);


	XL2=cont->x[0]/2;
	YL2=cont->x[1]/2;
	ZL2=cont->x[2]/2;
	double S[4000];
	int S_part[4000];
	twoPI_L_x = 2*3.1415/cont->x[0];
  	twoPI_L_y = 2*3.1415/cont->x[1];
  	twoPI_L_z = 2*3.1415/cont->x[2];

  	int order_x=50;
  	int order_y=50;
  	int order_z=50;

 	double del_q=0.1;

	for(int qq=0;qq<QMAX/del_q; qq++) {

	S[qq] = 0.0;
	S_part[qq]=0;
	}

	cout<<" "<<flush;
		for(int i=0; i<=order_x; i++){

			cout<<"*"<<flush;
			for(int j=-order_y; j<=order_y; j++) {

				for(int k=-order_z; k<=order_z; k++) {
					Q = i*i*(twoPI_L_x*twoPI_L_x)+ j*j*(twoPI_L_y*twoPI_L_y) + k*k*(twoPI_L_z*twoPI_L_z);
					Q=sqrt(Q);
					n=(int) (Q/del_q);
					if((n<QMAX/del_q)&&(n>3./del_q)){

						C_sum = S_sum = 0.0;
						for(int m=1; m<=NMOL; m++){

							double XIJ,YIJ,ZIJ,qr;
							XIJ=p[m].pos[1];YIJ=p[m].pos[2];ZIJ=p[m].pos[3];

							if(XIJ>XL2) XIJ=XIJ-2*XL2;if(XIJ<-XL2)XIJ=XIJ+2*XL2;
							if(YIJ>YL2) YIJ=YIJ-2*YL2;if(YIJ<-YL2)YIJ=YIJ+2*YL2;
							if(ZIJ>ZL2) ZIJ=ZIJ-2*ZL2;if(ZIJ<-ZL2)ZIJ=ZIJ+2*ZL2;

							qr = ( i*XIJ*twoPI_L_x + j*YIJ*twoPI_L_y + k*ZIJ*twoPI_L_z );

							C_sum+= cos(qr);S_sum+= sin(qr);
						}
						S[n]+=C_sum*C_sum + S_sum*S_sum;
						S_part[n]++;
					}	
				}
			}
		}

	ofstream dat_aus;
	string d1=".dat";
	ostringstream dateiname;
	dateiname<<"static_properties/S_(q)"<<d1;
	dat_aus.open(dateiname.str().c_str(),ios_base::out);

	if(!dat_aus){

		cout<<"Datei konnte nicht geöffnet werden RDF_Sq!";
		cout<<endl;

	}

	cout<<"\r"<<endl;
	dat_aus<<"# NMOL   "<<NMOL<<endl;
	dat_aus<<"# dichte   "<<rho<<endl;
	dat_aus<<"# q"<<"\t"<<"S(q)"<<endl;

	int z=0;
	for(int qi=1; qi<QMAX/del_q; qi++){

		if (S_part[qi]!=0){

			dat_aus<<qi*del_q<<"\t"<<S[qi]/(NMOL*S_part[qi])<<endl;
			if((S[qi]/(NMOL*S_part[qi])<S[qi+1]/(NMOL*S_part[qi+1]))&&(S[qi]/(NMOL*S_part[qi])>2.5)){

				z=1;
				qmax=(qi+1)*del_q;

			}

		}

	}

	return qmax;

}//**
