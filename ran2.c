#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <sstream>
#include <fstream>

#include "md_sim.h"
 
/* note #undef's at end of file */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static long idum;
static long idum2;
static long iy;
static long iv[NTAB];


void sran2(int seed){

   idum = (long)(seed);

   idum2= 123456;
   iy   = 0;
   memset( iv, 0, NTAB * sizeof(char));
}

float ran2()
{
    int  j;
    long k;
    float temp;

    if (idum <= 0) {
	if (-(idum) < 1) idum=1;
	else idum = -(idum);
	idum2=(idum);
	for (j=NTAB+7;j>=0;j--) {
	    k=(idum)/IQ1;
	    idum=IA1*(idum-k*IQ1)-k*IR1;
	    if (idum < 0) idum += IM1;
	    if (j < NTAB) iv[j] = idum;
	}
	iy=iv[0];
    }
    k=(idum)/IQ1;
    idum=IA1*(idum-k*IQ1)-k*IR1;
    if (idum < 0) idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}


#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX






float gauss(){

      double S=1.0;
      double AM=0.;
 

//
//      RETURNS A GAUSSIAN RANDOM NUMBER
//
      int i;
      double V=0.;
	for(i=0;i<12;i++)
	{

		V+=ran2();		
//		cout<<ran2(idum)<<endl;
	}
      V=(V-6.0)*S+AM;
      
      return V;
     }


void ini_v(Particle *p, int NMOL, int lala,double tag)
{

	double TOTX,TOTY,TOTZ,ENERGY;


	for(int i=1;i<=NMOL;i++)
	{
		
				
		for(int j=1;j<=3;j++)
		{ 	
				p[i].vel[j]=0.1*gauss();
		}
		p[i].time=tag;
		p[i].stoesse=0;
		p[i].wall=0;
		
	
	}
//	cout<<"  x  "<<x<<"   "<<NMOL<<endl;
	ENERGY=0.0;
	TOTX=0.0;
    	TOTY=0.0;
    	TOTZ=0.0;
      	
    	for(int j=1; j<=NMOL;j++)
	{
      		ENERGY=ENERGY+p[j].vel[1]*p[j].vel[1]+p[j].vel[2]*p[j].vel[2]+p[j].vel[3]*p[j].vel[3];
      		TOTX+=p[j].vel[1];
      		TOTY+=p[j].vel[2];
      		TOTZ+=p[j].vel[3];
	}
   //   std::cout<<"energy   "<<ENERGY<<std::endl;
        ENERGY=sqrt(ENERGY/(3*NMOL));
        TOTX=TOTX/(NMOL);
        TOTY=TOTY/(NMOL);
        TOTZ=TOTZ/(NMOL);
 //     cout<<TOTX<<"  "<<TOTY<<"  "<<TOTZ<<endl;
      	for (int j=1; j<=NMOL; j++)
		{
				p[j].vel[1]=(p[j].vel[1]-TOTX)/ENERGY;
		  		p[j].vel[2]=(p[j].vel[2]-TOTY)/ENERGY;
		  		p[j].vel[3]=(p[j].vel[3]-TOTZ)/ENERGY;
		}
/*	ENERGY=0.0;
	TOTX=0.0;
    TOTY=0.0;
    TOTZ=0.0;
     	
    for(int j=1; j<=NMOL;j++)
	{
			ENERGY=ENERGY+p[j].vel[1]*p[j].vel[1]+p[j].vel[2]*p[j].vel[2]+p[j].vel[3]*p[j].vel[3];
      		TOTX+=p[j].vel[1];
      		TOTY+=p[j].vel[2];
      		TOTZ+=p[j].vel[3];
     }	
//	cout<<"energy pp   "<<ENERGY/NMOL<<endl;
//	cout<<"tot_x  "<<TOTX<<endl;
//	cout<<"tot_y  "<<TOTY<<endl;
//	cout<<"tot_z  "<<TOTZ<<endl;
//	std::cout<<" # kin energy         "<<x<<"           "<<ENERGY/x<<std::endl;

	ENERGY=0.0;
	TOTX=0.0;
    TOTY=0.0;
    TOTZ=0.0;
     */ 	
}





void initial_r(Particle *p, double *Length, int NMOL){

//     F.C.C. COORDINATES
		int III,i,j,k,L;
		double D,DU,XIK,YIK,ZIK,R;


  
  
      III=exp(log(NMOL/4.0)*0.33333333)+0.01;
      

      if(4*III*III*III<NMOL) III=III+1;
      D=Length[1]/III;
      L=0;
      DU=D/2.0;
      for( i=1;i<=III;i++){
      for( j=1;j<=III; j++){
      for( k=1; k<=III; k++){
      if(L<NMOL){
      L=L+1;
      p[L].pos[1]=(i-1)*D+0.00001*gauss();
      p[L].pos[2]=(j-1)*D+0.00001*gauss();
      p[L].pos[3]=(k-1)*D+0.00001*gauss();
      L=L+1;
      p[L].pos[1]=(i-1)*D+DU+0.00001*gauss();
      p[L].pos[2]=(j-1)*D+DU+0.0001*gauss();
      p[L].pos[3]=(k-1)*D+0.00001*gauss();
      L=L+1;
      p[L].pos[1]=(i-1)*D+DU+0.00001;
      p[L].pos[2]=(j-1)*D+0.00001*gauss();
      p[L].pos[3]=(k-1)*D+DU+0.00001*gauss();
      L=L+1;
      p[L].pos[1]=(i-1)*D+0.00001*gauss();
      p[L].pos[2]=(j-1)*D+DU+0.00001*gauss();
      p[L].pos[3]=(k-1)*D+DU+0.00001*gauss();
      }
      }
      }
      }
   
     for (i=1;i<=NMOL;i++)
     {
     	for(k=1;k<=NMOL;k++)
     	{ 
      		if(i!=k)
      		{
      		XIK=p[i].pos[1]-p[k].pos[1];
      		YIK=p[i].pos[2]-p[k].pos[2];
      		ZIK=p[i].pos[3]-p[k].pos[3];
      		R=XIK*XIK+YIK*YIK+ZIK*ZIK;
      		if(R<1.)cout<<i<<"   "<<k<<"  "<<R<<endl;
      		}
     	}
     }	
//		write(*,*)'# STARTING CO-ORDINATES'
 
}




    
