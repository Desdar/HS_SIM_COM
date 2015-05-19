#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <sstream>
#include <fstream>
#include <time.h>
#include "md_sim.h"
#include <iomanip>

//*****************************************************************************

void read_in_Config(Particle *p, int NMOL , Cont* cont, int lala , string name , int r, int pinn,double tag)
{
	FILE* File = fopen( name.c_str(),"r");
	if(pinn)cout<<" # some particles will be pinned"<<endl;
	for(int  i = 1; i <=NMOL;i++)
	{
		if(r==6)fscanf(File,"%lf %lf %lf %lf %lf %lf \n", &p[i].pos[1], &p[i].pos[2], &p[i].pos[3],&p[i].vel[1], &p[i].vel[2], &p[i].vel[3]);
		else if(r==3)fscanf(File,"%lf %lf %lf \n", &p[i].pos[1], &p[i].pos[2], &p[i].pos[3]);
		else if(r==14)fscanf(File,"%lf %lf %lf %lf %lf %lf %*d %*d %*d %*d %*d %*d %*d %*d \n", &p[i].pos[1], &p[i].pos[2], &p[i].pos[3],&p[i].vel[1], &p[i].vel[2], &p[i].vel[3]);

		else if(r==4)
		{

			fscanf(File,"%lf %lf %lf %d \n", &p[i].pos[1], &p[i].pos[2], &p[i].pos[3], &p[i].weight);
			if(!pinn)	p[i].weight=0;
//			cout<<" hoho "<<i<<"  "<<p[i].weight<<endl;

		}
		else 
		{
			cout<<" no read in procedure defined in function read_in( .. ) . "<<endl;
			exit(8);
		}
		p[i].pos[1]+=cont->x[0]/2;
		p[i].pos[2]+=cont->x[1]/2;
		p[i].pos[3]+=cont->x[2]/2;

		p[i].pos_ini[1]=p[i].pos[1];
		p[i].pos_ini[2]=p[i].pos[2];
		p[i].pos_ini[3]=p[i].pos[3];

		p[i].pos_tot[1]=p[i].pos[1];
		p[i].pos_tot[2]=p[i].pos[2];
		p[i].pos_tot[3]=p[i].pos[3];

		p[i].diffusion=0.;
		p[i].label=i;
	}
	fclose(File);

	ini_v(p, NMOL,lala,tag);
	bool t = overlap(p, NMOL, cont);
//	sleep(8);
	return;
}


//*****************************************************************************
