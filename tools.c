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

using namespace std;


//-----------------------------------------------------------------------------

void clearClusters(int NMOL, Cluster *clusters,Particle *Part) {
//cout<<"Clearing rubble"<< endl;
for(int i=1; i<=NMOL;i++){

	clusters[i].Size = 0;
	clusters[i].npart = 0;
	clusters[i].firstParticle = 0;
	clusters[i].cm[0] = 0.0;
        clusters[i].cm[1] = 0.0;
        clusters[i].cm[2] = 0.0;

	Part[i].cluster = 0;
	Part[i].conneb = 0;
	Part[i].ncheck = 0;
	Part[i].nextCl = 0;
	Part[i].prevCl = 0;
	Part[i].distcm[0] = 0.0;
        Part[i].distcm[1] = 0.0;
        Part[i].distcm[2] = 0.0;
	}
//cout<<"Cleared bridgetower!"<<endl;
}//**

//-----------------------------------------------------------------------------


void ini_pos_meas(Particle *p, int NMOL, int lower, int *s, double t_intervalle){

	ostringstream dateiname;
	int m;
	double delta_t;

	for(int pp=0;pp<3;pp++){

		s[pp]=pp*pp;
		if(lower>0){

				dateiname.str("");
				dateiname<<"Initial/ini_pos_"<<pp<<".dat";

		// read postions for measurement

				FILE* File = fopen(dateiname.str().c_str(),"r");

				if(File==NULL){

					s[pp]+=lower;
					cerr<<dateiname.str()<<endl;
					cerr<<" no position file exists "<<pp<<"	"<<s[pp]<<endl;

				}

				else {

					m=fscanf(File,"%*s %d %lf \n",&s[pp], &delta_t);
					if(delta_t!=t_intervalle) {

						s[pp]=lower+2*pp;

					}

				for(int i=0;i<NMOL;i++){

					m=fscanf(File,"%lf %lf %lf \n", &p[i].pos_0[0][pp], &p[i].pos_0[1][pp], &p[i].pos_0[2][pp]);

					}
					fclose(File);
				}

		}

	}
	return;
}//**

//-----------------------------------------------------------------------------

void zelle(Box ***cell, Particle *p, Cont *cont, int NMOL, int *num){

	int i,k,l,m;

	for(k=0;k<num[0];k++) {

		for(l=0;l<num[1];l++){

			for(m=0;m<num[2];m++){

				cell[k][l][m].anzahl=0;

			}

		}

	}

	for(i=0;i<NMOL;i++){

		p[i].box[0]=(int)(p[i].pos[0]/cont->nxCell[0]);
                p[i].box[1]=(int)(p[i].pos[1]/cont->nxCell[1]);
                p[i].box[2]=(int)(p[i].pos[2]/cont->nxCell[2]);

		for(int j=0;j<3;j++){if (p[i].box[j]==num[j]) p[i].box[j]--;};

		k=p[i].box[0];
		l=p[i].box[1];
		m=p[i].box[2];

		cell[k][l][m].anzahl++;
		if(cell[k][l][m].anzahl>NMOL){cout<<"Cell overfilled:	"<<k<<"	"<<l<<"	"<<m<<endl;exit(42);}

		cell[k][l][m].members[cell[k][l][m].anzahl]=i;

	}
	return;
}//**

//-----------------------------------------------------------------------------

void pos_mes(Particle *p, int NMOL, double *Length,int s, int t, double t_interval){

	for(int j=0;j<NMOL;j++){

		p[j].pos[0][t]=p[j].pos_tot[0];
		p[j].pos[1][t]=p[j].pos_tot[1];
		p[j].pos[2][t]=p[j].pos_tot[2];

	}

	ostringstream dateiname;
	dateiname<<"Initial/initial_positions_"<<t<<".dat";
	ofstream dat_aus;
	dat_aus.open(dateiname.str().c_str(),ios::app);
	if(!dat_aus)
	{
		cout<<"Datei konnte nicht geöffnet werden in *positionen_einlese_fuer_messung()* ! ";
		cout<<endl;
	}

	dat_aus<<"#      "<<s<<"  "<<t_interval<<endl;
	for(int j=1;j<=NMOL;j++){

		dat_aus<<p[j].pos_0[0][t]<<"\t"<<p[j].pos_0[1][t]<<"\t"<<p[j].pos_0[2][t]<<endl;

	}
	return;
}//**


