#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include<sstream>
#include<fstream>
#include <time.h>
#include "md_sim.h"
#include <iomanip>

using namespace std;

//-----------------------------------------------------------------------------

bool overlap(Particle* p, int NMOL, Cont* cont){
bool t=false;
double del[0],tdel;
for(int i=0;i<NMOL;i++){
	
	for(int j=0;j<NMOL;j++){

		if((p[i].pos[j]<0.0)||(p[i].pos[j]>cont->x[j])) cout<<"Error pos:	"<<i<<"	"<<p[i].pos[j]<<"	"<<cont->x[j]<<endl;

	}

	for(int j=i+1;j<NMOL;j++){

		del[i]=p[i].pos[0]-p[j].pos[0];
                del[i]=p[i].pos[1]-p[j].pos[1];
                del[i]=p[i].pos[2]-p[j].pos[2];
		
		for(int k=0;k<3;k++){

			if(del[k]>cont->x[k]/2)del[k]-=cont->x[k];
			if(del[k]<cont->x[k]/2)del[k]-=cont->x[k];

		}
		tdel = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

		if(tdel<-0.0001){
			cout<<" ueberlapp "<<" I=  "<<i<<"  J=  "<<j<<"   "<<tdel<<endl;
			t=true;
		}
	}

}
 return t;
}//**
//-----------------------------------------------------------------------------

void setupList(Cont *cont, double *rccluster, Particle *p, int NMOL,Cell *cells){

	int icell;
	double rc = *rccluster;

	cont->update();

	for(int i=0;i<3;i++)
	{
		cont->nxCell[i] = int(cont->x[i]/rc);
		cont->xCell[i] = cont->x[i]/cont->nxCell[i];
	}
	cont->nCells = cont->nxCell[0]*cont->nxCell[1]*cont->nxCell[2];
	
	if ((cont->nxCell[0]<2)||(cont->nxCell[1]<2)||(cont->nxCell[2]<2))
	{
	cout<<"Error: Box dimension is less than 2*rc"<<endl;
	exit(8);	
	}	

	cells= new Cell[cont->nCells];

	for(int i=0;i<NMOL;i++){

		icell = int((p[i+1].pos[1])/cont->xCell[0]) + int((p[i+1].pos[2])/cont->xCell[1])*cont->nxCell[0] 
						+ int((p[i+1].pos[3])/cont->xCell[2])*cont->nxCell[0]*cont->nxCell[1];
		if (icell<cont->nCells) {

			p[i].insertToCell(cells[icell]);

		}

		else {

			cout<< "Error: Problem with cell assignment SetUP	" << i << endl;
			cout<<"icell:	"<<icell<<"	"<<"ncell: "<<cont->nCells<<endl;
			cout<<"System Size:	"<<cont->x[0]<<"	"<<cont->x[1]<<"	"<<cont->x[2]<<endl;

			cout<<"P pos:   "<<p[i].pos[0]<<"     "<<p[i].pos[1]<<"     "<<p[i].pos[2]<<endl;
			cout<<"P pos:	"<<p[i+1].pos[1]<<"	"<<p[i+1].pos[2]<<"	"<<p[i+1].pos[2]<<endl;
			cout<<"P pos:   "<<p[i+2].pos[1]<<"     "<<p[i+2].pos[2]<<"     "<<p[i+2].pos[2]<<endl;

			exit(8);

		}
	}

	setUpNeighbours(cells,cont);

}//**

//-----------------------------------------------------------------------------

void setUpNeighbours(Cell *cells,Cont *cont) {

int x,y,z,count,cl,cln;
  for (int k=0; k< cont->nxCell[2]; k++){
    for (int j=0; j< cont->nxCell[1]; j++){
      for (int i=0; i< cont->nxCell[0]; i++){

	count = 0;
	for (int zs = -1; zs <= 1; zs++){
	  // Modulokonstr. period. Randbed.
	  z = (k+zs + 100* cont->nxCell[2])%cont->nxCell[2]; 

	  for (int ys = -1; ys <= 1; ys++){
	    y = (j+ys + 100* cont->nxCell[1])%cont->nxCell[1];

	    for (int xs = -1; xs <= 1; xs++){
	      // don't count the middle cell

	      if (!((xs==0)&&(ys==0)&&(zs==0))){

		x = (i+xs + 100*cont->nxCell[0])%cont->nxCell[0];
		cl = i +cont->nxCell[0]*j + cont->nxCell[1]*cont->nxCell[0]*k;
		cln = z*cont->nxCell[1]*cont->nxCell[0] + y*cont->nxCell[0] + x;
		cells[cl].neighbours[count] = &cells[cln];
		count ++;

	      }
	    }
	  }
	}
      }
    }
  } 


}//**

//-----------------------------------------------------------------------------


