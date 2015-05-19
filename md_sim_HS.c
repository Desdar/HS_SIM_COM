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

int main(int argc, char* argv[]){

	int NMOL, num[3], num2[3], NcT;
	double rho;
	double t_intervalle = 10, tag = 0.0, int_leng = 2000;
	int s[10], interval = 0;
	double zaehler[10];
	long seed;
	double maxShape, fac, del = -1.0;
	double scale = -1.0, dichte, qmax=0.0;
	double delta_tV=0.1;
	double Press;
	double *rccluster = new double();
	Cont *cont = new Cont;
	Cell *cells;


	std::cerr << "\33[0;31m" << " *************************************************************"<< "\33[0m" << std::endl ;
	cout<<endl;
	
	int h,pinn=0;
	int r;

	h=system("mkdir -p static_properties");
	h=system("mkdir -p dyn_properties");
	h=system("mkdir -p dichteprofil");
	h=system("mkdir -p Initial");
	h=system("mkdir -p RDF");
	h=system("mkdir -p Almarza");

	for (int i; i<10;i++){zaehler[i]=0.0;}

	zaehler[0]=1.0;
	zaehler[1]=1.0;
	zaehler[2]=1.0;
	zaehler[3]=1.0;
	zaehler[4]=1.0;

//				
//*****************************************************************************
//				Parameter Input

string parameterFileName,	configFileName;

	 	if (argc < 2) 
	 	{
	    	cerr << "Too few arguments \n";
		}

    while ((argc > 1) && (argv[1][0] == '-')) {

	switch (argv[1][1]) {
	    case 'p' :
	    ++argv;
	    --argc;
		parameterFileName = argv[1];
 //       cerr << "Set control file: " << parameterFileName << endl;
		break;
	    case 'f':
	    ++argv;
	    --argc;
		configFileName    = argv[1];
  //      cerr << "Set input config file: " << configFileName << endl;
		break; 
		    default:
		cerr << "Bad option " << argv[1] << " " << argv[2] << '\n';
	//	Hsc.usage();
	}


	++argv;
	--argc;
    }

//				***************

	cout<<endl;
	cout<<"  parameterFileName 		:"<<parameterFileName<<endl;
	cout<<"  configfilename		:"<<configFileName<<endl;
	cout<<endl;

		ifstream infile;
    		
    	int check;
		char line[100];
    	char quantity[100];
    	char value[100];

		//initial  
		t_intervalle=0;


    	infile.open( parameterFileName.c_str() );
    	if ( !infile.good() ) 
		{ //open and test parameter file
        		cerr << "Parameter file:\""<<  parameterFileName << "\" not found\n";
        		exit(8);
    	}
		int nn=0;
    	
        //read parameters 
    	while (infile.peek() != EOF) 
		{
       		infile.getline(line,sizeof(line));
			nn++;
       		//allow whitespace or #-delimited comments
       	 	check = sscanf(line, "%s", quantity );
       		if( check <= 0 ) continue;
        	if( quantity[0] == '#' )continue;
	        if( quantity[0] == '\n' ) continue;
        	check = sscanf(line, "%s%s", quantity, value);
		//	
	if (check < 2) 
	{
		cerr << "Warning: wrong line format in " << parameterFileName << endl;
		cerr << "Line was: _" << line << "_\n";
	}	


if (strcmp(quantity,"NMOL")==0)			{ NMOL = atoi( value ); }
else if (strcmp(quantity,"tag")==0)		{ tag = atof(value); }
else if (strcmp(quantity,"rho")==0)		{ rho = atof(value); }

else if (strcmp(quantity,"L_x")==0)             { cont->x[0] = atof(value); }
else if (strcmp(quantity,"L_y")==0)     	{ cont->x[1] = atof(value); }
else if (strcmp(quantity,"L_z")==0)  		{ cont->x[2] = atof(value); }

else if (strcmp(quantity,"t_intervalle")==0)	{t_intervalle=atof( value );}
else if (strcmp(quantity,"seed")==0)		{ seed  = atoi( value ) ;}
else if (strcmp(quantity,"dichte")==0)		{ dichte = atof( value ) ;}

else if (strcmp(quantity,"RCCluster")==0)	{*rccluster=atof(value);}
else if (strcmp(quantity,"Pressure")==0)	{Press=atof(value);}
else if (strcmp(quantity,"Nc")==0)		{NcT=atof(value);}

else if (strcmp(quantity,"interval_length")==0)	{int_leng=atof(value);}
else if (strcmp(quantity,"MaxShapeChange")==0)	{maxShape=atof(value);}
else if (strcmp(quantity,"number_of_columns_ConfigFile")==0)  		{ r = atoi( value ) ;}

else	{cerr << "Error: Parameter line not recognised: " << line << endl;
	exit( 8 );}
}

if(seed==0) { seed=time(NULL); }
infile.close();

cont->update();

//				Parameter Input
//*****************************************************************************
//				Parameter Control

cout<<" # *************  Parameters  *****************************"<<endl;
cout<<endl;

cout<<" # number of hard spheres        "<<NMOL<<endl;
cout<<" # initial system size		"<<cont->x[0]<<"  *  "<<cont->x[1]<<"  *  "<<cont->x[2]<<endl;
cout<<" # initial number density        "<<NMOL/(cont->V)<<endl;
cout<<" # initial packing fraction      "<<3.14159/6.*NMOL/(cont->V)<<endl;
cout<<" # time interval  	  	"<<t_intervalle<<endl;
cout<<" # interval length		"<<int_leng<<endl;
cout<<" # initial time		        "<<tag<<endl;
cout<<" # desired density     	 	"<<dichte<<endl;
cout<<" # Used seed			"<<seed<<endl;
cout<<" # RCCluster			"<<*rccluster<<endl;
cout<<" # Pressure			"<<Press<<endl;




//				Parameter Control
//*****************************************************************************
//				More Definitions 

   	for(int i=0;i<3;i++)
	{
		num[i]=(int)(cont->x[i]/1.01);
		cont->nxCell[i]=cont->x[i]/num[i];
	}

Particle* p = new Particle[NMOL];
Cluster* clusters = new Cluster[NMOL];
clusters->Size = new long[NMOL];

int *Acc_TV = new int[1];
double *Rate_TV = new double[1];
int *Moves_TV = new int[1];

Acc_TV[0]=0;Acc_TV[1]=0;
Moves_TV[0]=0;Moves_TV[1]=0;
Rate_TV[0]=0.0;Rate_TV[1]=0.0;

//				****************

Box*** cell = new Box**[50];

for(int i=0;i<50;i++){
	cell[i] = new Box*[50];
	for(int j=0;j<50;j++){
		cell[i][j] = new Box[50];
	}
}

Timelist liste[NMOL];	
struct Node *n[2*NMOL];

for(int i=1; i<=2*NMOL;i++)
	{
		n[i]=(Node*) malloc(sizeof(Node));	
		n[i]->label=i;	
	}


//				More Definitions
//*****************************************************************************
//				Setup

read_in_Config( p, NMOL, cont, seed, configFileName ,r ,0,tag);

setupList(cont, rccluster, p, NMOL, cells);

for (int i=0;i<NMOL;i++){
	clusters->Size[i] = 0;
}

clearClusters(NMOL, clusters,p);

bool over = false;

	over = overlap(p,NMOL,cont);

int lower = (tag + 0.0001)/t_intervalle;
ini_pos_meas( p, NMOL,lower, s, t_intervalle);
int cluster = 10, solid;

zelle(cell, p, cont, NMOL, num);
if (qmax<0.01) qmax=S_q( p, NMOL, cont,0);

cout<<" # qmax	"<<qmax<<endl<<endl;
cout<<" # measurements will start at  : "<<s[0]<<"  "<<s[1]<<"  "<<s[2]<<"  "<<s[3]<<endl<<endl;
cout<<" # ************************************************************"<<endl;

//				Setup
//*****************************************************************************
//				Loop

for(int nn=0;nn<int_leng;nn++){

	if(nn==s[0]) pos_mes(p, NMOL, cont, s[0], 0, t_intervalle);
	if(nn==s[1]) pos_mes(p, NMOL, cont, s[1], 1, t_intervalle);
	if(nn==s[2]) pos_mes(p, NMOL, cont, s[2], 2, t_intervalle);
	if(nn==s[3]) pos_mes(p, NMOL, cont, s[3], 3, t_intervalle);	

	if(t_intervalle<0){exit(1);}

	bool overlap = false;
	double rijsq, rij[3];

	for(double t=0;t<t_intervalle;t=t+delta_tV){
		move(p, cell,delta_tV,NMOL, cont, num, liste,n,tag, s, zaehler, lower,1, qmax);

	}
}

//				Loop
//*****************************************************************************
//				x


}//**
