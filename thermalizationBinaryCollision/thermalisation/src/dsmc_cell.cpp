/*ADVANCED PIC 2015 Lesson 5

DSMC Collision Cell example

Use the following syntax to compile with g++

g++ -std=c++11 dsmc_cell.cpp

then run ./a.exe (Windows) or ./a.out (Linux)
The code will generate two files: 
vdf0.txt: initial velocity distribution function
vdf1.txt: final velocity distribution function
*/

#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <fstream>
using namespace std;

std::mt19937 mt_gen(0);		/*seed*/
std::uniform_real_distribution<double> rnd_dist(0, 1.0);
double rng(){return rnd_dist(mt_gen);}

const double AMU = 1.66053886e-27;	/*atomic mass*/
const double K = 1.3806488e-23;	
const double V_c=0.001*0.01*0.001;	/*cell volume*/
const double Fn = 5e8;				/*macroparticle weight*/

const double delta_t=2e-5;			/*time step length*/



//samples maxwellian velocity using the method of Birdsall
double maxw(double temp, double mass)
{
	double v_th = sqrt(2*K*temp/mass);
	
	const int M = 6;
	
	double sum=0;
	for (int i=0;i<M;i++)
	{
		sum+=rng();
	}
	return (sqrt(0.5)*v_th*(sum-0.5*M)/sqrt(M/12.0));
}

double PI = acos(-1.0);

/*data structure to hold particle data*/
struct Part
{
	double v[3];	/*3 components of velocity*/
	double mass;
};

/*selects random velocity with Maxwellian thermal distribution and a drift*/
void sampleParticle(Part *part, double v_drift, double temp, double mass)
{
	/*sample from Maxwellian*/
	double v_M = maxw(temp, mass);
	
	/*pick a random angle*/
	double theta = 2*PI*rng();
	 
	/*pick a random direction for n[2]*/
	double R = -1.0+2*rng();
	double a = sqrt(1-R*R);
	 
	part->v[0] = v_M*cos(theta)*a;
	part->v[1] = v_M*sin(theta)*a;
	part->v[2] = v_M*R + v_drift;
	
	part->mass = mass;
}



struct Hist
{
	double vel;
	int bin;
};

/*returns magnitude of a 3 component vector*/
double mag(double *v)
{
	return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

/* bins velocities into nbins between min and max */
void computeVDF(Hist *hist, double min, double max, int nbins, Part *part, int np)
{

	/*set delta range*/
	double delta = (max-min)/(nbins-1);
	
	/*shift left by half delta to center bins*/
	min -= 0.5*delta;
	max -= 0.5*delta;
	
	/*set initial values*/
	for (int i=0;i<nbins;i++) 
	{
		hist[i].vel = min+(max-min)*i/(double)(nbins-1);
		hist[i].bin = 0;
	}
	
	/*compute histogram*/
	for (int i=0;i<np;i++)
	{
		double vmag = mag(part[i].v);
		int b = (int)((vmag-min)/(delta));
		if (b<0 || b>=nbins-1) continue;
		hist[b].bin++;
	}
	
}

/*outputs histogram data to a file*/
void outputVDF(const string &file_name, Hist *hist, int nbins)
{
	ofstream out(file_name);
	
	for (int i=0;i<nbins;i++)
	{
		out<<hist[i].vel<<"\t"<<hist[i].bin<<"\n";
	}
	
	out.close();
}


/*DSMC algorithm*/

/*evaluates cross-section using a simple (non-physical) inverse relationship*/
double evalSigma(double rel_g) 
{
	return 1e-18*pow(rel_g,-0.5);
}

/* collides two particles*/
void collide(Part *part1, Part *part2)
{
	double cm[3];
	for (int i=0;i<3;i++)
		cm[i] = (part1->mass*part1->v[i] + part2->mass*part2->v[i])/(part1->mass+part2->mass);
	
	/*relative velocity, magnitude remains constant through the collision*/
	double cr[3];
	for (int i=0;i<3;i++)
		cr[i] = part1->v[i]-part2->v[i];

	double cr_mag = mag(cr);
	
	/*pick two random angles, per Bird's VHS method*/
	double theta = acos(2*rng()-1);
	double phi = 2*PI*rng();
	
	/*perform rotation*/
	cr[0] = cr_mag*cos(theta);
	cr[1] = cr_mag*sin(theta)*cos(phi);
	cr[2] = cr_mag*sin(theta)*sin(phi);
	
	/*post collision velocities*/
	for (int i=0;i<3;i++)
	{
		part1->v[i] = cm[i]+part2->mass/(part1->mass+part2->mass)*cr[i];
		part2->v[i] = cm[i]-part1->mass/(part1->mass+part2->mass)*cr[i];	
	}
}


/* performs a single iteration of DSMC*/	
double performDSMC(Part *part, int np, double sigma_cr_max)
{
	
	int nc = 0;	
	double sigma_cr_max_temp = 0;
	
	/*compute number of groups according to NTC*/
	double ng_f = 0.5*np*np*Fn*sigma_cr_max*delta_t/V_c;
	int ng = (int)(ng_f+0.5);	/*number of groups, round*/
	
	/*assumes at least two particles per cell*/
	
	for (int i=0;i<ng;i++)
	{
		
		int p1, p2;
		
		p1 = (int)(rng()*np);		/*returns some number between 0 and np-1 inclusive*/
		
		do {
		p2 = (int)(rng()*np);
		} while (p2==p1);
		
		/*compute relative velocity*/
		double cr_vec[3];
		for (int i=0;i<3;i++)
			cr_vec[i]=part[p1].v[i] - part[p2].v[i];
		double cr = mag(cr_vec);
		
		/*evaluate cross section*/
		double sigma = evalSigma(cr);
			
		/*eval sigma_cr*/
		double sigma_cr=sigma*cr;
		
		/*update sigma_cr_max*/
		if (sigma_cr>sigma_cr_max_temp)
			sigma_cr_max_temp=sigma_cr;
			
		/*eval prob*/
		double P=sigma_cr/sigma_cr_max;
		
		/*did the collision occur?*/
		if (P>rng())
		{
			nc++;
			collide(&part[p1],&part[p2]);
		}
		
	}	
	if (nc) return sigma_cr_max_temp;
	else return sigma_cr_max;	
}


/*****************MAIN*********************/
int main()
{
	/*load 1000 random particles*/
	const int np_spec = 1000;
	const int num_spec = 2;
	const int np = num_spec*np_spec;	/*total number of particles*/
	Part part[num_spec*np_spec];		/*static allocation*/
	// ---------------------------------------------------
	// double v_th=300;
	// void sampleParticle(Part *part, double v_drift, double temp, double mass)
	//----------------------------------------------------
	/*load population of Ne*/
	for (int i=0;i<np_spec;i++)
	{
		sampleParticle(&part[i], 0, 1000, 20*AMU);
	}
	
	/*load population for Ar*/
	for (int i=0;i<np_spec;i++)
	{
		sampleParticle(&part[np_spec+i], 0, 500, 40*AMU);
	}
	
	// ---------------------------------------------------
	/*allocate histogram data*/
	const int nbins = 25;
	Hist hist[nbins];
	
	computeVDF(hist, 0, 1000, nbins, part, np);
	outputVDF("vdf0.txt", hist, nbins);
	
	/*perform DSMC iterations*/
	double sigma_cr_max =1e-18;		/*initial guess*/
	
	for (int it=0;it<500;it++)
	{
		sigma_cr_max = performDSMC(part, np, sigma_cr_max);
		cout<<"Iteration: "<< it << endl;
		
	}
	// ----------------------------------------------------
	computeVDF(hist, 0, 600, nbins, part, np);
	outputVDF("vdf1.txt", hist, nbins);
	
	return 0;
}