/*
 spic_neut:
 * Stream of neutral particles are pushed 
   from one side of the domain to the other side.
 * DSMC collision algorithm is attached.
 * Two types of particles; A and B type. A is moving B is static. 
 * Calculations are done in un-normalized quantities.
 * Three velocity components. v[2] has been given a drift. 
 ****************************************************
 Developer:
 DR. RAKESH MOULICK, CPP-IPR, ASSAM, INDIA
 ****************************************************
 */

# include <iostream>
# include <cmath>
# include <cstdlib>
# include <vector>
# include <list>
# include <ctime>
# include <random>
# include <cstring>
# include <fstream>
# include <sstream>
# include <chrono>
# include <cassert>
extern "C"
{
	# include "iniparser.h"
}
using namespace std;

/***************** INIPARSER ***************/
int  parse_ini_file(char * ini_name);

/* Random Number Generator */
std::mt19937 mt_gen(0);
std::uniform_real_distribution<double> rnd_dist(0,1.0);
double rnd()
{
    return rnd_dist(mt_gen);
}

/* Define universal constants */
double EPS;
double K;
double chargeE;
double eV;
double e;
double AMU;
double massI;
double massA;
double massB;
double massE;
double EV_TO_K;
double pi;

/* Define Simulation Parameters*/
double neut_density;  // Neutral Density
double stepSize;      // Cell Spacing
double DT;            // Time steps
double tempA, tempB;         // Temperature of the neutral species in eV
double spwt; 
double sigma_cr_max; 

/*Simulation Normalizing Parameters*/
double ndA, ndB;

/* CHANGED TYPE FROM CONST TO VAR FOR INPUT DATA CONTROL  */
int NC;              // Total number of cells
int NUM_TS;          // Total Time steps (default)
int write_interval;
int dsmc_freq; 
int nParticlesA, nParticlesB; 

string output; // Open up the output file

/* Class Domain: Hold the domain paramassEters*/
class Domain
{
public:
    int ni;      // Number of nodes
    double x0;   // initial position
    double dx;   // cell spacing
    double xl;   // domain length
    double xmax; // domain maximum position
    int num_colls = 0; // number of collisions

    /* Field Data structures */
    double *phi; // Electric Potential
    double *ef;  // Electric field
    double *rho; // Charge Density
};

/* Class Particle: Hold particle position, velocity and particle identity*/
class Particle
{
public:
    double pos;  // particle position
    double vel[3]; // particle velocity
    int id;  // hold particle identity

    // Add a constructor
    Particle(double x, double v[]):pos(x){vel[0] = v[0]; vel[1] = v[1]; vel[2] = v[2];}
};

/* Class Species: Hold species data*/
class Species
{
public:
	int nc = 0;
    // Use linked list for the particles
    list<Particle> part_list;
    double mass;
    double charge;
    double spwt;
    string name;

    int NUM;
    double Temp;
    double *den;
    double *xvel;
    double *yvel;
    double *zvel;

    void add(Particle part)
    {
        part.id=part_id++;
        part_list.push_back(part);
    }

    // Add a constructor
    Species(string name, double mass, double charge, double spwt, int NUM, double Temp)
    {
        setName(name);
        setMass(mass);
        setCharge(charge);
        setSpwt(spwt);
        setNum(NUM);
        setTemp(Temp);
    }

    // Define the constructor functions
    void setName(string name){this->name = name;}
    void setMass(double mass){this->mass = mass;}
    void setCharge(double charge){this->charge = charge;}
    void setSpwt(double spwt){this->spwt = spwt;}
    void setNum (int NUM){this->NUM = NUM;}
    void setTemp(double Temp){this->Temp = Temp;}

private:
    int part_id = 0;
};

// Define Domain and File as the global variable
Domain domain;

// Open files to hold numerical data
FILE *file_res; // File to hold the data of overall data
FILE *f1; // file to hold the A-Type particle data
FILE *f2; // file to hold B-Type particle data
FILE *file_ke; // file to hold the kinetic energy data

// Define Helper functions
void ScatterSpecies(Species *species);
void ScatterSpeciesVel(Species *species);
void PushNeutrals(Species *species);
double SampleVel(double T, double mass); 
double ComputeKE(Species *species);
void WriteKE(double Time, Species *A, Species *B);

// [Write Outputs]
void Write_ts(int ts, Species *nA, Species *nB);
void WriteParticle(FILE *file, int ts, Species *species);
void SampleNeut(Species *species, double T, double mass);
void SampleBtype(Species *species, double den0, double temp, double x1, double x2);
double XtoL(double pos);
double gather(double lc, double *field);
void scatter(double lc, double value, double *field);
double Pos(double lc);
int Xtoi(double pos);

double eval_sigma(double rel_g);
double mag(double v[]);
void CollideCEX(double *vel1, double *vel2, double mass1, double mass2);
void Collide (double *v1, double *v2, double m1, double m2);
double CollideDSMC(vector<Species> &species_list, double dt, double sigma_cr_max);

// [Ini Parser File]
int parse_ini_file(char * ini_name)
{
    dictionary  *   ini;

    ini = iniparser_load(ini_name);
    if (ini==NULL) {
        fprintf(stderr, "cannot parse file: %s\n", ini_name);
        return -1 ;
    }
    //iniparser_dump(ini, stderr); // Comment out to fix issues with iniparser

    /*[file]*/
    output		  = iniparser_getstring(ini,"file:output",NULL);

	/*[constants]*/
	EPS = iniparser_getdouble(ini,"constants:EPS",-1.0);
	K   = iniparser_getdouble(ini,"constants:K",-1.0);
	eV  = iniparser_getdouble(ini,"constants:eV",-1.0);
	e   = iniparser_getdouble(ini,"constants:e",-1.0);
	AMU = iniparser_getdouble(ini,"constants:AMU",-1.0);
	EV_TO_K  = iniparser_getdouble(ini,"constants:EV_TO_K",-1.0);
	pi  = iniparser_getdouble(ini,"constants:pi",-1.0);
      
      
    /*[time]*/    
    NUM_TS     = iniparser_getint(ini,"time:NUM_TS",-1);
    write_interval = iniparser_getint(ini,"time:write_interval",-1);
    DT = iniparser_getdouble(ini,"time:DT",-1.0);
    
    /*[grid] */
    NC         = iniparser_getint(ini,"grid:NC",-1);
   	stepSize   = iniparser_getdouble(ini,"grid:stepSize",-1.0);
    
    /*[population]*/
    tempA = iniparser_getdouble(ini,"population:tempA",-1.0);
    tempB = iniparser_getdouble(ini,"population:tempB",-1.0);
    ndA   = iniparser_getdouble(ini,"population:ndA",-1.0);
    ndB   = iniparser_getdouble(ini,"population:ndB",-1.0);
    
    nParticlesA    = iniparser_getint(ini,"population:nParticlesA",-1);
    nParticlesB    = iniparser_getint(ini,"population:nParticlesB",-1);

    massA = iniparser_getdouble(ini,"population:massA",-1.0);
    massA = massA*AMU;
    massB = iniparser_getdouble(ini,"population:massB",-1.0);
    massB = massB*AMU;
    
    /*[collision]*/
    dsmc_freq = iniparser_getint(ini,"collision:dsmc_freq",-1);
    sigma_cr_max = iniparser_getdouble(ini,"collision:sigma_cr_max",-1.0);

    iniparser_freedict(ini);
    return 0;
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/********************* MAIN FUNCTION ***************************/
int main(int argc, char *argv[])
{
    /************* INIPARSER ********************/
    if(argc<2) {
      cout<<"ERROR, at least one argument expected (the input file)."<<endl;
      exit (EXIT_FAILURE);
    }
    parse_ini_file(argv[1]);
    /********************************************/
      
    double Time = 0;
    /*Construct the domain paramassEters*/
    domain.ni = NC+1;
    DT = DT; //wp*DT; // This is the normalized time interval
    domain.dx = stepSize; //stepSize/LD;
    domain.x0 = 0;
    domain.xl = (domain.ni-1)*domain.dx;
    domain.xmax = domain.x0 + domain.xl;  

    
    cout << "DX: " << stepSize << '\t' << "DT: " << DT << endl;
    

    /*Species Info: Create vector to hold the data*/
    vector <Species> species_list;
   
    /* Add Background Neutrals */
    /************************************************************************/
    double spwtA = (ndA*domain.xl)/nParticlesA;
    double spwtB = (ndB*domain.xl)/nParticlesB;
    
    /* Create the species lists*/
    species_list.emplace_back("A", massA, 0, spwtA, nParticlesA, tempA);
    species_list.emplace_back("B", massB, 0, spwtB, nParticlesB, tempB);     

    /*Assign the species list as ions and electrons*/    
    Species &A = species_list[0];
    Species &B = species_list[1];

    /*Initiate the species density fields*/    
    A.den = new double[domain.ni];
    B.den = new double[domain.ni];
    
    //nA.xvel = new double[domain.ni];
    //nA.yvel = new double[domain.ni];
    //nA.zvel = new double[domain.ni];      

    // Sample-A type particles
    SampleNeut(&A, tempA, A.mass); 
    // Sample-B type particles
    SampleNeut(&B, tempB, B.mass); 
    
    // Scatter the species
    ScatterSpecies(&A);
    ScatterSpecies(&B);   
    

    /*------------- Print Output ---------------*/
    /*create a folder named output and
	delete the previous output folder: print statement is just to show*/
    printf("Deleting the output folder ... \n");
    system(("rm -rf "+ output).c_str());
    /*create an output folder*/
    //system("mkdir output");
    system(("mkdir "+ output).c_str());
    system(("mkdir "+ output + "/files").c_str());
    char NAME[50];
    sprintf(NAME,"%s/files/results.txt",output.c_str());
    file_res = fopen(NAME,"w");
    sprintf(NAME,"%s/files/ke.txt",output.c_str());
    file_ke = fopen(NAME,"w");
    /*------------- ----------- ---------------*/

    /*MAIN LOOP*/
    clock_t start = clock();
    for (int ts=0; ts<NUM_TS+1; ts++)
    {        
        // ScatterSpecies(&nA);
        // ScatterSpecies(&nB);

        //Perform Collision
        if(ts%dsmc_freq == 0)
        {
            sigma_cr_max = CollideDSMC(species_list, DT*dsmc_freq, sigma_cr_max);	
        }
        
        // push the particles forward in time
	PushNeutrals(&A); // scatter algorithm is called automatically after push
        PushNeutrals(&B); // scatter algorithm is called automatically after push   

        //printf("Hi, end of domain, complete upto this!\n");

        //Write diagnostics
        if(ts%write_interval == 0)
        {                      
            //===========================================================================================
            sprintf(NAME, "%s/A%d.txt", output.c_str(),ts);
            f1 = fopen(NAME,"w");

            sprintf(NAME, "%s/B%d.txt", output.c_str(),ts);
            f2 = fopen(NAME,"w");
            
            /*print diagnostics to the screen*/
	        printf("TS: %i \t nA:%ld \t nB:%ld\n", ts, A.part_list.size(), B.part_list.size());

	        /*Write time evolution of plasma profiles and Kinetic energy*/
    	    Write_ts(ts, &A, &B);
            WriteKE(Time, &A, &B);
            
            WriteParticle(f1, ts, &A);
            WriteParticle(f2, ts, &B);                        
	       	    	
        }
        
        Time += DT;
    }
    clock_t end = clock();

    /*free up memory*/
    cout << "Time = " << ((end-start)/(double)CLOCKS_PER_SEC)/60 << " minutes" << endl;
    return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/********************* HELPER FUNCTIONS ***************************/
/*Initialize the particle data : initial positions and velocities of each particle of each species*/
void SampleNeut(Species *species, double T, double mass)
{
    double v_th = sqrt(2*K*T/mass); // normalized thermal velicty of the neutrals
    
    /*sample all the particles which got created at a particular time-step*/
	for (int p=0; p<species->NUM; p++)
	{		
		double pos = domain.x0 + rnd()*(domain.ni-1)*domain.dx;
		
		/*pick a random angle*/
        double theta = 2*pi*rnd();
 
        double R = -1.0+2*rnd();    /*pick a random direction for n[2]*/
        double a = sqrt(1-R*R);
 
        double n[3];
        n[0] = cos(theta)*a;
        n[1] = sin(theta)*a;
        n[2] = R;
                
        double vel[3];
        vel[0] = v_th*n[0]; 
        vel[1] = v_th*n[1]; 
        vel[2] = v_th*n[2]; 
		
		/*reverse if going in wrong direction*/
        if (vel[2] < 0) vel[2]=-vel[2];
		
		/*add to list*/
        species->add(Particle(pos,vel));			
	}
}
//*******************************************************
void PushNeutrals(Species *species)
{
    list<Particle>::iterator it = species->part_list.begin();
    while(it!=species->part_list.end())
    {
        Particle &part = *it;         
        
        // double lc = XtoL(part.pos);
        
        part.pos += part.vel[2]*DT;

        // Apply periodic boundary condition
        if (part.pos < domain.x0)
		{
			part.pos = part.pos + domain.xl;
		}
		else if(part.pos >= (domain.x0 + domain.xl))
		{
			part.pos = part.pos - domain.xl;
		}
			it++;             

    }
    // Update number density and velocity based on the new location of partivles
    ScatterSpecies(species);
}
//*******************************************************
/*Scatter the particles to the massEsh for evaluating densities*/
void ScatterSpecies(Species *species)
{
    /*grab a pointer to the species density data and change
     the density field using the pointer*/

    double *field = species->den;

    /*clear the field*/
    memset(field,0,sizeof(double)*domain.ni);

    /*scatter particles to the massEsh*/
    for(auto &p:species->part_list)
    {
        double lc = XtoL(p.pos);
        scatter(lc,species->spwt,field);
    }

    /*divide by cell volume*/
	/*Hint: we divide by the cell volume because we have scattered
	the spwt to the grid, which is just a number. But, number per
	unit volume is density. Hence we further divide the field by the cell volume*/
    for(int i=0; i<domain.ni; i++){
    	field[i] /=(domain.dx);}
	
    field[0] *=2.0;
    field[domain.ni-1] *= 2.0;
}
//*********************************************************
/*Covert the physical coordinate to the logical coordinate*/
double XtoL(double pos)
{
    double li = (pos-domain.x0)/domain.dx;
    return li;
}
//*******************************************************
/* Returns the integer out of lc */
int Xtoi(double pos)
{
    double lc=XtoL(pos);
    return (int) lc;
}
//*******************************************************
double Pos(double lc)
{
    return domain.x0 + domain.dx*lc;
}
//*******************************************************
/*scatter the particle data to the massEsh and collect the densities at the massEsh */
void scatter(double lc, double value, double *field)
{
    int i = (int)lc;
    double di = lc-i;
    field[i] += value*(1-di);
    field[i+1] += value*(di);
}
//*******************************************************
/* Gather field values at logical coordinates*/
double gather(double lc, double *field)
{
    int i=(int)lc;
    double di = lc-i;
    double val = field[i]*(1-di) + field[i+1]*(di);
    return val;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*Write the output with Time*/
void Write_ts(int ts, Species *nA, Species *nB)
{
    for(int i=0; i<domain.ni; i++)
    {
        fprintf(file_res,"%g \t %g \t %g\n", i*domain.dx, nA->den[i], nB->den[i]);

    }
    fflush(file_res);
}
/*******************************************************/
/* Write the Output results*/
void WriteParticle(FILE *file, int ts, Species *species)
{
    for(auto& p: species->part_list)
    {
        fprintf(file,"%g \t %g \t %g \t %g\n",p.pos, p.vel[0], p.vel[1], p.vel[2]);
    }
    fflush(file_res);
}

//*******************************************************
void SampleBtype(Species *species, double den0, double temp, double x1, double x2)
{
    // Find the logical coordinates of x1 and x2
    int i1 = (int) XtoL(x1);
    int i2 = (int) XtoL(x2);
    
    // Catch for any possible error
    if(i1<0) i1 = 0;
    if(i2>=domain.ni-1) i2 = domain.ni-1;
    
    // Fractional number of particles per cell
    double np_f = den0*domain.dx/species->spwt;
    double rem = 0;
    for(int i=i1; i<=i2; i++)
    {
        // number of macro-particles to create in this cell
        int np = (int)(np_f + rem); 
        // New reminder number of particles
        rem = np_f+rem - np;
        // Space available per particle in the cell
        double spacing = domain.dx/np;

        for(int p=0; p<np; p++)
        {
            // pick a random position
            double lc = i + p*spacing;
            double x = Pos(lc);
            double v[3];

            // pick a random angle
            double theta = 2*pi*rnd();
            double v_th = sqrt(2*K*temp/species->mass);
            double R = -1.0 + 2*rnd();
            double a = sqrt(1 - R*R);

            double n[3];
            n[0] = cos(theta)*a;
            n[1] = sin(theta)*a;
            n[2] = R;

            v[0] = v_th*n[0];
            v[1] = v_th*n[1];
            v[2] = v_th*n[2];
            species->add(Particle(x,v));
        }

    }
}
//*******************************************************
double eval_sigma(double rel_g)
{
    return 1e-16;
}
//*******************************************************
double mag(double v[])
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
//*******************************************************
void CollideCEX(double *vel1, double *vel2, double mass1, double mass2)
{
    double swap;
    for(int i=0; i<3; i++)
    {swap = vel1[i]; vel1[i] = vel2[i]; vel2[i] = swap;}
}
//*******************************************************
void Collide (double *v1, double *v2, double m1, double m2)
{
	double cm[3];
	for(int i=0; i<3; i++)
	{
		cm[i] = (m1*v1[i] + m2*v2[i])/(m1 + m2);
	}
	
	/*relative velocity*/
	double cr[3];
	for (int i=0; i<3; i++)
	{
		cr[i] = v1[i] - v2[i];
	}
	
	double cr_mag = mag(cr);
	/*Pick two random angles, as per Bird's VHS model*/
	double cos_chi = 2*rnd()-1;
	double sin_chi = sqrt(1 - cos_chi*cos_chi);
	double eps = 2*pi*rnd();
	
	/*perform rotation*/
	cr[0] = cr_mag*cos_chi;
	cr[1] = cr_mag*sin_chi*cos(eps);
	cr[2] = cr_mag*sin_chi*sin(eps);
	
	/*Post collision velocities*/
	for(int i=0; i<3; i++)
	{
		v1[i] = cm[i] + m2/(m1+m2)*cr[i];
		v2[i] = cm[i] - m1/(m1+m2)*cr[i];
	}
}
/*********************************************************************************/
double CollideDSMC(vector<Species> &species_list, double dt, double sigma_cr_max)
{
	/* Sorting particles to the cell*/
	//cout << "Entered DSMC module" << endl;
	
	/* Create a vector of pointers of type Particle */	
	vector<Particle*> *part_cell;
	
	/* Allocate the size of the array */
	part_cell = new vector<Particle*> [domain.ni];
	
	/* Use this array to store the species index for the particles */
	vector<int> *part_species_index;
	part_species_index = new vector<int> [domain.ni];
	
	/* sort particles to cell */
	for(size_t s=0; s<species_list.size(); s++)
	{
		Species &species = species_list[s];
		assert(species.spwt==species_list[0].spwt);
		
		for(Particle &part:species.part_list)
		{
			int i= Xtoi(part.pos);
			part_cell[i].push_back(&part);
			part_species_index[i].push_back(s);
		}
	}
	
	double sigma_cr_max_temp=0;
	double dV = domain.dx;
	double Fn = species_list[0].spwt;
	domain.num_colls = 0;
	
	/* Perform collision for all cells */
	for(int i=0; i<domain.ni; i++) 
	{
		vector<Particle*> &part = part_cell[i];
		
		/* count the number of particles in the cell */
		int np = part.size();
		
		/* if the number of particles in a cell is less than 2 then jump */
		if(np<2) continue;
		
		/* compute number of groups according to No time counter (NTC) */
		double ng_f = 0.5*np*np*Fn*sigma_cr_max*dt/dV;
		int ng = (int)(ng_f + rnd()); /* number of groups, rounded off to an integer */
		//cout << "ng: " << ng << endl;
		/* assumes at least two particles per cell */
		for(int g=0; g<ng; g++)
		{			
			int p1, p2;
			p1 = (int)(rnd()*np);
			
			do{
				p2 = (int)(rnd()*np);
			} while(p2==p1);
			
			/* compute relative velocity */
			double cr_vec[3];
			
			for(int i=0; i<3; i++)
				cr_vec[i]=part[p1]->vel[i] - part[p2]->vel[i];
			
			double cr = mag(cr_vec);
			
			/* evaluate cross section */
			double sigma = eval_sigma(cr);
			
			/*eval sigma_cr*/
			double sigma_cr=sigma*cr;
			
			/*update sigma_cr_max*/
			if (sigma_cr>sigma_cr_max_temp)
				sigma_cr_max_temp=sigma_cr;
			
			/*eval prob*/
			double P=sigma_cr/sigma_cr_max;
			
			/*did the collision occur?*/
			if (P>rnd())
			{
				domain.num_colls++;
				double m1 = species_list[part_species_index[i][p1]].mass;
				double m2 = species_list[part_species_index[i][p2]].mass;
				//collideCEX(part[p1]->vel, part[p2]->vel, m1, m2);
				Collide(part[p1]->vel, part[p2]->vel, m1, m2);				
			}		
		}
	}
	/* Memory Clean ups */
	delete[] part_cell;
	delete[] part_species_index;
	
	if(domain.num_colls) return sigma_cr_max_temp;
	else return sigma_cr_max;
}

void WriteKE(double Time, Species *A, Species *B)
{
    double ke_A = ComputeKE(A);
    double ke_B = ComputeKE(B);

    fprintf(file_ke,"%g \t %g \t %g\n",Time, ke_A, ke_B);
    fflush(file_ke);
}

double ComputeKE(Species *species)
{
    double ke = 0;
    for (auto &p:species->part_list)
    {
        // un-normalize the velocity by multiplying with the cold thermal velocity
        ke += (p.vel[0]*p.vel[0] + p.vel[1]*p.vel[1] + p.vel[2]*p.vel[2]);
    }
    /*Multiply 0.5*mass for all particles*/
    ke *= 0.5*(species->spwt*species->mass);   

    // Calculate the average kinetic energy   
    ke = ke/(species->NUM*species->spwt);
    return ke;
}
