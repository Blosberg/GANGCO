/*----- IMPLEMENTATION OF THE GILLESPIE ALGORITHM ON THE MOTION OF NUCLEOSOMES, TREATED AS 1-D PARTICLES ALONG DNA
// ---last updated on  Fri Oct 17 15:49:24 CEST 2014  by  ga79moz  at location  TUM , xenopus

//  changes from  Fri Oct 17 15:49:24 CEST 2014 : made this small_particles branch the master branch in Git. Git is now used as the change-tracker, check the history of the git repository from this point forward.

//  changes from  Tue Jul 1 16:59:36 CEST 2014 : set output folder to just local/muN-..., and added an output file for the void _density_ in addition to distribution

//  changes from  Tue Feb 18 18:27:29 CET 2014 : implemented simplified run for small particles. Tested, works.

//  changes from  Thu Jan 9 12:42:40 CET 2014 : additional collection series for the 2-point correlation function: gathering curves during transient filling process now in addition to during equilibrium.

//  STARTED FROM SCRATCH BUILDING CODE FOR SMALL PARTICLES -SOME FEATURES WERE KEPT, SOME WERE AXED.
*************************************************************************************/


#include <fstream>   
#include <math.h>
#include <iostream>  
#include <iomanip>  
#include <sstream>

#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_rng.h>
#include <cmath>
#include <queue>  //in order to stack the events of binding/unbinding
#include <dirent.h>

#include "GA_GC_NTF.h"
#include <GA_NTF_standards.h>
#include <GA_absolute_standards.h>
#include <bren_lib.h>

double const pi   = 3.14159265358979323846264;
//double const hbar = 1.0546E-34; // in units of Js.
double const kB   = 1.38065E-23;  // in units of J/k


using namespace std;
const int Np10 = 10; 		// Resolution of time-samples taken. this is probably high enough resolution on the log10 scale. 
				// if Np10 = 10, then this implies that we take ten data points for each base-ten time increment

const int nsample = 30; 	//----the minimum number of binding events we want to sample over for valid statistics.
const int charlength =400;	//----the length of the string kept in memory for file-path manipulation.

int VNN_SNG_calc_smallp(double * potential, const int NNRANGE, const int a, const double E0);
//---the SNG func has an extra input parameter for the LJ potential.
int VNN_LNG_calc_smallp(double * potential, const int a, const double E0); 
//---generate potential for the "L"inear potential for small particles.

int space_tpoints_logarithmic(const double t0, const double tf, const int Np10, const int nbins, double * tx);

bool should_observe_equilibrium_voiddist( const GAdata * P, const double tau, int ** void_histogram_equilibrium  ); 
//--- tests whether it is the right time to increment the void equilibrium distribution.

int calc_void_statistics( double * void_means, double * void_stddevs, double * rhocheck, int const* const * const void_histogram, const int total_obs_filling, const int Llim , const double CGF,  const int numtrials);

int pick_index_from_rand_norm_dist( double rand01, double * distr, int L );

int get_first_N_pvals_from_config_list( const vector < configuration >   Z_t, float * target, const int N);
// ---order the configurations along 'p', and then put the first N values of that 
// ---ordered list into the target array
//******************************************************************

int main(int argc, char *argv[])
{
int i,j,x,n, t_index;    
unsigned long int seed;	// --- random number generator seed.
int Llim;	// --- the length of the strand.

bool HNG;       // --- are we dealing with HNG particles? if false, then SNG
bool SNG;       // --- soft core particles?
bool LNG;       // --- linear-potential particles?

long int M;	// the number of transitions available at the moment;

double t=0, dt_inc=0;
double t0, tf;	// tf=final, or terminating time. t0=initial time at which we record filling.

const gsl_rng_type * T;
gsl_rng * r;
double r1, r2;	//  two random number places.

double E0=0.0;
int  footprint;		// footprint of the nuc, regardless of hard or soft. m=length of TF

int h[3]; h[0]=h[1]=h[2]=0;
double muN_input;	//---the original value (without coarse-graining)
double muN;		//--- this one MIGHT get changed.
double muTF0;

int F[2];

int numtrials;

double kS_N,  kA_N;   // BASIC RATE CONSTANTS FOR SLIDING AND ADDITION/REMOVAL for nucleosomes.
double kS_TF, kA_TF;  //     "		"		"		"      for transcription factors.	

bool krm_b;	//---  do we or do we not have active remodellers in the system.
double krm_val;	//---  if we do, how strongly do they act?

bool should_plot_snapshots;
bool should_plot_kymo;

string path, pathout, output_folder;
int Nplots2makeshort, Nplots2makelong;
int Nplots2make_kymo;

double dt_obs;
double dtau_plot;

//=====================GET COMMAND LINE PARAMETERS=============
int TASKID;
string NGtype;

//-------------  'command line parameters' ----------------------------------

if(argc <= 1)
	{
	goto command_line_garbage;
	}
NGtype = argv[1];

if( (NGtype!="SNG" && NGtype!="LNG") && (NGtype!="HNG") )
	{
	goto command_line_garbage;
	}

if (  ( (NGtype=="SNG"  || NGtype=="LNG" ) && argc != 5) || (NGtype=="HNG" && argc !=4)  )
	{
	
	command_line_garbage:
	cout << "\n error: unexpected number of input arguments; argc = " << argc << " \n and they were : ";
	for (j=1;j<argc;j++)
		{
		cout << argv[j] << ",  ";
		}
	cout << "\n exiting\n\n";
	exit(1);
	}

TASKID      = atoi(argv[2]);
muN_input   = atof(argv[3]); //----this should be input in "per footprint" convention; scaled to lattice site within the code.
 

if( (NGtype=="SNG"  || NGtype=="LNG" ) )
	{
	E0            = atof(argv[4]); //---atof is for floats and doubles.
	}

//---------------------------------------------------------------------------
if(NGtype=="SNG")
	{
	SNG=true; HNG=false; LNG=false;
	}
else if(NGtype=="HNG")
	{
	HNG=true; SNG=false; LNG=false;
	}
else if(NGtype=="LNG")
	{
	LNG=true; HNG=false; SNG=false;
	}
else
	{
	cout << "\n ERROR: ambiguous NG type. exiting. \n"; 
	exit(1);
	}

double max_tcorr;

int paritycheck=0;

string BZcond; //-----flag to determine which reactions are weighted by boltzmann factors.
double BZalpha; //-- denotes the degree to which boltzmann scaling is applied on addition or on removal.

		//-- r_on = exp(-V * BZalpha) , r_off = exp( V * (1-BZalpha) )
		// so alpha ==0 implies rewarding removal, while alpha==1 implies punishing addition.


double t_trans; //---- the 'transient' time period after which equilibrium properties are analyzed.
		//---- WHETHER TRANSIENT EFFECTS ARE REALLY GONE BY THIS TIME)

//-----------------------------------------------------------------------------------------------

//  ----  READ IN PARAMETERS FOR THE SIMULATION FROM THE IN FILE.-----------------

char cpath[charlength];
clear_charray(cpath, charlength );
sprintf(cpath, "./GA_GC_NTF.in"); //----input file should always be in the working directory.
ifstream datin;

datin.open(cpath);

i=0;
while( datin.fail() ) 
	{

	cout << "\n failed accessing input file " << i << "times.";	

	i++;
//	sleep(1); // wait for a second to avoid race condition.
	if( i >= 100)
		{
		break; // --- avoid an infinite loop.
		}	
	datin.open(cpath);
	
	}
if( datin.fail() )
	{
	cout << "\n ERROR, can't find input file. exiting \n";
	exit(1);
	}

datin  >> kS_N  >> kA_N >> kS_TF  >> kA_TF;
datin  >> Llim  ;  			//---Llim in LATTICE SITES, (after CG-ing), not bp.
datin  >> t0   >> tf >>  t_trans >> dt_obs >> dtau_plot >> max_tcorr; 
//--Llim - system size; dt_obs is now just the time we _start_ looking.
datin  >> footprint ;	//----in the HNG case, we just take 'w' to mean 'k'
datin  >> muTF0;

datin  >> krm_b                 >> krm_val;
datin  >> should_plot_snapshots >> Nplots2makeshort  >> Nplots2makelong;
datin  >> should_plot_kymo      >> Nplots2make_kymo;
datin  >> BZcond                >> BZalpha;

datin  >> output_folder;
datin  >> paritycheck;  //---this number is always 8888888 in the input file. 
// If it gets read as something different then somethings wrong with the I/O formatting.

datin  >> numtrials;
datin.close();

//-----------------------------------------------------------------------------------------------


//--------- @@@ THESE OPTIONS denote per-footprint input... JUST FOR K-MER TESTING------
muN = muN_input-gsl_sf_log(double(footprint)); //---- here convert per-footprint to per-lattice site muN
E0  = E0/double(footprint);
// ---------------------@@@ DELETE down to here -----------------------------


if( paritycheck != 88885888 )
	{
	cout << "\n ERROR in parity check -input parameters are disordered somehow. Check your input file.\n";
	exit(1);
	}
//----------------------------------------------------------

bool boltzmann_on_uphill          = false;
bool boltzmann_on_add             = false;
bool boltzmann_on_removal         = false; 
bool boltzmann_on_addrem_intermed = false; //-----this condition is now read in from file.


if (BZcond =="boltzmann_on_uphill")
	{
	boltzmann_on_uphill  = true;
	}
else if(BZcond == "boltzmann_on_add")
	{
	boltzmann_on_add     = true;
	}
else if(BZcond == "boltzmann_on_removal")
	{
	boltzmann_on_removal = true;
	}
else if(BZcond == "boltzmann_on_addrem_intermed")
	{
	boltzmann_on_addrem_intermed = true;
	}
else
	{
	cout << "\n ERROR: bolzmann criteria is either conflicting or undefined.\n";
	exit(1);
	}

if( (BZalpha<0.0) || (BZalpha>1.0))
	{
	cout << "\n ERROR: BZalpha is not between zero and one. exiting.\n";
	exit(1);
	}

bool output_patterns = false;	//---should we print out the 2-part correlation function?

//----- TAKE INPUT irho_target and use it to determine chemical potential ---------
//-----altered startup condition begins here.--------------------------------------------------------------------------


double tl,tu;	//---upper and lower bounds of the time necessary for convergence
		//---using punish/reward schemes respectively.
double ltf;	//---log of the final time. This is just a dummy temp variable.



/***************************************   THIS IS JUST FOR BOLTZMANN ON RELIEF CONDITIONS ********************************


if(boltzmann_on_removal)
	{
	if (TASKID == 2)  
		{tf=0.0008;}
	else if (TASKID == 4)  
		{tf=0.08;}
	else if( TASKID == 7 )
		{tf=0.8;}
	}
else if(boltzmann_on_addrem_intermed)
	{

	tu=tf;	//----upperbound fixed as per addition rate.

	if (TASKID == 2)  
		{tl=0.0008;}
	else if (TASKID == 4)  
		{tl=0.08;}
	else if( TASKID == 7 )
		{tl=0.8;}

	ltf = gsl_sf_log(tl) + BZalpha* ( gsl_sf_log(tu) - gsl_sf_log(tl) );
	tf  = gsl_sf_exp(ltf);
	cout << "\n With intermediate Boltzmann condition alpha = " << BZalpha << ", tf is set to " << tf << endl;

	} 
	//----else if boltzmann_on_add, then just leave tf as is -don't change anything.

**************************************************************************************************************************/

int const RMRANGE  = 2*footprint;

//-------------------------- setup the output directory and file --------------------------

clear_charray(cpath, charlength );

//   sprintf(cpath, "./   muN-%.2f_E0-%.2f/", muN, E0);
//   sprintf(cpath, "./%s_muN-%.2f_E0-%.2f_perfp_k-%d/", output_folder.c_str(), muN+gsl_sf_log(footprint), E0*footprint, footprint);
sprintf(cpath, "./%s_rhoeq-1.04_E0-%.2f_perfp_k-%d/", output_folder.c_str(), E0*footprint, footprint);


pathout = cpath;

if ( DirectoryExists( pathout.c_str()  ) )
	{
	cout << "\n output path already exists. \n";
	}
else
	{
	cout << "\n output path does not exist -creating it. \n";
	clear_charray(cpath, charlength );

	sprintf(cpath, "mkdir \"./%s\"",pathout.c_str());
	system(cpath);
	}

clear_charray(cpath, charlength );

sprintf(cpath, "%sGA_cjob.log",pathout.c_str());

ofstream *log = new ofstream(cpath, std::ofstream::out);

//------------------------------------------------------------------------------



ofstream * timestamps; 
if(should_plot_snapshots )
	{
	clear_charray(cpath, charlength );
	sprintf(cpath, "%ssnapshot_timestamps.txt",pathout.c_str());
	timestamps = new ofstream;
	(*timestamps).open(cpath);
	}

//-----------------  NOW   SEED   THE   RNG ------------------------------------
clear_charray(cpath, charlength );
sprintf(cpath, "./rngSEED.in");
ifstream fseedin(cpath);

if(fseedin.fail())
	{
	*log << "\n cannot find rngseed file... exiting.\n";
	cout << "\n cannot find rngseed file... exiting.\n";
	(*log).close();
	exit(1);
	}
else
	{
	fseedin >> seed;
	fseedin.close();
	}

cout << "TASKID        = " << TASKID  << endl;
cout << "original seed = " << seed    << endl;
seed=(seed/TASKID);
cout << "new seed = " << seed  << endl;

gsl_rng_env_setup();
      
T = gsl_rng_default;
r = gsl_rng_alloc (T);

gsl_rng_set(r,seed);

r1=gsl_rng_uniform(r); 
r2=gsl_rng_uniform(r);

//------------------------- ASSIGN VALUES FOR THE REACTION ARRAY -----------------------

t=0;
int R=0; // ----the index of the reaction of interest.
int Rtype;

M = 8*Llim;	// the number of possible reactions -- most of which will have 0 amplitude
cout << endl;

//-------------------get arrays for t-filling--------------

//-------------------- linear -----------------------------

int total_obs_eq;
double * tpoints_eq    = NULL;
double ti=dt_obs, tip1 = 0.0;

if( tf > t_trans + dt_obs) //-----check if our "post-transient" period is at least long enough for one obervation.
	{

	total_obs_eq = floor((tf -t_trans)/dt_obs); //--total number of observations for equilibrium observations
	tpoints_eq= new double[total_obs_eq];
	for(i=0;i<total_obs_eq;i++)
		{
		tpoints_eq[i]=0.0;
		}
	space_tpoints_linear(t_trans, dt_obs, tf, total_obs_eq, tpoints_eq );
	}
else
	{
	total_obs_eq=0;
	}
	
	
	
int    * void_histogram_equilibrium = NULL;

if (get_voiddist_equilibrium)
	{
	void_histogram_equilibrium  = new int[Llim+1];
	for(j=0;j<=Llim;j++)
		{
		void_histogram_equilibrium[j]=0;
		}
	}
//-------------------- logarithmic -------------------------

int total_obs_filling = ceil(   Np10* log10( tf/t0)  );	//----'filling' in this context just means transient.

double tpoints_filling[total_obs_filling];

space_tpoints_logarithmic(t0, tf, Np10, total_obs_filling, tpoints_filling );

//---------------------------------------------------------------------

int    ** void_histogram  = new int*[total_obs_filling];
//-----------------   allocate and initialize overall population arrays ----------------

double OavN_eq[Llim];
double OavTF_eq[Llim];
double output2partcorr_eq[Llim];


init_array( OavN_eq           , Llim, 0.0 );
init_array( OavTF_eq          , Llim, 0.0 );
init_array( output2partcorr_eq, Llim, 0.0 );

double * output2partcorr_ti[total_obs_filling]; //----now the transient quantity.
double * OavN_ti[total_obs_filling];
double * OavTF_ti[total_obs_filling];


for(t_index=0; t_index<total_obs_filling; t_index++)
	{
	OavN_ti[t_index]            = new double[Llim];
	init_array( OavN_ti[t_index], Llim, 0.0);

        OavTF_ti[t_index]           = new double[Llim];
	init_array( OavTF_ti[t_index], Llim, 0.0);

	output2partcorr_ti[t_index] = new double[Llim];
	init_array( output2partcorr_ti[t_index], Llim, 0.0);

	}

//-------------------------------------------------------------------------------------

config_set_t Z_all_t[total_obs_filling]; //--- an array of size for the number of time-observation points.
					 //--- this is used for the entropy calculation.



//-------------------------------------------------------------------------------------


double *  void_means      = new double[total_obs_filling];
double *  void_stddevs    = new double[total_obs_filling];
double *  Ncheck          = new double[total_obs_filling];
double *  fillingfrac 	  = new double[total_obs_filling]; 	
							//--- this is the _TOTAL_ filling frac array. 
							//--- The one stored in the simdat struct is 
							//--- just for a single run, and gets 
							//--- overwritten each time.

// 2-D array for the histogram of void sizes at the various time points. 
// ---Now allocate memory:
for(i=0;i<total_obs_filling;i++)
	{
	fillingfrac[i]   = 0.0;
	void_means[i]    = 0.0;
	void_stddevs[i]  = 0.0;
	Ncheck[i]        = 0.0;

	//------------------------------------
	void_histogram[i] = new int[Llim+1];
	for(j=0;j<=Llim;j++)
		{
		void_histogram[i][j]=0;
		}

	//------------------------------------

	Z_all_t[i].tval          = tpoints_filling[i];
	Z_all_t[i].tpoint_passed = false;
	Z_all_t[i].Nave = 0.0;

	Z_all_t[i].S    = 0.0;
	Z_all_t[i].H    = 0.0;
	Z_all_t[i].Htot = 0.0;

	//----always declare this so that there's only one version of "should_observe"
	//----if calculate_entropy is false, then we just won't do anything with this.

	}

//-----------get arrays for two-TF correlation------------------------------


int nbins = floor((max_tcorr)/dtau_plot);

//--------- GET TWO BODY POTENTIALS: -----------------------------------------


int NNRANGE;
if(SNG)
	{
	NNRANGE = Llim-1;	// Lennard-Jones potential across the entire array.
	}
else if(LNG)
	{
	NNRANGE = footprint;
	}
else if(HNG)
	{
	NNRANGE = footprint;
	}

double VNN[NNRANGE+1];
double xcoarse[NNRANGE+1];

for(i=0;i<=NNRANGE;i++)
	{
	VNN[i]     = 0.0;
	xcoarse[i] = i;
	}
if(SNG)
	{
	VNN_SNG_calc_smallp(VNN, NNRANGE , footprint, VLJ_rm, E0);
	}
else if(LNG)
	{
	VNN_LNG_calc_smallp( VNN, footprint, E0 );
	}
else if(HNG)
	{
	for(i=0;i<=NNRANGE;i++)
		{
		VNN[i]     = 1.0/0.0;
		xcoarse[i] = i;
		}
	}
else
	{
	cout << "\n ERROR: no NGtype selected.\n";
	exit(1);
	}
 //--coarse-grain the 2-body interaction potential into a smaller system.
//--------------------------------SETUP INITIALIZATION PARAMETERS---------------------
int     sizes_n_ranges[4];
sizes_n_ranges[0] = footprint;
	//-------sizes_n_ranges[1] = kHNG; // kHNG *is* 'footprint'
sizes_n_ranges[1] = RMRANGE;
sizes_n_ranges[2] = Llim;
sizes_n_ranges[3] = NNRANGE;
//--------------------------------------
double  k_rates[5];
k_rates[0]=kS_N;
k_rates[1]=kA_N;
k_rates[2]=kS_TF;
k_rates[3]=kA_TF;
k_rates[4]=krm_val;
//--------------------------------------
double times[4];
times[0] = tf;
times[1] = dt_obs;
times[2] = dtau_plot;
times[3] = t_trans;
//--------------------------------------
bool flags[7];
flags[0] = HNG;
flags[1] = SNG;
flags[2] = LNG;
flags[3] = boltzmann_on_uphill ;
flags[4] = boltzmann_on_add    ;
flags[5] = boltzmann_on_removal; 
flags[6] = boltzmann_on_addrem_intermed; //-----this condition is now read in from file.
//--------------------------------------
int observations[6];
observations[0] = Nplots2makeshort;
observations[1] = Nplots2makelong;
observations[2] = Nplots2make_kymo;
observations[3] = total_obs_filling;
observations[4] = total_obs_eq;
observations[5] = nbins;
//--------------------------------------
double energetics[4];
energetics[0] = muN;	//--- input the actual chemical potential 
			//--- the constructor handles all aspects of coarse-graining.
energetics[1] = muTF0;
energetics[2] = E0;
energetics[3] = BZalpha;//----the boltzmann alpha value that defines the ratio between addition and removal.

//--------------------------------------

double num_binds_per_simulation_F0 =0.0;
double num_binds_per_simulation_F1 =0.0;

double avg_F0_occupation = 0.0;
double avg_F1_occupation = 0.0;

// these will be taken each round from the simdat data variables of the same name

//---------------------------------------------------------------------------------------
string BZ;
if( ( int(boltzmann_on_add)+int(boltzmann_on_removal) + int(boltzmann_on_uphill) + int(boltzmann_on_addrem_intermed) ) != 1)  
	{
	cout << "\n ERROR: the boltzmann conditions are not unique and/or sufficient.\n";
	*log << "\n ERROR: the boltzmann conditions are not unique and/or sufficient.\n"; 	  (*log).close();
	exit(1);
	}

if( boltzmann_on_uphill )
	{
	BZ="uphill";
	}
else if ( boltzmann_on_add )
	{
	BZ="add";
	}
else if ( boltzmann_on_removal )
	{
	BZ="remove";
	}
else if ( boltzmann_on_addrem_intermed )
	{
	BZ="addrem_intermed";
	}
else 
	{
	cout << "\n not sure what the BZ condition is. Exiting. \n";
	*log << "\n not sure what the BZ condition is. Exiting. \n"; (*log).close();
	exit(1);
	}
//------------------------ SEND DOCUMENTATION TO THE LOG FILE ---------------------------

cout  << "\n bind_irrev    = "    << bind_irrev << endl ;

*log  << "\n NGtype     = "       << NGtype;

*log  << "\n bind_irrev  = "       << bind_irrev; 
*log  << "\n fixed_ref   = "       << fixed_ref;
*log  << "\n TFs_allowed = "       << TFs_allowed;
*log  << "\n debugging   = "       << debugging;
*log  << "\n BZ    = "             << BZcond;
*log  << "\n BZalpha = "           << BZalpha;

*log  << "\n RUN beginning with the following parameters : \n";
*log  << "\n t_trans = " << t_trans << ", tf="            << tf      << endl;
*log  << "Llim="         << Llim    << endl;
*log  << "E0 = "         << E0      << ", footprint = "   << footprint      << endl ;
*log  << "h0 = "         << h[0]    << ", h1= "   << h[1]   << ", h2= "   << h[2] << endl ;
*log  << "muN0 = "       << muN     << endl ;
*log  << "uTF0 = "       << muTF0   << endl ;
*log  << "krm_b = "      << krm_b   << ", krm_val = "  << krm_val  << endl ;
*log  << " with kA_N = " << kA_N    << endl ;
*log  << " with kA_N = " << kA_N    << endl ;
*log  << " with kS_N = " << kS_N    << endl ;
*log  << " total_obs_filling = "    << total_obs_filling << endl;
*log  << " VLJ_rm =    " << VLJ_rm  << endl;


//-----------  SETUP KYMO FILES   ----------------------------

ofstream * fkymo_outN;
ofstream * fkymo_outTF; 


if(should_plot_kymo)
	{

	fkymo_outN  = new ofstream;
	fkymo_outTF = new ofstream;

	clear_charray(cpath, charlength );
	sprintf(cpath, "%skymograph_t_Npos.txt",pathout.c_str());
	(*fkymo_outN).open(cpath);

	clear_charray(cpath, charlength );
	sprintf(cpath, "%skymograph_t_TFpos.txt",pathout.c_str());
	(*fkymo_outTF).open(cpath);

	}

//=================================================================================================
//-------------------------    BEGIN numtrials SIMULATIONS HERE:    ----------------------------

for(i=0;i<numtrials;i++)
	{
	//-------setup the data structure------------------

	GAdata * simdat = new GAdata(energetics, observations, k_rates, times, sizes_n_ranges, h,F, krm_b,  path, log, timestamps, flags);

	//------- ask if we should start with some fixed initial configuration. ---------
	if( set_fixed_initial )
		{
		*log << "\n set_fixed_initial is TRUE : setting initial ones.\n";
		(*simdat).set_fixed_initial_particles( path, r );
		}
	//------------------------------------------------------------------------------

	(*simdat).tpoints_filling   = tpoints_filling; 
	(*simdat).tpoints_eq        = tpoints_eq; 

	//------------------------------------

	if(fixed_ref)	//-----start by adding the +1 at zero *IF* fixed ref
	  {
	  x      = h[0]+h[1];
	  Rtype  = 1;
	  react(x,Rtype,*simdat);
	  }


	//-----start at t=0-------
	while( (*simdat).t<tf)
		{

		if( fixed_ref && (*simdat).pos[h[0]+h[1]].state != 1)
		{
		cout << "\nERROR -lost our +1 reference";
		*log << "\nERROR -lost our +1 reference";
		exit(1);
		}

		if( !TFs_allowed && (*simdat).TFnum > 0)
		{
		cout << "\nERROR -we got some TFs in here";
		*log << "\nERROR -we got some TFs in here";
		exit(1);
		}

		dt_inc=(1.0/(*simdat).a0)*gsl_sf_log(1.0/gsl_rng_uniform(r)); //---this is how much we step forward in time this time step.
	

		(*simdat).should_observe_filling(dt_inc, void_histogram, Z_all_t ) ; 
		// updates the filling and void-histograms.

		if (get_voiddist_equilibrium)
			{//----DETERMINES WHETHER WE SHOULD DETERMINE THE EQUILIBRIUM VOID DISTRIBUTION.
			simdat->obs_count_eq_vdist += (*simdat).should_observe_equilibrium_voiddist( dt_inc, void_histogram_equilibrium  );
			}

		if ( output_patterns) 
			{
			(*simdat).should_observe_patterns(dt_inc);
			}

		(*simdat).t+=dt_inc;

		if( (*simdat).t>tf) 
			break;

	//----necessary for deselection of certain reactions----
	reroll:
	//--------------------

		r2=gsl_rng_uniform(r);
		R=(*simdat).choose_reaction(r2);	// returns the index of the Rx that should occur
					// --- both position and type.

		if(R < M)
			{
			x=floor(float(R)/float(8));
			Rtype = R-8*x;

			if((x==h[0]+h[1]) && ( (Rtype==0) || (Rtype==2) || (Rtype==3) ) && fixed_ref )
				{
		//		cout << "\n blocked a passive N0 removal";
				goto reroll;
				}

			if(  (Rtype ==5 )   && !TFs_allowed  )
				{
		//		cout << "\n blocked a TF addition";
		//		goto reroll;
				*log << "\n ERROR: TF addition attempted -but TFs are not allowed -this is an error!\n";
				(*log).close();
				exit(1);
				}

			react(x,Rtype,*simdat);
			}
		else
			{
			(*simdat).remodel(R,r);
			}


		/*  ---THIS IS TO ADD LINES TO THE FILE FOR THE KYMOGRAPH----- */

		if( (*simdat).printtime_kymo() && should_plot_kymo && i==0) // WE ONLY DO IT FOR THE FIRST RUN.
		   {
		    (*simdat).plot_snapshot_kymo( fkymo_outN , fkymo_outTF );
		   }



		//----------REDUNDANT ERROR CHECKING: -WE'RE CONFIDENT BY NOW -------------
		// ---this should only be done if/when there is a discrepency discovered.
		//	 (*simdat).check_rates(); 	
		//	 (*simdat).check_states();
		//-------------------------------------------------------------------------

		(*simdat).counter++;

	}//--end the "while(t<tf)" -end of this simulation.	

	/*---------------DONE THE INDIVIDUAL RUN--------------------------------*/
  			 
	//-----------UPDATE OVERALL AVERAGES ----------------------------
	if ( output_patterns) //---should we bother calculating the 2pc every time?
		{

		(*simdat).normalize_2_part_corr( );

		for(x=0;x<Llim;x++)	//----- EQUILIBRIUM VALUES------
			{
			output2partcorr_eq[x] += (1.0/float(numtrials))*(*simdat).two_part_corr_eq[x];

			OavN_eq[x]  += ((1.0/float(numtrials))*( (1/float(total_obs_eq))* float((*simdat).onepoint_histocc_N_eq[x])) );
			OavTF_eq[x] += ((1.0/float(numtrials))*( (1/float(total_obs_eq))* float((*simdat).onepoint_histocc_TF_eq[x]))); 
			}

		//----------------
		for(t_index=0; t_index<total_obs_filling; t_index++)
			{
			for(x=0;x<Llim;x++)	//----- TRANSIENT VALUES------
				{
				output2partcorr_ti[t_index][x] += (1.0/float(numtrials))*(*simdat).two_part_corr_ti[t_index][x];

				OavN_ti[t_index][x]  += (1.0/float(numtrials))*float( (*simdat).onepoint_histocc_N_ti[t_index][x] );
				OavTF_ti[t_index][x] += (1.0/float(numtrials))*float( (*simdat).onepoint_histocc_TF_ti[t_index][x] );
				}
			}

		}
      //-------------------------------------------------------
      (*simdat).reset();    
      (*simdat).N_remod=0; //this has to be done manually, since reset is designed for
			// other implementations and therefore cannot refer to N_remod.
			// It only exists here, and must be dealt with uniquely here.

	//-------UPDATE OVERALL FILLING FRACS------------------
    for(j=0;j<total_obs_filling;j++)
      {
	fillingfrac[j] += (1.0/float(numtrials))*((*simdat).filling_frac[j])  ;
	//fillingfrac is what we keep in scope in main. filling_frac is a member of class GAdata
      }


    if( floor(numtrials/(i+1)) <= 100 && i%(int(0.1*numtrials)) == 0)
	{
	*log << "\n at the end of run " << i << ", there were " << (*simdat).Nucnum << " particles in the system\n";

	cout << "\n completed trial number " << i << " of " << numtrials << "; " << (100*(i)/numtrials) << "% done." ;	
	*log << "\n completed trial number " << i << " of " << numtrials << "; " << (100*(i)/numtrials) << "% done." ;	
	}

    delete simdat;

    if(calculate_entropy)
	{	
	for(j=0;j<total_obs_filling;j++)
		{ //----reset configuration evaluation points to 'false' 
		  //----(i.e. haven't looked at it yet.)
		Z_all_t[j].tpoint_passed = false;
		}
	}		
    }//---this closes the for(i=0..numtrials) loops


//============================================================================================


//-------------------------DONE ALL THE TRIALS-----------------------

//! -----this was for nucleosomes---- : calc_void_statistics(  void_means, void_stddevs, Ncheck, void_histogram, total_obs_filling, Llim , CGF, numtrials);
//! removed coarse-graining functionality 
calc_void_statistics(  void_means, void_stddevs, Ncheck, void_histogram, total_obs_filling, Llim , 1.0, numtrials);

	num_binds_per_simulation_F0 = 	num_binds_per_simulation_F0 /(double (numtrials) );
	num_binds_per_simulation_F1 = 	num_binds_per_simulation_F1 /(double (numtrials) );

	avg_F0_occupation = avg_F0_occupation / ( double (numtrials) );
	avg_F1_occupation = avg_F1_occupation / ( double (numtrials) );

//--------------------------------------------
double avg_as=0.0;
int num_samp=floor(Llim/3);
int start = num_samp;

for(i=0;i<num_samp;i++)
	{
	x=i+start;
	avg_as+= (OavN_eq[x]/(double(num_samp)));	//----collect the central third population density average.
	}
//--------------------------------------------
//======================  OUTPUT THE ENTROPY (as necessary)  ===============
ofstream *fentropyout;

int ii=0;
int size_of_Zt=0;
double p=0.0;
double pnorm=0.0;


//-------------allocate and initialize the number of configuration counts--------
int get_first_N=100;
float * num_config_counts[total_obs_filling];

for(i=0;i<total_obs_filling;i++)
	{
	num_config_counts[i] = new float[get_first_N];;
	}

for(i=0; i< total_obs_filling ; i++)
	{
	for(j=0;j<get_first_N;j++)
		{
		num_config_counts[i][j]=0.0;
		}
	}

//-----------------------------------------------------------------------------
bool output_config_strings =true;
ofstream configs_ti_t;	// ---- file stream to a t_index vs. tval for configstrings  
int Size_current_tvec;

if(calculate_entropy)
	{
	//-------------------------------------------------------------------

	for(i=0;i<total_obs_filling;i++)
		{ 
		pnorm=0.0;
		size_of_Zt = Z_all_t[i].Z_t.size() ;

		//----------- renormalize the probability -----------------------
		for(ii=0; ii<size_of_Zt; ii++)
			{
			pnorm += Z_all_t[i].Z_t.at(ii).pcount;
			}

		//------------calculate the entropy for this time point----------
 		for(ii=0; ii<size_of_Zt; ii++)
			{
			p = double(Z_all_t[i].Z_t.at(ii).pcount) / pnorm;
			
			Z_all_t[i].S -= p*gsl_sf_log(p);
			}
		//--------------average H for this time point  over all runs ------

		Z_all_t[i].H = (Z_all_t[i].H/float(numtrials)) ;

		//----------now add Nave*mu to get Htot for this time point ------

		Z_all_t[i].Nave = Z_all_t[i].Nave/float(numtrials);
		Z_all_t[i].Htot = Z_all_t[i].H - Z_all_t[i].Nave*muN;


		//--------- NOW SAVE THE MOST COMMON N CONFIGURATIONS PER TIME POINT ------
		get_first_N_pvals_from_config_list ( Z_all_t[i].Z_t, num_config_counts[i], get_first_N );
		//----------------------------------------------------------------
		}

	//---------------NOW OUTPUT THE LIST OF CONFIGURATION COUNTS--------------
	clear_charray(cpath, charlength );
	sprintf(cpath, "%sthermostuff_configuration_sampling_tcols.txt",pathout.c_str());
	
	fentropyout = new ofstream(cpath);

	for(j=0; j<get_first_N; j++)
		{

		for(i=0;i<total_obs_filling;i++)
			{ 
			*fentropyout << num_config_counts[i][j]  << " \t " ;
			}

		*fentropyout << endl ;
		}	
	(*fentropyout).close();
	delete fentropyout;

	//---------------NOW OUTPUT THE RESULTS OF THAT CALCULATION---------------	

	clear_charray(cpath, charlength );

	sprintf(cpath, "%sthermoquantities_t_S_H_Htot_Nave.txt",pathout.c_str());
	fentropyout = new ofstream(cpath);

	for(i=0;i<total_obs_filling;i++)
		{ 
		*fentropyout << Z_all_t[i].tval  << " \t " << (1.0/float(Llim))*Z_all_t[i].S << " \t " << (1.0/float(Llim))*Z_all_t[i].H << "\t" << (1.0/float(Llim))*Z_all_t[i].Htot << "\t" << (1.0/float(Llim))*Z_all_t[i].Nave << endl;
		//----normalized by Llim now to imply energetic quantities PER LATTICE SITE.
		}

	(*fentropyout).close();
	delete fentropyout;

	//------------OUTPUT CONFIG STRINGS?-----------------------------

	if( output_config_strings)
		{
		clear_charray(cpath, charlength );
		sprintf(cpath, "%sconfigstrings_ti_vs_t.txt",pathout.c_str());
		configs_ti_t.open(cpath);

		fentropyout = new ofstream;

		for(i=0;i<total_obs_filling;i++)
			{ 
			configs_ti_t << i << "\t" <<  tpoints_filling[i] << endl;
 
			clear_charray(cpath, charlength );
			sprintf(cpath, "%sconfigstrings_ti_%d.txt",pathout.c_str(),i);

			(*fentropyout).open(cpath);
			Size_current_tvec = Z_all_t[i].Z_t.size();

			for(j=0; j < Size_current_tvec ; j++)
				{
				*fentropyout << Z_all_t[i].Z_t.at(j).description << "\t" ;
				*fentropyout << Z_all_t[i].Z_t.at(j).pcount << endl;
				}
			(*fentropyout).close();
			
			}

		delete fentropyout; 

		}
	//--------------------------------------------------------------

	}
//======================   OUTPUT two-body potential =======================
bool output_v2=true;
ofstream *fVNNout;

if(output_v2)
	{
	clear_charray(cpath, charlength );
	if(!krm_b)
		{
		sprintf(cpath, "%sv2%s_fp-%d_E0-%.3f_krm_0.txt",pathout.c_str(),NGtype.c_str(), footprint, E0);
		}
		else
		{
		sprintf(cpath, "%sv2%s_fp-%d_E0-%.3f_krm_%lf.txt",pathout.c_str(),NGtype.c_str(),footprint,E0,krm_val);
		}
	fVNNout = new ofstream(cpath);

	for(x=0; x<NNRANGE; x++)
		{
		*fVNNout << xcoarse[x] << " \t " << VNN[x] << endl;
		}
	(*fVNNout).close();
	delete fVNNout;
	}
//======================   OUTPUT FIXED REF AVERAGING   =======================
// bool output_fixed_ave=false;
ofstream *favgout;

double CGF = 1.0;

if(output_patterns)	//---"patterns" refers to both 1-point and 2-point functions.
	{
	clear_charray(cpath, charlength );

	//--------------------FIRST THE EQUILIBRIUM VALUE ------------------
	sprintf(cpath, "%sFixed-occ-hist-eq_x_N_TF.txt",pathout.c_str() );
	favgout = new ofstream(cpath);

	for(x=0;x<Llim;x++)
		{
		*favgout << (x-(h[0] +h[1]))*CGF << " \t " << (OavN_eq[x])/CGF << " \t " << (OavTF_eq[x]/CGF) << endl;
		}
	(*favgout).close();
	delete favgout;
	//--------------------  NOW THE TRANSIENT VALUE --------------------
	
	clear_charray(cpath, charlength );
	sprintf(cpath, "%sFixed-occ-hist-ti_x_N.txt",pathout.c_str() );
	favgout = new ofstream(cpath);

	for(x=0;x<Llim;x++)
		{
		for( t_index=0;  t_index<total_obs_filling;  t_index++)
			{
			*favgout << (OavN_ti[t_index][x]/CGF) << "\t";
			}
		*favgout << endl;
		}
	(*favgout).close();
	delete favgout;
		
	//------now for the TF case------
	clear_charray(cpath, charlength );
	sprintf(cpath, "%sFixed-occ-hist-ti_x_TF.txt",pathout.c_str() );
	favgout = new ofstream(cpath);

	for(x=0;x<Llim;x++)
		{
		for( t_index=0;  t_index<total_obs_filling;  t_index++)
			{
			*favgout << (OavTF_ti[t_index][x]/CGF) << "\t";
			}
		*favgout << endl;
		}
	(*favgout).close();
	delete favgout;
	
	//------------------------------------------------------------------

	}

//======================   OUTPUT TWO-PART correlation  ======================
//output_patterns is defined above so that we can save on calculation time.

if(output_patterns)
	{

	//---------------   FIRST, EQUILIBRIUM DISTRIBUTIONS -------------
	clear_charray(cpath, charlength );
	if(!krm_b)
		sprintf(cpath, "%stwopartcorr_eq_%s.txt",pathout.c_str(),NGtype.c_str());
	else	
		sprintf(cpath, "%stwopartcorr_eq_%s.txt",pathout.c_str(),NGtype.c_str());

	favgout = new ofstream(cpath);
	for(x=0;x<Llim;x++)
		{
		*favgout << x*CGF << " \t " << (output2partcorr_eq[x]/CGF) << endl;
		}
	(*favgout).close();

	//---------------   NEXT, TRANSIENT DISTRIBUTIONS -------------

	clear_charray(cpath, charlength );
	if(!krm_b)
		sprintf(cpath, "%stwopartcorr_ti_%s.txt",pathout.c_str(),NGtype.c_str());
	else	
		sprintf(cpath, "%stwopartcorr_ti_%s.txt",pathout.c_str(),NGtype.c_str());

	(*favgout).open(cpath);


	for(x=0;x<Llim;x++)
		{
		for( t_index=0;  t_index<total_obs_filling;  t_index++)
			{
			*favgout << (output2partcorr_ti[t_index][x]/CGF) << "\t";
			}
		*favgout << endl;
		}

	(*favgout).close();
	delete favgout;
	}

//======================   OUTPUT FILLING DENSITY  ======================
bool output_filling=true;

if(output_filling)
	{
	clear_charray(cpath, charlength );

	sprintf(cpath, "%sfilling%s_t_v_rho_lines-are-tpoints-for-Vdistscoverage.txt",pathout.c_str(),NGtype.c_str());

	favgout = new ofstream(cpath);

	for(j=0;j<total_obs_filling;j++)
		{
		*favgout << tpoints_filling[j] << " \t " << fillingfrac[j] << endl;
		//--- NOTE: get_filling_frac() already scales out the filling frac by CGF upon computation. 
		//--- Hence, we  Don't need to do it again here.
		}
	(*favgout).close();
	delete favgout;
}
//======================   OUTPUT mean/std.dev STATISTICS  ======================
bool output_meanstddevvtime=false;

if(output_meanstddevvtime)
	{
	clear_charray(cpath, charlength );
	sprintf(cpath, "%svoidstats%sBZ%s_t_mean_stddev_Ncheck.txt",pathout.c_str(),NGtype.c_str(),BZ.c_str() );

	favgout = new ofstream(cpath);
	(*favgout).precision(10);

	for(j=0;j<total_obs_filling;j++)
		{
		*favgout << tpoints_filling[j] << " \t " << void_means[j] << " \t " << void_stddevs[j] << " \t " << Ncheck[j] << endl;
		}
	(*favgout).close();
	delete favgout;
}

//======================   OUTPUT void DISTRIBUTION STATISTICS AT EACH TIME POINT  ======================
bool output_vdist_timepoints=true;

double * Nvoidstot;
ofstream * favgout2;

if(output_vdist_timepoints)
	{

	Nvoidstot = new double[total_obs_filling];
	for (i=0;i<total_obs_filling;i++)
		{
		Nvoidstot[i] =0.0;
		}

	for (i=0;i<total_obs_filling;i++)
		{
		for(j=0;j<=Llim;j++)
			{
			Nvoidstot[i] += void_histogram[i][j]; //--'i' is the time point, so sum over 'j' (void sizes)
			}
		}


	//--------------------------------------
	clear_charray(cpath, charlength );
	sprintf(cpath, "%svoiddists%s_tps.txt",pathout.c_str(),NGtype.c_str());
	favgout = new ofstream(cpath);

	clear_charray(cpath, charlength );
	sprintf(cpath, "%svoiddens%s_tps.txt",pathout.c_str(),NGtype.c_str());
	favgout2 = new ofstream(cpath);
	//----------------------------------------

	for(j=0;j<Llim;j++)
		{
		//---we simply don't output the last "L" void value to the data file as it really just means an empty array.
		// we want to sum over these data to get the density of particle so to ensure one-to-one correspondence between number of voids counted and number of particles present, we omit the last void value from the file, which really indicates no particles at all.

		for (i=0;i<total_obs_filling;i++)
			{
			*favgout << (double(void_histogram[i][j]) /Nvoidstot[i])  << "\t";
			*favgout2 << (double(void_histogram[i][j]) / ( double(Llim*numtrials) ))  << "\t";
			}
		*favgout  << endl;	//----produces columns of void sizes with common time points ('i') at each. 	
		*favgout2 << endl;	//----the rows are then constant void size at increasing t. 
		}


	(*favgout).close();
	delete favgout;

	(*favgout2).close();
	delete favgout2;

	delete [] Nvoidstot;
	}


//===============   OUTPUT void DISTRIBUTION STATISTICS across equilibrium  ================
double Nvoidstotal_eq=0.0;

if (get_voiddist_equilibrium)
	{
	for(j=0;j<=Llim;j++)
		{
		Nvoidstotal_eq += void_histogram_equilibrium[j]; //--it's averaged over all time, so only a 1-D array.
		}
	}

if (get_voiddist_equilibrium && Nvoidstotal_eq >=1 )
	{
	clear_charray(cpath, charlength );
	sprintf(cpath, "%svoideq_dens_dist.txt",pathout.c_str() );
 
	favgout = new ofstream(cpath);


	for(j=0;j<=Llim;j++)
		{
		*favgout << fillingfrac[total_obs_filling-1]*(double(void_histogram_equilibrium[j]) /Nvoidstotal_eq) <<  " \t " << (double(void_histogram_equilibrium[j]) /Nvoidstotal_eq)  << endl;
		}

	delete favgout;
	}
//======================   OUTPUT OVERSHOOT DATA     ============================
bool output_overshoot_phase_diag=true;
//---try to use the explicitly equilibrium averages later.

double maxrho  = get_array_max( NULL, fillingfrac , total_obs_filling);
double rho_eq  = get_local_array_avg( fillingfrac, (total_obs_filling-1) -1 , 1); //--averages over the last 3 entries 

// double rho_eq2 = (double(Nvoidstotal_eq)/double(total_obs_eq));

//--- this will take the last 5 data points averaged. 

if(output_overshoot_phase_diag)
	{
	clear_charray(cpath, charlength );
	sprintf(cpath, "%sovershootdat-%sBZ%s_muN-%.2f_E0-%.2f_maxrho_finalrho.txt",pathout.c_str(),NGtype.c_str(),BZ.c_str(), muN, E0 );

	favgout = new ofstream(cpath);

	

	*favgout << muN << " \t " << E0 << " \t " << maxrho << " \t " << rho_eq << endl;


	(*favgout).close();
	delete favgout;
}


//======================  --  HOUSE KEEPING  --  ===================================================

//-------PLANT THE SEED FOR THE NEXT RUN-----------

seed=gsl_rng_get(r);

/*------- DONT REPLANT SEED FOR BATCH JOBS--------
clear_charray(cpath, charlength );
sprintf(cpath, "%srngSEED.in",path.c_str());
ofstream fseedout(cpath);

 fseedout.precision(18);
 fseedout << seed;
 fseedout.close();
----------------------------------*/

//----------------- CLEAR UP ALLOCATED MEMORY -------------------------------

delete [] fillingfrac;

for (i=0;i<total_obs_filling;i++)
	{
	delete [] void_histogram[i];
	}
delete [] void_histogram;

if (get_voiddist_equilibrium)
	{
	delete [] void_histogram_equilibrium; //--only a 1-D array
	}

if( tf > t_trans + dt_obs)
	{
	delete [] tpoints_eq;
	}

delete []  void_means     ;
delete []  void_stddevs   ;
delete []  Ncheck         ;

gsl_rng_free (r);


for(t_index=0; t_index<total_obs_filling; t_index++)
	{
	delete [] num_config_counts[t_index];
	delete [] output2partcorr_ti[t_index];

	delete [] OavN_ti[t_index];
	delete [] OavTF_ti[t_index];

	}



//====================================================================================================

if(should_plot_snapshots )
	{	
	(*timestamps).close();
	}

if(should_plot_kymo)
	{
	(*fkymo_outN).close();
	(*fkymo_outTF).close();
	}

//(*fout).close();
//delete fout;

//---------------    FINISH AND EXIT   -----------------------------------

*log  << "\n measurement successfully completed with the following parameters: \n\n" ;

*log  <<  "\n average num_binds_per_simulation_F0 = " << num_binds_per_simulation_F0 ;
*log  <<  "\n average num_binds_per_simulation_F1 = " << num_binds_per_simulation_F1 << endl ;



*log  << "\n numtrials = " << numtrials << endl ;
(*log).close();
delete log;

cout << "\n\n program completed successfully.\n ";
cout << " the asymptotic average of the central third of the system was: " << avg_as << endl;
return 0;
}

