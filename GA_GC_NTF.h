/*  -----GA_GC_NTF.h - defines the classes and structs and declares some of the functions for the GC GA
// ---last updated on  Tue Feb 18 18:27:29 CET 2014  by  Brendan.Osberg  at location  th-ws-e537

//  changes from  Tue Feb 18 18:27:29 CET 2014 : implemented simplified run for small particles. Tested, works.

//  changes from  Thu Jan 9 12:42:40 CET 2014 : additional collection series for the 2-point correlation function: gathering curves during transient filling process now in addition to during equilibrium.


//  STARTED FROM SCRATCH BUILDING CODE FOR SMALL PARTICLES -SOME FEATURES WERE KEPT, SOME WERE AXED.
-------------------------------------------------------------------------------------------------*/


#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>

using namespace std;


const bool TFs_allowed     = false; 	//--are TFs allowed in this simulation?
const bool debugging       = false; 	//-- are we just doing a dummy run for debugging purposes?

const bool bind_irrev      = false;	// do particles bind irreversibly? if so, then k_off always = 0.

const bool choose_carefully = true;	//do we sum over the rates again every time to check the an/a0 ratio?

const bool get_voiddist_equilibrium = true;


const bool fixed_ref       = false;	//--Is the +1 nucleosome hard-coded at its position?
					// BELOW: ARE THERE *OTHERS* that are also fixed in place?

const bool set_fixed_initial = false;   // -- if true, then we set an initial set of particles 
					//  to be fixed in place without the possibility of removal,
					//  other particles coming after that bind reversibly.

const bool calculate_entropy = false; 	// if true, then count configurations at every  
					// time point, for every iteration. 

/********************************************************************/
struct bindevent
{
double t; 	// the time that this occured.
int   species; 	// the type of binding object, 1=Nucl, 2=TF
int   onoroff; 	// +1 = binding, -1 = unbinding.

};
/********************************************************************/
struct site
{
int state;		// the current occupation state of this site: 0=empty,
			// 1=nucleosome, 1..m=positions along TF

bool permanent;		// if (1) -then at this site, the particles never come off.

int part_right; 	// the location of the nearest particle to the right
int part_left;		// the location of the nearest particle to the left


double * a_removeN;
double * a_addN;	// delta-Energy of removal from system
double * a_sNL;		// delta-Energy of sliding left one position
double * a_sNR;		//	"	"	"	" right "	"

double * a_removeTF;
double * a_addTF;
double * a_sTFL;			//delta-Energy of sliding left one position
double * a_sTFR;			// "	"	"	" right "	"
//---- ALL OF THE ABOVE POINTERS ARE JUST SET TO POINT TO MEMORY 
//---- ALREADY ALLOCATED FOR THE Rx array BY THE GAdata CONSTRUCTOR THEY DO NOT NEED TO BE DELETED.
//---- THESE POINTERS ARE NEVER ALLOCATED MEMORY IN THEIR OWN RIGHT.


double t_occ_TF;		// the TOTAL time that this site has been
double t_occ_N;			// occupied by a TF, Nuc, respectively.

double t_bind_TF;		// the time at which the TF/Nuc currently there
double t_bind_N;		// binded

queue<bindevent> event_list;
// n.b. we eliminated an array of binding events to avoid confusion.

};
//***********************************************************************
struct particle
{
int type;
int width;

site * location;
particle *Right;
particle *Left;
};
//***********************************************************************
class configuration
{
public:
string description;	//-----characterizes one full "state" of the system.
			//---- this string captures the whole configuration, 
			//---- eg 26-27-100 has those three voids in that order.

int pcount;		//--- the actual number of observations of this configuration.

//---------------------- THINGS WE DECIDED WE DON'T NEED ANYMORE -------------
// float p;		//--- the total probability of this configuration.
// vector<int>  dists; 	//--- size of N -one int for each N denoting its position.
// int N;		//--- the number of particles in the system.
			//--the number of voids is then also N
//----------------------------------------------------------------------------
};
//***************************************************************************
class config_set_t 
{ //----this is the set of configurations at a particular time 't' 
  //--- also includes thermodynamic data at that time.

public:

double tval;
bool   tpoint_passed;  // --- array[i] indicates whether 'i'th t-point has been considered yet.
float    Nave;

double   S;		// entropy

double   H;		// here is JUST the neighbour interaction, i.e. 
			// sum of 'v' between neighbours.
double Htot;		// here is the above H, PLUS, <N> times mu

vector < configuration >   Z_t;	//----Z_t is always just the set of configurations.

};
//**************************************************************************
class GAdata
{

private:
bool HNG; 	// are we dealing with Hard core particles? 
bool SNG;	// Soft-core?
bool LNG;	// Linear potential?


bool boltzmann_on_uphill ;
bool boltzmann_on_add    ;
bool boltzmann_on_removal; 
bool boltzmann_on_addrem_intermed;

double BZalpha;
//-----this condition is now read in from file.
int h0, h1, h2;
int F0, F1;  // F0, F1 are now the locations of the two transcription factors directly. 
			//--see figure in main program.

double muN0;
double *muNarray; //----the array of ALL values of the chemical potential.
double muTF0;
int Llim;
int M;
int footprint;

double *VNN;		//---effective VNN, after coarse-graining.
double *VNTF;

double *xcoarse;

public: 
bool krm_b;	// 1- remodellers are present
		// 0- "    "    "   not "
double krm_val; // the actual value of the remodelling complex 'k'.
int    RMrange; // the range over which the remodeller is able to reach.
double N_remod; //-  number of pairs of nucleosomes that are adjacent, and 
		//-  close enough together to be remodelled.

int NNRANGE;	//----the total range of the N-N interaction
int NTFRANGE;	//----" 	"	"    N-TF interaction

int min_dist;	//---the minimum possible distance between adjacent particles. This is 1 for LNG and SNG, but 147 for HNG.

double CGF;

double* tpoints_filling; //--- these are time points used to take the filling fraction through the transient period 
			 //--- (should usually by logarithmically spaced.)
double* tpoints_eq;	 //--- these are time points used to time-average equilibrium values (should usually be linearly spaced)


int flag;
int counter; // records the iteration we're currently on.

int recentRx;

//------------------------------ SETTING UP THE GAdata	INITIALLY: ------------------------
GAdata(const double *energetics, const int* observations, const double* k_rates,const double* times,const int* sizes_n_ranges, const int *h, const int * F, bool krm_B, string path, ofstream * log_in , ofstream * timestamps_in, const bool * flags);	//-- the constructor

int set_fixed_initial_particles( const string path , gsl_rng * r); //-- possible function to implement some initially fixed particles.

~GAdata();
//-----------------------------------------------------------------------------------------


double muN(int x );
double muTF(int x );

int distance(const int a, const int b ); 	//---number of spaces away, not the number of intervening sites.
int remod_cand(const int a, const int b);	// returns either 1 or 0, to denote whether or not this pair is a candidate for remodelling.	

int left(const int x );
int right(const int x );

long double a0;
double * Rx;	// 0-8 indicates the most recent reaction undertaken.
// 0 - remove_Nuc
// 1 - add_Nuc
// 2 - slide_Nuc_left
// 3 - slide_Nuc_right
// 4 - remove_TF
// 5 - add_TF
// 6 - slide_TF_left
// 7 - slide_TF_right

site * pos;

int choose_reaction(double ranvar);

double interaction_dEadd(int type, int x, int *neighbours   );
double interaction_NN(int x);
double interaction_NTF( int x);
double interaction_TFTF( int x);

double dEsL(int type, int x, int * neighbours );
double dEsR(int type, int x, int * neighbours );

double k_E(const double dE);
//---------------the actual reaction functions----------------

int remove_Nuc( int x);
int add_Nuc(  int x);
int slide_Nuc_left( int x);
int slide_Nuc_right( int x);

int remove_TF(  int x);
int add_TF(  int x);
int slide_TF_left(  int x);
int slide_TF_right( int x);

int remodel(int R, gsl_rng * r);

int calc_rates( const int min, const int max, const int n);
int get_HNG_rates(  const int type, double * HNG_outputrates, const int * neighbours, const int x ); // subroutine called in cases of HNG

void process_error( int x);
void interaction_error(int type, int x, int *neighbours );
int check_rates( void );
int check_states( void );
int reset_prevnext_pointers();

bool already_warned;

//------------ the observables -------------------------------

	int increment_void_histogram( int ** void_hist);
	int increment_void_histogram_equilibrium(int* void_hist);

//----FILLING RATES-----------

	double * filling_frac;	// the array of filling fractions at various time points
	int get_filling_frac(void);

	bool should_observe_filling(double tau, int ** void_histogram, config_set_t * Z_all_t );
	int  should_observe_equilibrium_voiddist( const double tau, int * void_histogram_equilibrium  ); 
	int should_observe_patterns(double tau);

	int grab_current_configuration( config_set_t & C_t);

//----------------------------
	int increment_2_part_corr( int timepoint ); // done many times during the run.
	int normalize_2_part_corr(void );	// done at the end of each run.

	double *  two_part_corr_eq;		// two-particle correlation at equilibrium.
	int       num_times_2pc_eq_incremented; // which has been incremented this many times.
	double ** two_part_corr_ti;		// two-particle correlation per 't'.
	int    *  num_times_2pc_ti_incremented; // each of which have been 
						// incremented this many times.
	//-------
	int increment_1_point_hist( int timepoint );
	int * onepoint_histocc_N_eq;
	int ** onepoint_histocc_N_ti;

	int * onepoint_histocc_TF_eq;
	int ** onepoint_histocc_TF_ti;	//==== fixed-average histograms of N/TF positions.
	//-------


	int printout_states( ofstream *fout);
	int printout_avgs( ofstream *fout);

	//-------  TF time-averaged occupancy stuff   ----------
	int bindeventnum_F0;
        int bindeventnum_F1;

        bindevent * event_array_F0;
        bindevent * event_array_F1;

        bool built_array_of_binding_events; //--have we or have we not constructed such an array from the stack?

	int obs_count_filling;	// the current register of observation snapshot that have been taken.

	double avg_F0_occupation;
	double avg_F1_occupation;

//-------DELETE THIS ----------
	int observe_TF_occ_lag( void );
	int get_site_tpoint_occ( const double t  , const int site);

	double testing_avg_F0_occupation;
	double testing_avg_F1_occupation;

	double * testing_olap_y;
	double dtau_obs;
//------DOWN TO HERE ----------

	int get_tarray_from_stack(void);

//------DELETE THIS:
	int get_delayed_tcorr(const double * time_corr_t, double *time_corr_y, double *raw_overlap, const int nbins);
//-----DOWN TO HERE.

	double get_overlap(const double  tau );
	double get_frac_occupied(const double  t1, const double t2 , const int site);

	double increment_beginning_at( const int i, const int site );
	double increment_ending_at( const int i, const int site );

//------variables to do with printing out the current density------
int Nplots2makeshort, Nplots2makelong;
int Nplots2make_kymo;

int nbins; // ---the number of discrete points we consider in time-lag correlation.

// bool printtime(void);
bool printtime_kymo(void);

int  plot_snapshot(void);
int  plot_snapshot_kymo( ofstream * foutN, ofstream * foutTF );

int  plotnum; //---the current plot number
int  plotnum_kymo;

string path;
//-------------------------------------------------------------

int Nucnum, initialized_Nucnum;
int TFnum;
int partnum;	// total number of particles.

double t;  	//--- the current time.
double tf; 	//--- the termination (final) time
double t_trans; //--- the transient time before we start averaging on steady state.

double ka_N,  ks_N;
double ka_TF, ks_TF;
//! int w, m; -----don't need w anymore!
int m;
double E0;

double dt_obs;	// the amount of time elapsing in between snapshots (assuming linear spacing).
int obs_count_eq_2pc;	// the current register of 2pc observations that have been taken.
int obs_count_2pc_ti;
int obs_count_eq_vdist;	// the current register of eq. void dist. observation snapshot that have been taken.


int total_obs_filling;	// the total number of observation snapshots that will be taken.
int total_obs_eq;	// the total number of observation snapshots that will be taken.


ofstream * log;
ofstream * timestamps;


int cleanup_occ( void);		// TIDYING UP MEMORY
int reset( void);
};
//*********************************************************************************

int VNTF_calc(double * potential,const int w, const double E0);

double interact_NN(double * VNN, int w, int x);
double interact_NTF(double * VNTF, int w, int x);

double k_E(const double dE);
//--------------------------------------------------------------

int react(int x,int Rtype, GAdata &P);


//--------------observables--------
int space_tpoints_linear(const double ti, const double dt_obs, const double tf, const int nbins, double * tx);
//================================

double min(const double a, const double b);
double max(const double a, const double b);


