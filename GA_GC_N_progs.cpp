/*----- GA_GC_NTF_progs.cpp: functions/subroutines for the grand-canonical nucleosome gillespie algorithm.-----
// ---last updated on  Tue Nov 4 18:40:24 CET 2014  by  ga79moz  at location  TUM , xenopus

//  changes from  Tue Nov 4 18:40:24 CET 2014 : resolved some min/max wrap-around issues for long-range interactions. Also updated the SNG potential for Lennard-Jones

//  changes from  Fri Oct 17 15:49:24 CEST 2014 : made this small_particles branch the master branch in Git. Git is now used as the change-tracker, check the history of the git repository from this point forward.

//  changes from  Tue Jul 1 16:59:36 CEST 2014 : set output folder to just local/muN-..., and added an output file for the void _density_ in addition to distribution

//  changes from  Wed Apr 23 15:32:48 CEST 2014 : removed VNN_S/LNG_calc_smallp definition from this file -now the only definition is in bren_lib.h

//  changes from  Tue Feb 18 18:27:29 CET 2014 : implemented simplified run for small particles. Tested, works.

//  changes from  Thu Jan 9 12:42:40 CET 2014 : additional collection series for the 2-point correlation function: gathering curves during transient filling process now in addition to during equilibrium.

//  STARTED FROM SCRATCH BUILDING CODE FOR SMALL PARTICLES -SOME FEATURES WERE KEPT, SOME WERE AXED.

===========================================================================================
-----------------------------------------------------------------------------------------*/

#include <fstream>   
#include <iostream>  
#include <iomanip>  
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <queue>
#include <math.h>
#include "GA_GC_N.h"
//#include <bren_lib.h>

#include <gsl/gsl_rng.h>
int VNN_LNG_calc_smallp(double * potential, const int a, const double E0);
int VNN_SNG_calc_smallp(double * potential, const int NNRANGE, const int a, const double rm, const double E0);
//---the SNG func has an extra input parameter for the LJ potential.

int coarse_grain(const double * Vin_full, double * xcoarse, double * Vout_coarse, const int L, const int p, const double CGF); //--coarse-grain the 2-body interaction potential into a smaller system.

//-------these two functions are from bren_lib, but are getting declared again, because we can only inlude the library once
int clear_charray(char* c, const int charlength );
int pick_index_from_rand_norm_dist( double rand01, double * distr, int L );
double interpolate_general( const double * x1, const double * y1, const int size_in, const double * x2, double * y2_out, const int size_out);

bool compare_dists ( const vector<int> A, const vector<int> B , const int N);
//----take two vectors of ints and see whether the entire list is the same or not.


int get_first_N_pvals_from_config_list( const vector < configuration >   Z_t, float * target, const int N);
// ---order the configurations along 'p', and then put the first N values of that 
// ---ordered list into the target array


string bren_itoa( const int x );

//-------------------------------------------------------------------------------------------------------------------------

//*******************************************************************
//-------------------------- THE BIG CONSTRUCTOR --------------------
GAdata::GAdata(const double *energetics, const int* observations, const double* k_rates,const double* times, const int* sizes_n_ranges, const int *h, bool krm_B, string pathin, ofstream * log_in , ofstream * timestamps_in, const bool * flags)  //constructor
{  // constructor
int x,i,j;

//-----------SIZES AND RANGES-----------
   footprint      = sizes_n_ranges[0];
 	//---   kHNG        *is now just 'a'*      = sizes_n_ranges[1];
 	//the range in units of uncoarsened lattice sites.
   RMrange            = sizes_n_ranges[1];
   Llim               = sizes_n_ranges[2];
   NNRANGE	      = sizes_n_ranges[3];	
 
   NTFRANGE=0;

   M    =  reactions_per_site*Llim;	// the number of possible reactions -- most of which will have 0 amplitude

//------------RATES---------------------
   ks_N       =  k_rates[0];
   ka_N       =  k_rates[1];
   krm_val    =  k_rates[2];

//-------------TIMES--------------------
   tf       = times[0];
   dt_obs   = times[1];
   t_trans  = times[2];
//------------FlAGS---------------------
   HNG = flags[0];
   SNG = flags[1];
   LNG = flags[2];
   boltzmann_on_uphill          = flags[3] ;
   boltzmann_on_add             = flags[4]   ;
   boltzmann_on_removal         = flags[5]; //-----this condition is now read in from file.
   boltzmann_on_addrem_intermed = flags[6]; //-----this condition is now read in from file.
//------------OBSERVATIONS--------------
   Nplots2makeshort   = observations[0];
   Nplots2makelong    = observations[1];
   Nplots2make_kymo   = observations[2] ;
   total_obs_filling  = observations[3];	//----total number of times for transient observations
   total_obs_eq       = observations[4];	//---- ""	""	""	equilibrium 	""
   nbins              = observations[5];
//------------ENERGETICS----------------
   muN0          =  energetics[0];   //---this is the genuine muN 

   E0            =  energetics[1];
   BZalpha       =  energetics[2];   //--- the boltzmann alpha value 
				     //--- that defines the ratio 
				     //--- between addition and removal.

//--------------------------------------
if( int(boltzmann_on_uphill) + int(  boltzmann_on_add ) + int(boltzmann_on_removal) + int(boltzmann_on_addrem_intermed)  != 1)
	{
	cout    << "\n ERROR: sum of bolzmann conditions should be 1 (i.e. one should be 1, all others zero.) somethings wrong here, exiting. \n";
	*log_in << "\n ERROR: sum of bolzmann conditions should be 1 (i.e. one should be 1, all others zero.) somethings wrong here, exiting. \n";
	(*log_in).close();
	exit(1);
	}

//---------------------THESE ARE PARAMETERS THAT GET USED ALL THE TIME ----------------------------

counter = 0;
flag = 0;
CGF=1.0;

//---the number of times that we have made observations of :--------------------------
obs_count_eq_2pc	=0;	// ---- the equilibrium 2 particle correlation funciton.
obs_count_2pc_ti	=0;	// ---- the transient 2 particle correlation function.
obs_count_eq_vdist	=0;	// ---- the equilibrium void distribution
//------------------------------------------------------------------------------------



plotnum_kymo=0;
Nucnum  = 0;
partnum = 0;
a0=0.0;
t=0.0;

N_remod=0;  //----number of possible available remodelling sites.

num_times_2pc_eq_incremented=0;
two_part_corr_eq 	= new double [Llim];
onepoint_histocc_N_eq	= new int [Llim];

for(i=0;i<Llim;i++)
	{
	two_part_corr_eq[i]       = 0.0;
	onepoint_histocc_N_eq[i]  = 0;
	}


two_part_corr_ti             = new double * [total_obs_filling];
num_times_2pc_ti_incremented = new int    [total_obs_filling];

onepoint_histocc_N_ti	= new int* [total_obs_filling];

for(i=0;i<total_obs_filling;i++)
	{
	num_times_2pc_ti_incremented[i]  = 0.0;
	two_part_corr_ti[i]              = new double [Llim];

	onepoint_histocc_N_ti[i]         = new int[Llim];

	for(j=0;j<Llim;j++)
		{
		two_part_corr_ti[i][j]       = 0.0;
		onepoint_histocc_N_ti[i][j]  = 0;
		}
	}

//----NB: the twpoints_filling/eq arrays are assigned by pointer value in the main function!

already_warned = false;


//--------------  THESE ARE ALL FILLING DYNAMICS PARAMETERS ---------------------

filling_frac = new double[total_obs_filling];
for(i=0; i<total_obs_filling; i++)
	{
	filling_frac[i]=0.0;
	}

obs_count_filling=0;	//--- the index of the time-snapshot that is 
			//--- currently being considered for filling evaluation.


//----------------------------------------

   recentRx = -12;
   krm_b    = krm_B;

//------------------- NOW GET THE POTENTIALS ----------------------------------

double * VNN_full = new double[NNRANGE+1];
xcoarse           = new double[Llim]; for(i=0;i<Llim;i++){xcoarse[i]=0.0;}
VNTF=NULL;

if(SNG)
	{
	min_dist=1;	//---minimum distance that _can_ exist between 2 adjacent particles.
			//---in this case it's 1 (in principle), but that close is energetically costly.
	VNN_SNG_calc_smallp(VNN_full, NNRANGE, footprint, VLJ_rm, E0);
	}
else if(LNG)
	{
	min_dist=1;
	VNN_LNG_calc_smallp(VNN_full, footprint, E0);
	}
else
	{
	// must be HNG
	if (!HNG)
		{
		cout    << "\n ERROR: the NGtype is not defined.\n";
		*log_in << "\n ERRN_OR: the NGtype is not defined.\n";
		(*log_in).close();
		exit(1);
		}

	min_dist= footprint;	//--  The minimum possibile distance between adjacent particles. 
				//--  Prevent remodellers from violating this. 

	for(i=0;i<=NNRANGE;i++)
		{
		if(i<(footprint-1)) //---the coarse-graining comes later.
			{
			VNN_full[i]=(1/0.0);
			}
		else	
			{
			VNN_full[i]=0.0;
			}
		}

	}
//------- the last entry in this array should be zero.

//------- whatever the type, we now have the potential. the next step (coarse-graining) 
//------- is independent of type.

VNN       = new double[Llim]; //--most of these will just be zero, and never get used.
for(i=0;i<Llim;i++)
	{
	if (i< NNRANGE )
		{
		VNN[i] = VNN_full[i];
		}
	else
		{
		VNN[i] = 0.0;
		}
	}
delete [] VNN_full;


//---------------------------------------------------------------------------------

path    = pathin;
plotnum = 0;	// the "current" plot starts at zero, naturally.

log = log_in;
timestamps=timestamps_in;


//------------------ NOW SET YOUR POINTERS PROPERLY ----------------------------

Rx = new double[reactions_per_site*Llim]; //---changed this to 4
pos = new site[Llim];
muNarray = new double[Llim]; ///------for now we just set this to muN0

for(x=0;x<Llim;x++)
	{
	pos[x].a_removeN = &Rx[x*reactions_per_site + 0 ];
	pos[x].a_addN    = &Rx[x*reactions_per_site + 1 ];
	pos[x].a_sNL     = &Rx[x*reactions_per_site + 2 ];
	pos[x].a_sNR     = &Rx[x*reactions_per_site + 3 ];

	pos[x].permanent  = false;	//all sites should be set to allow for removal of particles once there.

	muNarray[x] = muN0;
	}


//------HERE WE SET THE INITIAL REACTION RATES----------------
for(x=0;x<Llim;x++)	
	{
	
	pos[x].state=0;		// initialize each lattice pos to 0 (nothing)
	pos[x].part_right = x;
	pos[x].part_left  = x;


	*pos[x].a_removeN = 0.0;
	*pos[x].a_sNL = 0.0;
	*pos[x].a_sNR = 0.0;

	//----------nucleosome addition rate:
	*pos[x].a_addN  = ka_N*k_E( -muN(x) ); 	//--- ADDING NUCLEOSOMES (always weighted from initially empty)  -------------

	a0+=*pos[x].a_addN;

	//-----------------AND SET AVERAGE OCCUPANCY TO ZERO INITIALLY------------
	pos[x].t_occ_N  = 0.0;
	pos[x].t_bind_N  = -1.0;

	}//---end for loop initializing through the array of positions.
	

//---------ALWAYS EXPORT AS EMPTY FROM HERE, IF (set_fixed_initial) is true, then we'll add the fixed ones in main.

}
//*******************************************************************
int GAdata::set_fixed_initial_particles( const string path , gsl_rng * r )
{ //---- takes an equilibrium distribution and assigns some initial, permanent particles that never leave.

int charlength=400; //--the number of characters in the string for our path.
char cpath[charlength];
ifstream fin_eq;

clear_charray(cpath, charlength );
// sprintf( cpath, "%svoid_dist_fixed_seed.txt",path.c_str());
 sprintf( cpath, "%svoid_dist_fixed_seed_RSAdiff.txt",path.c_str());

// sprintf( cpath, "%svoid_dist_fixed_seed_postRSA.txt",path.c_str());


fin_eq.open(cpath);

if( fin_eq.fail() ) 
	{
	*log << "\n ERROR in set_fixed_initial_particles : FAILED to load equilibrium distribution file.\n\n";
	cout << "\n ERROR in set_fixed_initial_particles : FAILED to load equilibrium distribution file.\n\n";
	(*log).close();
	exit(1);	
	}

//---------get length l------------------
int l=0;
double dummy=0.0;
fin_eq >> dummy;

while( !fin_eq.eof() )
	{
	l++;
	fin_eq >> dummy;
	}
fin_eq.close();
//---------------------------------------
double eqdist_full[l];	//----the array of probabilities.
double x_full[l];

double eqdist_CG[Llim];

int i=0;
fin_eq.open(cpath);
while( !fin_eq.eof() )
	{
	x_full[i]  =  i;
	fin_eq    >> eqdist_full[i];

	i++;
	}

fin_eq.close();

//-----------------------NOW COARSE-GRAIN THE PROBABILITY DISTRIBUTION---------

double eq_CG_norm = interpolate_general( x_full, eqdist_full, l, xcoarse,  eqdist_CG, Llim);
for(i=0;i<Llim;i++)
	{
	eqdist_CG[i] = eqdist_CG[i]/eq_CG_norm; //---normalize it back to unity.
	}

//--------------------now we've got our distribution fill up the array -------

double rand01=0.0;
int binddist=0;
int x=0;
i=0;

while(1)	
	{

	pos[x].permanent  = true;	//this must be done first so that the off rate is set to zero.
	add_Nuc(x);

	i++;

	rand01=gsl_rng_uniform(r);
	binddist = 1 + pick_index_from_rand_norm_dist( rand01, eqdist_CG, Llim );

	x = x+binddist;

	if ( x >=Llim )
		{
		break;
		}	
	//----else, keep going on and fill in the next 'x'
	}


initialized_Nucnum = Nucnum;
}
//*******************************************************************
GAdata::~GAdata()  //destructor
{
int x=0,i=0,j=0;

delete [] VNN;
delete [] Rx;

	
delete [] muNarray;


//	delete [] filling_frac;
delete [] xcoarse ;
delete [] filling_frac;


//-------------------------------
delete [] two_part_corr_eq;
delete [] onepoint_histocc_N_eq;

for(i=0;i<total_obs_filling;i++)
	{
	delete [] two_part_corr_ti[i];
	delete [] onepoint_histocc_N_ti[i];

	}
delete [] two_part_corr_ti;
delete [] num_times_2pc_ti_incremented;

delete [] onepoint_histocc_N_ti;

//-------------------------------

}

//************************************************************************************
//------     HERE ARE THE DYNAMIC FUNCTIONS -  MOVING NUCLEOSOMES AROUND     --------
//************************************************************************************

int GAdata::choose_reaction(double ranvar)
{
int i=0;
int j=0;
double an=0.0;
bool found=false;
bool remod=false;
int result;

//----------------ERROR CHECK --- DELETE THIS ---
static int choose_carefully_counter = 0;
if (choose_carefully) {
	if (choose_carefully_counter >= choose_carefully_checking_interval) {
		for (i=0; i< M ;i++)
			an+=Rx[i]; 
		an += 2.0 * N_remod * krm_val;

		double local_errtest = (an-a0)/a0;

		if (fabs(local_errtest) > 1E-10) {
			if (!already_warned) {
				// *log << "\n WARNING, reaction rate sums don't agree\n recalibrating \n";
				// cout << "\n WARNING, reaction rate sums don't agree\n recalibrating at counter = " << counter << endl;
				already_warned = true;
			}
	
			if(fabs(local_errtest) > 1E-7) {
				cout << "\n WARNING: significant discrepency discovered. \n";
				*log << "\n WARNING: significant discrepency discovered. \n";
				// exit(1);
			}
			check_rates();
		}
		choose_carefully_counter = 0;
	} 
	choose_carefully_counter++;
}
an=0.0;
//-------------------------------DOWN TO HERE ----

for(i=0; i< M ;i++)
	{
	an+=Rx[i]; 
	if(an > ranvar*a0) 
		{
		found = true;
		break;
		}
	}


if( !found && krm_b) //--error check here.
  {

   double errtest = fabs( ( (an+(2.0*N_remod*krm_val))-a0)/a0); 

   if( errtest >1E-7 )
	{
 	flag=103112;
	process_error(0);
	cout << "\n exiting because errtest>1E-7"  << endl;
	*log << "\n exiting because errtest>1E-7"  << endl;
	exit(1);
	}
  }


else if( !found && (!krm_b ||i!=M) )
  {
 	flag=163071;
	recentRx=-11;
	cout << "\n exiting upon ranvar=" << ranvar << endl;
	*log << "\n exiting upon ranvar=" << ranvar << endl;
	(*log).close();

	process_error(1);
	exit(1);
  }



if(!found && i==M && krm_b )
	{
	for(j=0; j< N_remod ;j++)
	    {
	    an+= 2.0* krm_val;
		if(an > ranvar*a0) 
		{
		found = true;
		remod = true;
		break;
		}
	    }
	}


if(!found)
	{
	if( fabs((an-a0)/a0) >1E-12 )
	   {
 	   flag=103112;
	   recentRx=-11;
	   process_error(0);
	   exit(1);
	   }
	else
	   {
	   j=N_remod-1;
	   found=true;
	   }

	}
if(found)
	{
	result=i+j;
	}
else
  {
 	flag=12357;
	recentRx=-11;
	process_error(3);
	exit(1);
  }

	return result;
}

//***************************************************************************************
int  GAdata::remodel(int R, gsl_rng * r)
{
int Rtype;
int x=0;
int first_part, last_part;
int pL,pR;
bool found=false;

if(R < M )	// check if the reaction index 'R' is less than the total number of possible reactions.
   {
   flag      =995;
   recentRx  =-20;
   process_error(x);	
   }
else
   {	
   Rtype = R-M;
   }

first_part = pos[Llim-1].part_right;
last_part  = pos[0].part_left;  // either way.

x=first_part;
int n    = 0;
int Rpos = 0;

double r1;

if(Nucnum < 2)
   {
   flag=1295;
   recentRx=-20;
   process_error(x);	
   exit(1);
   }

while(1)
  {
   pR = pos[x].part_right;
   if( pos[x].state ==1  && pos[pR].state ==1 && remod_cand(x,pR) )	//---ARE THEY BOTH NUCL'S AND ARE THEY CANDIDATES FOR REMODELLING?
	{
	if(n==Rtype && !found)
	   {
	   found=true;
	   Rpos=x;
	   }

	n++;//always increment here.
	}
   if(pR == first_part) 
	{
	if(!found) //--if we've come around and not found anything
	   {
	   flag=12365; recentRx=-224; process_error(x); exit(1); 
	   }
	else
	   break;
	}
   else
	{
	x=pR; //keep going.
	}

  }

if(n!=N_remod)	//----  having counted all the way around, n should now be identical to the number of 
   {		//----  pairs that are eligible for remodelling.
   
   flag=105;
   recentRx=-20;
   process_error(x);	
   exit(1);
   }


x  = Rpos;
pR = pos[x].part_right;

//--------don't remodel away the +1 Nucl. if fixed_ref is true-----------
//---------down to here---------

r1=gsl_rng_uniform(r); 
if(r1 <0.5) 
	slide_Nuc_right( x); 
else 
	slide_Nuc_left( pR);


if( N_remod < 0 )
	{
	cout << "\n ERROR: N =" << N_remod << " at t =" << t << endl;
	*log << "\n ERROR: N =" << N_remod << " at t =" << t << endl;
	(*log).close();
	exit(1);
	}

return n;
}
//************************************************************************************
int GAdata::remove_Nuc( int x )
{
if(bind_irrev || pos[x].permanent)
	{cout << "\n ERROR, binding is irreversible, and yet remove_Nucl is being called. exiting.\n";
	*log << "\n ERROR, binding is irreversible, and yet remove_Nucl is being called. exiting.\n";
	(*log).close();
	exit(1);
	}


int pR = pos[x].part_right;
int pL = pos[x].part_left;
int i, min, max;

bindevent Q;  Q.t = t;  Q.species = 1;  Q.onoroff = -1;

//------------------REMODELLING SITE INCREMENT-------------------------
int n=0;

if( remod_cand(x,pR)  && pos[pR].state == 1 )//--IN EITHER CASE REMOVING THIS 'N' ELIMINATES ONE RM 
	n--;
if( remod_cand(pL,x)  && pos[pL].state ==1 ) //--POTENTIAL REMODELLING SITE LOST.
	n--;
if( remod_cand(pL,pR) && pos[pL].state ==1 && pos[pR].state ==1) //--BUT IT COULD CREATE A NEW ONE
	n++;
//---------------------------------------------------------------------


//--------------CHANGE STATE AND INCREMENT OCCUPATION TIME---------

if(pos[x].state != 1)	
	{
	*log << "\n time = " << t << endl;
	cout << "\n time = " << t << endl;

	flag=103;
	process_error( x );
	}
pos[x].state = 0;

if( (pos[x].t_bind_N <0) || (pos[x].t_bind_N > t) )
	{ 
	flag = 4001; 
	process_error(x);
	}
else
	{
	pos[x].t_occ_N += (t - pos[x].t_bind_N);
	pos[x].t_bind_N = -1;
	}

//----------------------------------------------

if( Nucnum > 1 )
    {
    i= right(x);
    while(1)
	{
	pos[i].part_left=pL;	// make this position's next left particle 
				// the former left particle of the current one.
	if( (pos[i].state == 1) || ( pos[i].state == 2) )
		break;

	i=right(i);	// increment position under consideration.
	}
    max=i;

    i=left(x);
    while(1)
	{
	pos[i].part_right=pR;	// make this position's next left particle the former left particle of the current one.
	if( (pos[i].state == 1) || (pos[i].state == 2) )
		break;	

	i=left(i);
	}
    min=i;

	//-----------Optimization (causes wrap-around problems when NNRANGE is larger than 1/2 Llim)------------------
	if ( distance(x,max) > (NNRANGE))
		max = ((x+(NNRANGE))%Llim);
	if ( distance(min,x) > (NNRANGE))
		min = x-(NNRANGE);
	if(min<0)
		min=min+Llim;
	//----------------------------------------
    }
else if( Nucnum == 1 )
	{
	min = x;
	max = x;
	
	reset_prevnext_pointers();
	}
else
	{		//----  pairs that are eligible for remodelling.
	flag=235623; //----UNCLEAR NUMBER OF PARTICLES.
	recentRx=-20;
	process_error(x);	
	exit(1);
	}

Nucnum--;
partnum--;

recentRx=0;


return calc_rates( min, max,n); 
}
//*****************************************************************

int GAdata::add_Nuc( int x)
{
int i, min=0, max=0;
bindevent Q; Q.t = t;  Q.species = 1;  Q.onoroff = 1;

int pR = pos[x].part_right;
int pL = pos[x].part_left;

if(pos[x].state != 0)
	{
	flag=104;
	process_error( x );
	}

//------------------REMODELLING SITE INCREMENT-------------------------
int n=0;

if( remod_cand(x,pR)  && pos[pR].state ==1 )//--IN EITHER CASE ADDING THIS 'N' CREATES ONE
	n++;
if( remod_cand(pL,x)  && pos[pL].state ==1)//--POTENTIAL REMODELLING SITE.
	n++;
if( remod_cand(pL,pR) && pos[pL].state ==1 && pos[pR].state ==1)//--BUT IT COULD DESTROY A PREVIOUS ONE
	n--;
//---------------------------------------------------------------------


//----------UPDATE STATE AND START OCCUPATION TIMER-------
pos[x].state=1;

if( pos[x].t_bind_N > 0  )
	{ 
	flag = 4002; 
	process_error(x);
	}
else
	{
	pos[x].t_bind_N = t;
	}

//------------------------------------------------------------------------

i=right(x);
while(1)
	{
	pos[i].part_left=x;	
	// i's next left particle to the left is now x.

	if( ( pos[i].state == 1) || ( pos[i].state == 2) )
		{
		break;
		}
	else
		{
		i=right(i);	// increment position under consideration.
		}
	}
max=i;

i=left(x);
while(1)
	{
	pos[i].part_right=x;	// make this position's next left particle the former left particle of the current one.
	if( ( pos[i].state == 1) || ( pos[i].state == 2) )
		{
		break;	
		}
	else
		{
		i=left(i);
		}
	}
//------------------------------------------------------------------------
min=i;

Nucnum++;
partnum++;
recentRx=1;

return calc_rates(min, max,n); 
}
//*****************************************************************
int GAdata::slide_Nuc_left( int x)
{
int xp = left(x);	//---xp denotes where it's going to:
int pR = pos[x].part_right;
int pL = pos[x].part_left;

int i, min, max;

bindevent Q;  Q.t = t;  Q.species = 1;  Q.onoroff = -1;

//------------------REMODELLING SITE INCREMENT MOVING x LEFT TO xp-------------------------
//-----bookkeeping for whether the number of eligible remodeller pairs has changed through this reaction (on either side of x).
int n=0;

if( pos[pR].state ==1)
	n += ( remod_cand(xp,pR)  - remod_cand(x,pR) ) ; 

if(  pos[pL].state ==1)
	n += (  remod_cand(pL,xp) - remod_cand(pL,x) ) ;



//----------UPDATE STATE AND START OCCUPATION TIMER-------
if( ( pos[xp].state != 0) || (pos[x].state != 1) )	// error flag
	{
	flag=105;
	process_error( x );
	}
pos[x].state=0;
pos[xp].state=1;

if( (pos[xp].t_bind_N > 0 ) || ( pos[x].t_bind_N < 0  ) || (pos[x].t_bind_N > t)  )
	{ 
	flag = 4003; 
	process_error(x);
	}
else
	{
	pos[xp].t_bind_N = t;
	
	pos[x].t_occ_N  += (t - pos[x].t_bind_N );
	pos[x].t_bind_N = -1;
	}


//------------------------------------------------------------------------
i=x;
while(1)
	{
	pos[i].part_left=xp;	// make this position's next left particle the former left particle of the current one.
	if( (pos[i].state == 1) || (pos[i].state == 2) )
		break;
			// that there is only one particle in the whole system.
	i=right(i);	// increment position under consideration.
	}

max=i;

i=left(xp);
while(1)
	{
	pos[i].part_right=xp;	// make this position's next left particle the former left particle of the current one.
	if( (pos[i].state == 1) || (pos[i].state == 2) )
		break;	

	i=left(i);
	}

if(partnum ==1)
	pos[xp].part_right=xp; 
else	
	pos[xp].part_right=pR; 	//whatever was to the right of the old x

//------------------------------------------------------------------------
min=i;

recentRx=2;

return calc_rates( min, max,n); 
}

//*****************************************************************
int GAdata::slide_Nuc_right( int x)
{
int xp = right(x);	//---xp denotes where it's going to:
int pL = pos[x].part_left;
int pR = pos[x].part_right;
int i, min, max;

bindevent Q;  Q.t = t;  Q.species = 1;  Q.onoroff = -1;

//------------------REMODELLING SITE INCREMENT-------------------------
//-----bookkeeping for whether the number of eligible remodeller pairs has changed through this reaction (on either side of x).

int n=0;

if( pos[pR].state ==1)
	n += (  remod_cand(xp,pR)  - remod_cand(x,pR) ) ;
if( pos[pL].state ==1 )
	n += (  remod_cand(pL,xp)  - remod_cand(pL,x) );


//----------UPDATE STATE AND START OCCUPATION TIMER-------
if( ( pos[xp].state != 0) || ( pos[x].state != 1))	// error flag
	{
	flag=106;
	process_error( x );
	}
pos[x].state=0;
pos[xp].state=1;
	
if( (pos[xp].t_bind_N > 0 ) || ( pos[x].t_bind_N < 0  ) || (pos[x].t_bind_N > t)  )
	{ 
	flag = 4004; 
	process_error(x);
	}
else
	{
	pos[xp].t_bind_N = t;
	
	pos[x].t_occ_N  += (t - pos[x].t_bind_N );
	pos[x].t_bind_N = -1;
	}

//------------------------------------------------------------------------
i=right(xp);
while(1)
	{
	pos[i].part_left=xp;	// make this position's next left particle the former left particle of the current one.
	if( (pos[i].state == 1) || (pos[i].state == 2) )
		break;	

	i=right(i);
	}

max=i;
//------------------------------------------------------------------------
i=x;
while(1)
	{
	pos[i].part_right=xp;	// make this position's next left particle the former left particle of the current one.
	if( ( pos[i].state == 1) || ( pos[i].state == 2) )
		break;

	i=left(i);	// increment position under consideration.
	}
if(partnum ==1)
	pos[xp].part_left=xp; 
else	
	pos[xp].part_left=pL; 	//whatever was to the left of the old x
min=i;

recentRx=3;


return calc_rates( min, max,n); 
}
//*****************************************************************
int GAdata::calc_rates(  const int minpos, const int maxpos, const int n) //---n is the number of new possible remodelling site pairs
{
long double da0 = 0.0;

int temp;

int dright, dleft;

double new_a_addN;		//delta-Energy of removal from system
double new_a_removeN;
double new_a_sNL;		//delta-Energy of sliding left one position
double new_a_sNR;		// "	"	"	" right "	"

int i,x;
bool crossed_middle=false;
bool passed_min_once; 

//-------------------------INITIALIZE SCAN PARAMETERS-------------------------
x=minpos;
passed_min_once = false;

int neighbours[4];

double Boltzdiffuse; //----the boltzmann factor used only as a temp variable for diffusion.

// neighbours[1]= location left
// neighbours[3]= location right

neighbours[1] = pos[x].part_left;  //locL -location of particle to the left
neighbours[0] = pos[neighbours[1]].state; //typeL -type of particle to the left
neighbours[3] = pos[x].part_right; //--locR
neighbours[2] = pos[neighbours[3]].state; // typeR should always be 1... ok, as long as deADD energy func always gets this value
//----  DON'T CHANGE THE ABOVE ORDER   ------------------------------------------------

if ( (partnum > 0) && ( neighbours[0] != 1   ||  neighbours[2] != 1 ))
	{
	flag=114;
	process_error( x );
	}
//--------------------------START SCAN----------------------------------------
double HNGrates[reactions_per_site];
for(i=0;i<reactions_per_site;i++)
	HNGrates[i] =0.0;
//--- 'x' is initialized to minpos above.

passed_min_once = false;
while(1)
{

//------FOR OPTIMIZATION: KEEP TRACK OF WHAT CHANGES HAVE ALREADY TAKEN PLACE 
//------SO AS NOT TO HAVE TO ADD UP ALL POSSIBLE NEW TRANSITIONS.
switch(pos[x].state )
	{
	
	case 0:	//---unoccupied lattice site.------------------------------------------------------------------------------------------------------
	  if ( HNG )
		{
		get_HNG_rates(0,HNGrates, neighbours, x);
		new_a_removeN  = HNGrates[0]; new_a_addN  = HNGrates[1]; new_a_sNL  = HNGrates[2]; new_a_sNR  = HNGrates[3];
		}
	  else
		{ //--- SNG or LNG------------------------------------------- :


		if( boltzmann_on_add || boltzmann_on_uphill )
			{
			new_a_addN  = ka_N*k_E( -muN(x)  + interaction_dEadd(1, x, neighbours ) );//--neighbours will always 
												  // get "1" for nucleosomes.
			}
		else if( boltzmann_on_removal )
			{
			new_a_addN  = ka_N*k_E( -muN(x) ); ///----mu is ALWAYS weighted on addition.
			}
		else if( boltzmann_on_addrem_intermed )
			{
			new_a_addN  = ka_N*k_E( -muN(x) + BZalpha * interaction_dEadd(1, x, neighbours )  ); 
			//---here the addition reaction takes the full mu, but only takes account of BZalpha*the interaction energy.
			}
			//--k_E takes as an argument the energy change that will result.
			//--interaction already accounts for the offset of the outside interactions

		new_a_removeN  = 0.0;
		new_a_sNL  = 0.0;
		new_a_sNR  = 0.0;

		}

	break;
	
	case 1: //---Nucleosome here -----------------------------------------------------------

		neighbours[3] = pos[x].part_right; //--location of the next particle to the right
		// neighbours[2] just stays 1 as always 

	   if(HNG)
		{		
		get_HNG_rates(1,HNGrates, neighbours, x);
		new_a_removeN  = HNGrates[0]; new_a_addN  = HNGrates[1]; new_a_sNL  = HNGrates[2]; new_a_sNR  = HNGrates[3];
		}
	   else
		{ //--- SNG or LNG----------:	

		if (boltzmann_on_add)
			{
			new_a_removeN  = ka_N*1.0;
			}

		else if ( boltzmann_on_removal || boltzmann_on_uphill )
			{
			new_a_removeN  = ka_N*k_E( - interaction_dEadd(1, x, neighbours) );
			}
		else if ( boltzmann_on_addrem_intermed )
			{
			new_a_removeN  = ka_N*k_E( (1.0-BZalpha)*(-interaction_dEadd(1, x, neighbours)) );
			}
		else
			{
			flag=41267; cout <<"\n ERROR 1 : not sure when to apply boltzmann factor.\n";
			process_error( x );
			}

		if(pos[left(x)].state ==0 && ks_N > 0.0 )
			{
			new_a_sNL = ks_N*gsl_sf_exp( -1.0 * max ( dEsL(1, x, neighbours ),0  )  ); 
			//---unlike addition/removal, diffusion is ALWAYS weighted by the boltzmann 
			//---factor at UPHILL steps (i.e. punished.)
			}
		else
			{
			new_a_sNL = 0.0;
			}

		if(pos[right(x)].state ==0 && ks_N > 0.0 )
			{
			new_a_sNR = ks_N * gsl_sf_exp( -1.0 * max ( dEsR(1, x, neighbours ),0  )  );
			//---unlike addition/removal, diffusion is ALWAYS weighted by the boltzmann 
			//---factor at UPHILL steps (i.e. punished.)
			}
		else
			{
			new_a_sNR = 0.0;
			}

		new_a_addN  = 0.0;
		
		}

	   if(bind_irrev || pos[x].permanent) //--- do this after previous calculation and overwrite the  
		{	  //--- previously determined rate if the binding is irreversible.
		 
		new_a_removeN  = 0.0; 
		}

	   neighbours[1] = x;  // locL
	   // ---neighbours[0] stays at 1; // typeL	//----(from here on, shift the next neighbour to being THIS one)

	break;


	default: //----here is somewhere inside a TF -all reactions are zero here.

		cout << "\n ERROR: encountering a state neither 1 nor 0. Exiting." << endl;
		exit(1);
	break;

	}  //----end the switch-case structure to determine what the state is 

	//-------  KEEP TRACK OF HOW MUCH a0 IS CHANGING------
	 da0 += (new_a_addN -    *pos[x].a_addN) ;
	 da0 += (new_a_removeN - *pos[x].a_removeN) ;
	 da0 += (new_a_sNL -     *pos[x].a_sNL) ;
	 da0 += (new_a_sNR -     *pos[x].a_sNR) ;

	//-------  NOW UPDATE THE a VALUES   ------------------

	*pos[x].a_addN     = new_a_addN;
	*pos[x].a_removeN  = new_a_removeN ;
	*pos[x].a_sNL      = new_a_sNL;
	*pos[x].a_sNR      = new_a_sNR;

	//----NOW DO IT AGAIN FOR THE TF REACTIONS, BUT ONLY TFs ARE PRESENT to save efficiency
	
	if ( (x==maxpos && minpos != maxpos) || (x==maxpos && passed_min_once) ) /// @@@ 
		{
		break;
		}
	else
		{
		if(x == minpos)
			{
			passed_min_once=true;
			}
		x=right(x);
		}
   }// ---- end the scan running through the xmin->xmax segment.

//------------------------------------------------------------------------------------
if( krm_b  ) //--if we have active remodelling
  {
  N_remod += n;
  da0     += 2.0*n*krm_val;
  }
//------------------------------------------------------------------------------------


a0+=da0;

return 1;
}//----end the calc_rates 
//*****************************************************************
int GAdata::get_HNG_rates(  const int type, double * HNG_outputrates, const int * neighbours, const int x ) //----subroutine for calculating the HNG transition rates.
{

double new_a_removeN;
double new_a_addN;		//delta-Energy of removal from system
double new_a_sNL;		//delta-Energy of sliding left one position
double new_a_sNR;		// "	"	"	" right "	"


switch( type)
	{
	case 0:
		//---------------------|<------------------> <------>
		//---------------------|==============|     x       |==============|
	   if (  (distance(pos[x].part_left,x) < footprint &&  neighbours[0] ==1 ) || ( (distance(x,pos[x].part_right) < footprint) )   )
		{
		new_a_addN  = 0.0;
		}
	   else
		{ 
		new_a_addN  = ka_N*k_E( -muN(x) ); //---always weight mu on addition.
		}

	new_a_removeN  = 0.0;

	new_a_sNL  = 0.0;
	new_a_sNR  = 0.0;

	break;

	case 1: //------------------nucleosome --------------------

	new_a_removeN  = ka_N*1.0;

	if ( distance(pos[x].part_left,x) < footprint || distance(x,pos[x].part_right) < footprint )
		{
		flag=131429;
		process_error( x );
		}


	if( (distance(pos[x].part_left,x) <= footprint &&  neighbours[0] ==1 ) ||  pos[left(x)].state !=0 )
		{
		new_a_sNL = 0.0;
		}
	else
		{
		new_a_sNL = ks_N*gsl_sf_exp( -1.0 * max ( muN(left(x))- muN(x) , 0  )  );
		}

	if(  distance(x,pos[x].part_right) <= footprint)    
		{
		new_a_sNR = 0.0;
		}
	else
		{
		new_a_sNR = ks_N*gsl_sf_exp( -1.0 * max ( muN(right(x))- muN(x) , 0  )  );
		}

	new_a_addN  = 0.0;

	break;
	

	default://---------------------ERROR CATCH ------------------------------
		cout <<"\n HNG analysing senseless case!!! \n ";
		*log <<"\n HNG analysing senseless case!!! \n ";
		(*log).close();
		exit(1);
	}//---END OF "case" STRUCTURE

/* order is:
0 removeN 
1 addN    
2 sNL     
3 sNR     

*/

//----------output is always the same---------------
HNG_outputrates[0]=new_a_removeN;
HNG_outputrates[1]=new_a_addN;
HNG_outputrates[2]=new_a_sNL;
HNG_outputrates[3]=new_a_sNR;

return type;
}



//******************************************************************************
//------     HERE ARE THE OBSERVABLES  AND ERROR CHECKS  -----------
//******************************************************************************

bool GAdata::should_observe_filling( const double tau, int ** void_histogram, config_set_t * Z_all_t ) //--- for the purposes of determining the filling rates of the system
{
bool result=false;
int repetition =0; 
int j;

// for(j=0;j<total_obs_filling;j++)
 for(j=obs_count_filling;j<total_obs_filling;j++) //---this should be the same, but faster.
	{
	if(  ( tpoints_filling[j] > t)  &&  ( tpoints_filling[j] <= (t+tau))  ) // if this time-step (t->t+tau) straddles the observation point tpoints_filling[j] (can happen for more than one j).
		{ 
		if(result) //-----  have we already observed at this time step?
			{  //-----  if so, then this time step straddels more than one time step.
			   //-----  THIS CAN HAPPEN A FEW TIMES AND ITS OK - it just means we are
			   // coarse-graining.
			repetition++;
//			*log << "\n WARNING -repeating observation " << repetition << " times in should_observe_filling.\n";
			}
		else
			{
			result=true;
			}

		get_filling_frac( );
		increment_void_histogram(void_histogram);
		
		if (calculate_entropy)
			{
	 		grab_current_configuration( Z_all_t[j] );
			Z_all_t[j].tpoint_passed = true;
			}

		if( j != obs_count_filling )
			{//------    this should really NEVER happen-------
			*log << "\n error! faulty observation count number.\n";
			cout << "\n error! faulty observation count number.\n";
			exit(1);
			}
		obs_count_filling++;


		}  

	if( tpoints_filling[j] >  (t+tau) )
		{
		break; //@@@ also add an initializer of current time plot to start j there.
		}

	}

return result;
}
//******************************************************************************
int GAdata::should_observe_patterns( const double tau ) 
{//--- tests whether it is the right time to update the two-body corr. and the one-point histograms

int result=0;
int repetition =0; 
int j;

//---------   FIRST TEST FOR WHETHER WE ARE AT ONE OF THE EQUILIBRIUM POINTS:  ----------------
for(j=0;j<total_obs_eq;j++)
	{
	if(  ( tpoints_eq[j] > t)  &&  ( tpoints_eq[j] <= (t+tau))  ) 
		{ 
		//---- If this time-step straddles an observation point.

		increment_2_part_corr( -1 ); 
		increment_1_point_hist( -1 ); //--THE '-1' MEANS EQUILIBRIUM.
		
		if( j != obs_count_eq_2pc)
			{//------    this should really NEVER happen-------
			*log << "\nerror! faulty observation count number.\n";
			exit(1);
			}

		obs_count_eq_2pc ++;
		}  

	}

//---------    NOW TEST FOR WHETHER WE ARE AT ONE OF THE TRANSIENT POINTS:   -----------------
for(j=0;j<total_obs_filling;j++)
	{
	if(  ( tpoints_filling[j] > t)  &&  ( tpoints_filling[j] <= (t+tau))  ) 
		{
		//--- If this time-step straddles an observation point.

		increment_2_part_corr( j );
		increment_1_point_hist( j ); //--THE -1 MEANS EQUILIBRIUM.		

		if( j != obs_count_2pc_ti )
			{//------    this should really NEVER happen-------
			*log << "\nerror! faulty observation count number.\n";
			exit(1);
			}

		obs_count_2pc_ti++;
		}  

	}
//------------------------------------------------------------------------------------------


return result;
}

//******************************************************************************
int GAdata::should_observe_equilibrium_voiddist( const double tau, int * void_histogram_equilibrium  ) 
{ //--- tests whether it is the right time to update the two-body corr.
int result=0;
int j;

for(j=0;j<total_obs_eq;j++)//---total_obs_eq is the number of observations we are making at equilibrium time points.
	{
	if(  ( tpoints_eq[j] > t)  &&  tpoints_eq[j] <= (t+tau) ) 
		{ // if this time-step straddles an observation point, and there are no TFs present.

		increment_void_histogram_equilibrium( void_histogram_equilibrium  );
		
		
		if( j != obs_count_eq_vdist + result)
			{//------    this should really NEVER happen-------
			cout << "\nerror! faulty observation count number.\n";
			*log << "\nerror! faulty observation count number.\n";

			(*log).close();

			exit(1);
			}

		result++;
		}  

	}

return result;
}
//******************************************************************************
int GAdata::grab_current_configuration( config_set_t & C_t)
{
int i;
int Size_current_tvec = C_t.Z_t.size();	//---the number of unique configurations already collected.
int x;
int result;	//---returns the index of the matched configuration, or (for a new one) -1.

int offset=0;   //---the offset from x=0 to the first particle observed. with this, the following string of separations is a complete configuration. 

vector<int>  dists;

//-----initialize a new one:-------------

configuration  Q_temp; 
// Q_temp.N = Nucnum;
Q_temp.pcount = 1;
Q_temp.E=-Nucnum*muN0;

//---------------------------------------
string gapstr;

if(Nucnum == 0)
	{
	// do nothing here: there should be no voids to consider.
	//------	Q_temp.dists.push_back(Llim);

	Q_temp.description="0";	//---simply "empty"

	}
else 
	{

	offset=pos[Llim-1].part_right;	
	x=offset;


	Q_temp.description= bren_itoa(Nucnum)  + string("_") + bren_itoa ( offset ) + string("_") ;

	for(i=0;i<Nucnum;i++)	
		{
		if(LNG || SNG)//---@@@
			{
			Q_temp.E += VNN[distance(x,pos[x].part_right)-1];	
			}
		dists.push_back( distance(x,pos[x].part_right) );
		x=pos[x].part_right;	

		}

	//----after this, 'x' should end up back where it started, at the first particle position.
	if( x != offset )	//---check that after running the loop it comes back to offset .
		{
		flag=1351244;
		cout << "\n ERROR in grab_current_configuration, not aligned back with original x.";
		*log << "\n ERROR in grab_current_configuration, not aligned back with original x.";
		exit(1);
		}
	//----------------done check--------------------
	}

//---------------SET UP THE STRING-----------------------


for (i=0; i< Nucnum ; i++) //---N.B. if N==0, then this will just automatically return 'true', as it should
	{	  

	if ( i ==  Nucnum -1)
		{
		gapstr = bren_itoa ( dists.at(i) ) ;
		}
	else
		{
		gapstr = bren_itoa ( dists.at(i) ) + string("-");
		}

	Q_temp.description = Q_temp.description  + gapstr ;

	}


//------we now have our full current configuration. Now look for a match in the current Z.
bool match_found=false;

i=0;
while(i < Size_current_tvec && !match_found)
	{
	if( Q_temp.description ==   C_t.Z_t.at(i).description )	
		{	// if the current configuration matches this 
		match_found = true;
		C_t.Z_t.at(i).pcount+=1;	//--- increment the probability 
						//--- of this configuration by one. 
		result=i;
		}		
	i++;
	}


int d=0; //---distance between adjacent nucl's


if( Size_current_tvec==0 || !match_found)	//--- DIDN'T FIND A MATCH IN THE EXISTING PART. FUNC. 
	{					//--- ADD A NEW ENTRY
	
	C_t.Z_t.resize(Size_current_tvec+1); //--- extend the size by one in order to add one new 
					     //--- configuration to the partition func for this time.
/*
	C_t.Z_t[Size_current_tvec].pcount = 1;	//----- set the pcount of the new one that 
						//----- just got added to p1
	C_t.Z_t[Size_current_tvec].description = Q_temp.description;	
*/

	C_t.Z_t[Size_current_tvec] = Q_temp; //----- just copy the whole thing at once.

	result = -1;	//----the '-1' signifies that this was a new configuration.
	}

//------------------------ now add the H contribution together -------------------
x=pos[Llim-1].part_right;

for(i=0;i<Nucnum;i++)	
	{
	C_t.U += interaction_NN( distance(x,pos[x].part_right)) ;
	x=pos[x].part_right;		
	}

//--------------- NOW INCREMENT Nave----------------------------------------------
	C_t.Nave += Nucnum;
//--------------------------------------------------------------------------------
return result;
}

//******************************************************************************

int get_first_N_pvals_from_config_list( const vector < configuration >   Z_t, float * target, const int N)
{ //----takes the configuration list Z_t, and orders the vector entries by p-value to output into the target array.

int i,j;
int Size_current_tvec = Z_t.size();	//---the number of unique configurations already collected.
int x;
int result;	//---returns the index of the matched configuration, or (for a new one) -1.

int temp[Size_current_tvec];
double temp_ordered[Size_current_tvec];

int current_val;
int dummy=0;


//----------initialize arrays----------------
for (i=0;i<Size_current_tvec;i++)
	{
	 temp[i]         = Z_t.at(i).pcount;
 	 temp_ordered[i] = 0;

	 if (i<N)//----only initialize the target array within its size.
		{
		target[i] = 0.0;
		}
	}


//----------- put in the ordered list ---------

for (i=0; i<Size_current_tvec; i++)
	{
	current_val=temp[i];
	//------------------------------------------
	
	for (j=0; j<i; j++)
		{
		 if( current_val > temp_ordered[j] )
			{
			dummy             = temp_ordered[j];
			temp_ordered[j]   = current_val;
			current_val = dummy; 
			}
		 //----else do nothing, just keep current_val the same.
		}

	temp_ordered[j] = current_val;
	}

//--------TAKE THE FIRST N VALUES AND PUT THEM INTO TARGET
for(i=0;i<N;i++)
	{
	if (i >= Size_current_tvec)
		{
		break;  //---  catch, in case we have less than N distinct states, 
			//--- so we don't read off end of array.
		}
	target[i]=temp_ordered[i];
	}

//---- DONE ----------------------------------------------------------------------

return 1;
}

//******************************************************************************

bool compare_dists ( const vector<int> A, const vector<int> B , const int N)
{
int i;
if( A.size() != N || B.size() != N )
	{
	cout << "\n ERROR in compare_dists: sizes don't match.\n";
	exit(1);
	}

bool result=true;

for (i=0;i<N;i++) //---N.B. if N==0, then this will just automatically return 'true', as it should
	{	  //---zero particles means zero voids, and any two such configurations must be the same.
	if( A.at(i) != B.at(i) )
		{
		result = false;
		break;
		}
	}
return result;	
}
//******************************************************************************

int GAdata::increment_2_part_corr( int timepoint )
{//---calculates the distributed 2-particle correlation function for the entire array.

 //--- timepoint specifies the index of the temporal positions at 
 //--- which the current distribution should be taken. '-1' implies the equilibrium value.

 //--- we then assign the local temporary pointer to the address of the corresponding 
 //--- distribution and proceed as usual.


//-----------------------------------------------------------
double * two_part_corr;	//----temporary array that can point to either of the sets of 2pc arrays we want.
if( timepoint == -1)
	{
	two_part_corr = two_part_corr_eq;
	}
else
	{
	two_part_corr = two_part_corr_ti[timepoint];
	}
//-----------------------------------------------------------


int N_so_far =0;
int ref_point; 
	
int x=0;
int i=0;

int	first_part = pos[Llim-1].part_right;
int	last_part = pos[0].part_left;

	
if(Nucnum < 1 )
	{
	goto 	exiting_incrementing_2partcorr;
	}


while( pos[first_part].state ==2 )
	{
	first_part=pos[first_part].part_right;
	}
	
ref_point = first_part;
	
while( ! (N_so_far >Nucnum) )
	{
	if(pos[ref_point].state != 1 )
		{
		flag=8193546;
		process_error(ref_point);
		exit(1);
		}
		
	x=ref_point;
		
	for(i=0;i<Llim;i++) //'i' refers to nothing, we just want to make sure to do this Llim times.
		{
		if(pos[x].state == 1 )
			{
			two_part_corr[i] += 1.0;
			}
		x=right(x);
		}
		
	//-----	CHECK FOR WHERE WE'VE ENDED UP ------	
	if(x!=ref_point)
		{
		flag = 54729;
		process_error(x);
		exit(1);
		}
	//---------------------------------------------
		
		
	//-------   NOW INCREMENT THE NUMBER OF PARTICLES CONSIDERED AND MOVE ON TO THE NEXT ONE   -------
	N_so_far ++;
		
	if( ref_point == last_part)
		break;
	
	else
		{
		ref_point = pos[ref_point].part_right;
		
		while( pos[ref_point].state ==2 )
			{
			if( ref_point == last_part)
				{
				goto 	exiting_incrementing_2partcorr;
				}				
			else
				ref_point=pos[ref_point].part_right;
			}

		}
		
	}
	
	exiting_incrementing_2partcorr:
	
	if (N_so_far != Nucnum)
	{
	flag = 1004562;
	process_error(1);
	exit(1);	
	}
	
//-----------------------------------------------------------
if( timepoint == -1)
	{
	num_times_2pc_eq_incremented += Nucnum ;	// if it's the equilibrium value, 
							// just one array counted.
	}
else
	{
	num_times_2pc_ti_incremented[timepoint] += Nucnum ;	// if it's one of the transient time
								// points, then we need to increment
	}							// just that one.
//-----------------------------------------------------------


	
return 1;

}//----end the calc_rates 

//*****************************************************************
int GAdata::increment_1_point_hist( int timepoint ) 
{//--- increments the fixed 1-point histogram function for the entire array.

 //--- timepoint specifies the index of the temporal positions at 
 //--- which the current distribution should be taken. '-1' implies the equilibrium value.

 //--- we then assign the local temporary pointer to the address of the corresponding 
 //--- distribution and proceed as usual.


//-----------------------------------------------------------
int * one_point_hist_N;	//----temporary array that can point to either of the sets of 2pc arrays we want.
if( timepoint == -1)
	{
	one_point_hist_N  = onepoint_histocc_N_eq;
	}
else
	{
	one_point_hist_N  = onepoint_histocc_N_ti[timepoint];
	}

//-----------------------------------------------------------

int x;

for(x=0;x<Llim;x++)
	{
	if( pos[x].state == 1 )	//---if there is a nuc here at this time.
		{
		one_point_hist_N[x]++;
		}
	}

	
return 1;

}//----end the calc_rates 


//*****************************************************************

int GAdata::normalize_2_part_corr(void ) //---n is the number of new possible remodelling site pairs
{
// --------------WE DON'T ACTUALLY USE THIS IN THE GCE
if( num_times_2pc_eq_incremented == 0  )
	{
	cout << "\n ERROR: taking average, but never incremented.\n";
	*log << "\n ERROR: taking average, but never incremented.\n";
	exit(1);
	}	

int i,j;
//-------------------   FIRST THE EQUILIBRIUM CASE :   ---------------------------
for(i=0;i<Llim;i++)
	{
	two_part_corr_eq[i] = 	two_part_corr_eq[i] / float(num_times_2pc_eq_incremented);
	}

//-------------------   NOW THE TRANSIENT CASE :    -----------------------------
for(i=0; i<total_obs_filling; i++)
	{
	if( num_times_2pc_ti_incremented[i] > 0)
		{
		for(j=0;j<Llim;j++)
			{
			two_part_corr_ti[i][j] = two_part_corr_ti[i][j] / float(num_times_2pc_ti_incremented[i]);
			}
		}
	}
//----------------------------------------------------------------------------

//---- this now sets the two-point correlation functions to their num-count averages, and 
//---- we increment this onto the stack from all runs in main. runs with greater nuc numbers
//---- are not given more weight in this average.

}

//*****************************************************************
int GAdata::get_filling_frac(  ) //---n is the number of new possible remodelling site pairs
{
double rhoi=0.0;

rhoi= (double(Nucnum))/(double(Llim)*CGF);   // we just want the number of particles per unit length.

/*!!!
if(fabs( t - obs_count_filling*dt_obs) > dt_obs )
	{
	*log << "\n WARNING in filling_frac -dt_obs large enough to be on scale with individual reactions.\n";
	}
!!!*/
filling_frac[obs_count_filling] = rhoi;
}
//****************************************************************************
int GAdata::increment_void_histogram(int** void_hist) //---- this just counts the incidence of voids, doesn't normalize for time or voidnum.
{
int N_so_far =0;
int ref_point; 
	
int x=0;
int i=0;

int	first_part = pos[Llim-1].part_right;
int	last_part = pos[0].part_left;

int p2p; //---interparticle distance (from the same landmark on the particle).
int void_size;	//-self-explanatory.

//-----------------------------------------------------------------------
if(Nucnum < 1 )
	{
	void_hist[obs_count_filling][Llim] ++;
	goto 	exiting_increment_void_histogram;
	}
else
   {
   ref_point = first_part;
   p2p = distance(ref_point,pos[ref_point].part_right );

   if(HNG)
     {void_size=p2p-footprint;}
   else if(SNG||LNG)
     {void_size=p2p-1;}

   if(void_size<0)
	{
	flag=926740;
	process_error(void_size);
	exit(1);
	}
   void_hist[obs_count_filling][void_size] ++; 
   N_so_far++;

   while( pos[ref_point].part_right != first_part )
	{

	ref_point = pos[ref_point].part_right;
 	p2p       = distance(ref_point,pos[ref_point].part_right );

	if(HNG)
	    {void_size=p2p-footprint;}
	else if(SNG||LNG)
	    {void_size=p2p-1;}
		
	void_hist[obs_count_filling][void_size] ++; 
	N_so_far ++;						
	}

  }//----end the 'else' from if there are no particles in the system.	

	exiting_increment_void_histogram:
	
	if (N_so_far != Nucnum)
	{
	flag = 1213462;
	process_error(Nucnum);
	exit(1);	
	}
	

return 1;
}
//****************************************************************************
int GAdata::increment_void_histogram_equilibrium(int* void_hist)
{
int N_so_far =0;
int ref_point; 
	
int x=0;
int i=0;

int	first_part = pos[Llim-1].part_right;
int	last_part = pos[0].part_left;

int p2p; //---interparticle distance (from the same landmark on the particle).
int void_size;	//-self-explanatory.

//-----------------------------------------------------------------------
if(Nucnum < 1 )
	{
	void_hist[Llim] ++;
	goto 	exiting_increment_void_histogram_equilibrium;
	}
else
   {
   ref_point = first_part;
   p2p = distance(ref_point,pos[ref_point].part_right );

   if(HNG)
     {void_size=p2p-footprint;}
   else if(SNG||LNG)
     {void_size=p2p-1;}

   if(void_size<0)
	{
	flag=926739;
	process_error(void_size);
	exit(1);
	}
   void_hist[void_size] ++; 
   N_so_far++;

   while( pos[ref_point].part_right != first_part )
	{

	ref_point = pos[ref_point].part_right;
 	p2p       = distance(ref_point,pos[ref_point].part_right );

	if(HNG)
	    {void_size=p2p-footprint;}
	else if(SNG||LNG)
	    {void_size=p2p-1;}
		
	void_hist[void_size] ++; 
	N_so_far ++;						
	}

  }//----end the 'else' from if there are no particles in the system.	

	exiting_increment_void_histogram_equilibrium:
	
	if (N_so_far != Nucnum)
	{
	flag = 1213462;
	process_error(Nucnum);
	exit(1);	
	}
	

return 1;
}

//************THIS NEXT FUNCTION TAKES THE STATISTICS FROM THE HISTOGRAM******
int calc_void_statistics( double * void_means, double * void_stddevs, double * rhocheck,int const* const * const void_histogram, const int total_obs_filling, const int Llim, const double CGF , const int numtrials)
{
int tp=0; //---'tp'=timepoint, when this is occuring
int m=0;  //---the void size
int N=0; //---temp variable holding the number of voids at each time-point.
double mean,var;

for(tp=0;tp<total_obs_filling;tp++)	//-----loop over time points.
	{
	N    =0;
	mean =0.0;
	var  =0.0;

	for(m=0;m<Llim;m++)	//-----loop over void sizes -exclude L just for now.
		{
		N+=void_histogram[tp][m];
		}

	rhocheck[tp] = (double(N))/(double(Llim)*CGF*numtrials);//---check the density as a function of time. 
	N            +=void_histogram[tp][Llim]; //---now add on the last one (which doesn't count towards rhocheck)
	//----now we have 'N' -the total number of voids that occured at timepoint tp.

	for(m=0;m<=Llim;m++)
		{
		mean              += m*CGF*(void_histogram[tp][m] / (double(N))) ; //---m*CGF is the void size, void_hist/N is the prob.
		}
	void_means[tp] =mean;

	for(m=0;m<=Llim;m++)
		{
		var  += (  (double(m)*CGF-mean)*(double(m)*CGF-mean) ) * ( void_histogram[tp][m] /double(N) ) ; 
		//---calculate the std.dev.
		}
	void_stddevs[tp] = sqrt( var );
	}

return 1;
}



//****************************************************************************

int GAdata::check_rates( void )
{
long double temp=0.0;
int i = 0;
int x = 0;
int result;

int neighbours[reactions_per_site];

double check_a_addN=0.0;	//delta-Energy of removal from system
double check_a_removeN=0.0;
double check_a_sNL=0.0;		//delta-Energy of sliding left one position
double check_a_sNR=0.0;		// "	"	"	" right "	"



double HNGrates[reactions_per_site];
for(i=0;i<reactions_per_site;i++)
	HNGrates[i] =0.0;

	//--------------NOW SCAN THROUGH THE ARRAY-------------

for(x=0;x<Llim;x++)
   {

   neighbours[1] = pos[x].part_left;  // locL
   neighbours[0] = pos[neighbours[1]].state; // typeL
   neighbours[3] = pos[x].part_right; //--locR
   neighbours[2] = pos[neighbours[3]].state; // typeR

   //------FOR OPTIMIZATION: KEEP TRACK OF WHAT CHANGES HAVE ALREADY TAKEN PLACE 
   //------SO AS NOT TO HAVE TO ADD UP ALL POSSIBLE NEW TRANSITIONS.
   switch(pos[x].state )
	{

	case 0:	//---unoccupied lattice site.

	  if ( HNG )
		{
		get_HNG_rates(0,HNGrates, neighbours, x);
		check_a_removeN  = HNGrates[0]; check_a_addN  = HNGrates[1]; check_a_sNL  = HNGrates[2]; check_a_sNR  = HNGrates[3];
		}
	  else
		{ //--- SNG||LNG------------------------------------------- :

		if( boltzmann_on_add || boltzmann_on_uphill )
			{
			check_a_addN  = ka_N*k_E( -muN(x)  + interaction_dEadd(1, x, neighbours ) );
			}
		else if( boltzmann_on_removal )
			{
			check_a_addN  = ka_N*k_E( -muN(x) );
			}
		else if( boltzmann_on_addrem_intermed )
			{
			check_a_addN  = ka_N*k_E( -muN(x) + BZalpha * interaction_dEadd(1, x, neighbours )  );
			//---here the addition reaction takes the full mu, but only takes account of BZalpha*the interaction energy.
			}


		//--k_E takes as an argument the energy change that will result.
		//--interaction already accounts for the offset of the outside interactions

		check_a_removeN  = 0.0;


		check_a_sNL  = 0.0;
		check_a_sNR  = 0.0;

		}
	break;
	
	case 1: //------  Nucleosome here   ------------------
		

	   if(HNG)
		{		
		get_HNG_rates(1,HNGrates, neighbours, x);
		check_a_removeN  = HNGrates[0]; check_a_addN  = HNGrates[1]; check_a_sNL  = HNGrates[2]; check_a_sNR  = HNGrates[3];
		}
	   else
		{ //--- SNG   or LNG  ----------:	


		if (boltzmann_on_add)
			{
			check_a_removeN  = ka_N*1.0;
			}

		else if ( boltzmann_on_removal || boltzmann_on_uphill )
			{
			check_a_removeN  = ka_N*k_E( - interaction_dEadd(1, x, neighbours) );
			}
		else if ( boltzmann_on_addrem_intermed )
			{
			check_a_removeN  = ka_N*k_E( (1.0-BZalpha)*(-interaction_dEadd(1, x, neighbours)) );
			}
		else
			{
			flag=411357; cout <<"\n ERROR 4: not sure when to apply boltzmann factor.\n";
			process_error( x );
			}

		
		if(pos[left(x)].state == 0 && ks_N > 0.0)
			check_a_sNL = ks_N * k_E( max( dEsL(1, x, neighbours), 0.0) );
		else
			check_a_sNL = 0.0;

		if(pos[right(x)].state == 0  && ks_N > 0.0)
			check_a_sNR = ks_N * k_E( max( dEsR(1, x, neighbours), 0.0) );
		else
			check_a_sNR = 0.0;

		check_a_addN  = 0.0;

		}

	   if(bind_irrev || pos[x].permanent) //---do this at the end, thus overwriting previous assignments intentionally
		{ 
		check_a_removeN  = 0.0; 
		}

	break;


	default: //----here is somewhere inside a TF -all reactions are zero here.
		cout << "\n ERROR: checking outside of range of reactions.\n";
		exit(1);

	break;

	}  //----end the switch-case structure to determine what the state is 

//----CHECK IF ANYTHING IS NEGATIVE----
   if( (check_a_addN < 0.0 ) ||(check_a_removeN < 0.0 ) ||(check_a_sNL < 0.0 ) ||(check_a_sNR < 0.0 )  ) 
	{ 
	flag = -1234;
	process_error(x);
	}
//----EVERYTHING MUST BE POSITIVE ----

//-----a0----------
   if( check_a_removeN == 0.0 ) 
	{
	   if( *pos[x].a_removeN != 0.0 )	   
		{ 
		flag = 7000; 
		process_error(x); 
		}
	}
   else if( gsl_sf_log(fabs(check_a_removeN/(*pos[x].a_removeN))) > 0.00001 ) 
	{ 
	flag = 7101; 
	*log << "\n rate discrepancy in a0 at position x=" << x << ", exiting.\n";
	*pos[x].a_removeN = check_a_removeN;
	process_error(x);
	}
//-----a1----------
   if( check_a_addN == 0.0 ) 
	{
	   if( *pos[x].a_addN != 0.0 )	   
		{ 
		flag = 7200; 
		process_error(x); 
		}
	}
   else if( gsl_sf_log(fabs(check_a_addN/(*pos[x].a_addN))) > 0.00001 ) 
	{ 
	flag = 7102; 
	*log << "\n rate discrepancy in a1 at position x=" << x << ", exiting.\n";
	*pos[x].a_addN = check_a_addN;
	process_error(x);
	}
//-----a2----------
   if( check_a_sNL == 0.0 ) 
	{
	   if( *pos[x].a_sNL != 0.0 )	   
		{ 
		flag = 7300; 
		process_error(x); 
		}
	}
   else if( gsl_sf_log(fabs(check_a_sNL/(*pos[x].a_sNL))) > 0.00001 ) 
	{ 
	flag = 7103; 
	*log << "\n rate discrepancy in a2 at position x=" << x << ", exiting.\n";
	*pos[x].a_sNL = check_a_sNL;
	process_error(x);
	}
//-----a3----------
   if( check_a_sNR == 0.0 ) 
	{
	   if( *pos[x].a_sNR != 0.0 )	   
		{ 
		flag = 7400; 
		process_error(x); 
		}
	}
   else if( gsl_sf_log(fabs(check_a_sNR/(*pos[x].a_sNR))) > 0.00001 ) 
	{ 
	flag = 7104; 
	*log << "\n rate discrepancy in a3 at position x=" << x << ", exiting.\n";
	*pos[x].a_sNR = check_a_sNR;
	process_error(x);
	}
//----------------------

	}//---end for x=0..Llim


//+++++++++++++++++++++++++++++DOWN TO HERE+++++++++++++++++++++++++++++++++


//--------    CHECK THE SUM:  ----------------
temp=0.0;
for(i=0;i<M;i++)
	{
	temp += Rx[i];
	}
if(krm_b)
	{
	temp +=2.0*N_remod*krm_val;
	}

a0=temp;
//-------  DONE CHECKING THE SUM   -------------


return result;
}

