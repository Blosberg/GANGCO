/*----- IMPLEMENTATION OF THE GILLESPIE ALGORITHM ON THE MOTION OF NUCLEOSOMES, TREATED AS 1-D PARTICLES ALONG DNA
// ---last updated on  Fri Oct 17 15:49:24 CEST 2014  by  ga79moz  at location  TUM , xenopus

// ---- written to process void size configuration strings to infer correlations.
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

// -- #include "GA_GC_NTF.h"
#include <bren_lib.h>

double const pi   = 3.14159265358979323846264;
double const kB   = 1.38065E-23;  // in units of J/k
const int charlength =400;	//----the length of the string kept in memory for file-path manipulation.

using namespace std;

int voidsize_corr_get_configs_from_file (string pathin, string prefix, const int time_index, gsl_matrix * M2, gsl_vector * v, const int numruns, const int L);
int voidsize_corr_inc_states_from_string(string config ,int repnum, gsl_matrix * M2, gsl_vector * v1, const int L );

//******************************************************************

int main(int argc, char *argv[])
{
int num_t_points=0, L=0, num_runs=0, ti=0, x=0;
string file_pathin, file_prefix, pathout;

//  ---- TASKID      = atoi(argv[2]); // -int
//  ---- muN_input   = atof(argv[3]); // -double
//  ---- READ IN PARAMETERS FOR THE SIMULATION FROM THE IN FILE.-----------------


char cpath[charlength];
clear_charray(cpath, charlength );
sprintf(cpath, "./GA_voidsize_corr.in"); //----input file should always be in the working directory.
ifstream datin;

datin.open(cpath);
if( datin.fail() )
	{
	cout << "\n ERROR, can't find input file. exiting \n";
	exit(1);
	}

int parity_check;

datin  >> L >> num_runs  >> num_t_points;
datin  >> file_pathin >> file_prefix >> pathout;
datin  >> parity_check;
datin.close();



if (parity_check != 885588)
	{
	cout << "\n ERROR: variables are disordered in input file. Check input file layout.\n";
	exit(1);
	}

double C11[num_t_points]; init_array( C11, num_t_points, 0.0 );
double C12[num_t_points]; init_array( C12, num_t_points, 0.0 );
double C22[num_t_points]; init_array( C22, num_t_points, 0.0 );

double MF_11=0.0,MF_12=0.0,MF_22=0.0;
double rho;
double sf = 1.0/(double(L*num_runs));

for(ti=0; ti<num_t_points; ti++)
	{
	rho=0.0;
	gsl_vector* v1 = gsl_vector_alloc(L+1);
	gsl_vector_set_zero(v1);
	gsl_matrix* M2 = gsl_matrix_alloc(L+1,L+1);
	gsl_matrix_set_zero(M2); 

	voidsize_corr_get_configs_from_file (file_pathin, file_prefix, ti, M2, v1,  num_runs, L);
	for(x=0; x<L; x++)
		{
		rho+=gsl_vector_get(v1,x);
		}
	rho=rho*sf;

	gsl_vector_scale (v1, sf);
	gsl_matrix_scale (M2, sf);

	//------------------------------------------------------------
	
	if( rho > 0.01*sf*sf) //---check if there's at least one.
		{
		MF_11 = gsl_vector_get(v1,2)*gsl_vector_get(v1,2)/rho;
		MF_12 = gsl_vector_get(v1,2)*gsl_vector_get(v1,3)/rho;
		MF_22 = gsl_vector_get(v1,3)*gsl_vector_get(v1,3)/rho;
	
		if(MF_11 > 0.01*sf*sf)
			{
			C11[ti] = (gsl_matrix_get(M2,2,2)- MF_11)/MF_11;
			}
		if(MF_12 > 0.01*sf*sf)
			{
			C12[ti] = (gsl_matrix_get(M2,2,3)- MF_12)/MF_12;
			}
		if(MF_22 > 0.01*sf*sf)
			{
			C22[ti] = (gsl_matrix_get(M2,3,3)- MF_22)/MF_22;
			}

		}
	else
		{
		MF_11 = MF_12 = MF_22 =0.0;
		}

	}
//------------------------ NOW GET THE TIME POINTS -------------------------


clear_charray(cpath, charlength );
sprintf(cpath, "%sthermoquantities_t_S_H_Htot_Nave.txt",file_pathin.c_str());


ofstream fout(pathout.c_str());

for(ti=0; ti<num_t_points; ti++)
	{
	fout << C11[ti]  << "\t" << C12[ti] << "\t" << C22[ti] << endl; 
	}

fout.close();

return 0;
}

