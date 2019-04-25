//stock-pile of various small functions that are useful in a variety of different versions of software
// ---last updated on  Fri Nov 21 17:55:51 CET 2014  by  ga79moz  at location  TUM , xenopus

//  changes from  Fri Nov 21 17:55:51 CET 2014 : testing synchronization across machines

//  changes from  Thu Oct 23 13:44:15 CEST 2014 : added functions for processing void data from the strings of configuration files

//  changes from  Fri Apr 25 11:12:14 CEST 2014 : added initialization functions to automatically set all entries in an array to 'val'

//  changes from  Sun Apr 13 23:45:29 CEST 2014 : added VNN_gen_smallp for the SNG and LNG interactions so that the void equations and GA can import from the same source.

//  changes from  Thu Nov 21 15:05:27 CET 2013 : fixed coarse_grain function to read only up to x<L, not x<=L; This resolved the previously observed segfault

//  changes from  Wed Nov 6 15:37:25 CET 2013 : Added interpolate general for arbitrary input arrays, and pick_index_from_norm_dist. Both do what might be expected from the name
//----------------------------------------------------------------------------

 #ifndef __brenlib_STANDARDS  //---check whether it's been defined already so we don't do it twice. 
 #define __brenlib_STANDARDS 

#include <fstream>   
#include <iostream>  
#include <iomanip>  
#include <sstream>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <queue>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_pow_int.h>
#include <sstream>

//*******************************************************************************************
using namespace std;

//--------- HERE IS A SET OF FUNC'S TO INITIALIZE AN ARRAY TO A CERTAIN VALUES  ------
//---------------- updating a single line to test git fetch on Nov. 21-----------------------------------------------------


//*************************************************************************

int clear_charray(char* c, const int charlength )
{
int i;
for(i=0;i<charlength;i++)
	{
	c[i]='\0';
	}
return 1;
}

//-----------------------------------
int init_array( int * A, const int length, const int val )
{
int i;
if (length <= 0)
	{
	cout << "\n ERROR: initializing array with zero size.\n";
	exit(1);
	}

for(i=0; i<length ; i++)
	{
	A[i]=val;
	}

return A[length-1];
}
//-----------------------------------
int init_array( double * A, const int length, const double val )
{
int i;
if (length <= 0)
	{
	cout << "\n ERROR: initializing array with zero size.\n";
	exit(1);
	}

for(i=0; i<length ; i++)
	{
	A[i]=val;
	}

return A[length-1];
}

//-----------------------------------
int init_array( bool * A, const int length, const bool val )
{
int i;
if (length <= 0)
	{
	cout << "\n ERROR: initializing array with zero size.\n";
	exit(1);
	}

for(i=0; i<length ; i++)
	{
	A[i]=val;
	}

return A[length-1];
}
//**************************************************************************
bool DirectoryExists( const char* pzPath )
{
    if ( pzPath == NULL) return false;

    DIR *pDir;
    bool bExists = false;

    pDir = opendir (pzPath);

    if (pDir != NULL)
    {
        bExists = true;    
        (void) closedir (pDir);
    }

    return bExists;
}
//**************************************************************************
double interpolate_mu_from_rhoi(const double irho_target,const int w, const double E0, const int kHNG, const string NGtype, const string CGF_feedin_string)
{
int j;
int charlength=400; //--the number of characters in the string for our path.
char cpath[charlength];
ifstream fin;
double result=0;

clear_charray(cpath, charlength );

if(NGtype == "HNG")
	{
	sprintf(cpath, "mu_v_irhoHNG_k%d_CGF_%s.txt",kHNG,CGF_feedin_string.c_str());
	}
else if(NGtype == "SNG")
	{
	sprintf(cpath, "mu_v_irhoSNG_w%d_E%.4f_CGF_%s.txt", w, E0,CGF_feedin_string.c_str());
	}
else if(NGtype == "LNG")
	{
	sprintf(cpath, "mu_v_irhoLNG_w%d_E%.4f_CGF_%s.txt", w, E0,CGF_feedin_string.c_str());
	}
else
	{
	cout << "\n ERORR: undefined NGtype in interpolate_mu_from_rhoi. Exiting.\n\n";
	exit(1);
	}

fin.open(cpath);

if( fin.fail() )
	{
	cout << "\n ERROR: cannot locate interpolation file at path = " << cpath << "\n. exiting.\n";
	exit(1);
	}

//-----------------------------------------------------
int i=0, N=0;
double dummy1, dummy2;
fin >> dummy1 >> dummy2;

while(!fin.eof())
	{
	N++;
	fin >> dummy1 >> dummy2;
	}
fin.close();
fin.open(cpath);
//-----------------------------------------------------
double muarr[N];
double irhoarr[N];

for(i=0;i<N;i++)
	{
	fin >> muarr[i] >> irhoarr[i];
	}
fin.close();
//-----------------------------------------------------

bool found=false;
for(i=0;i<N-1;i++)
	{
	
	if( irho_target == irhoarr[i] )
		{
		result = muarr[i];
		found=true;
		break;			
		}
	else		
		{	
		if( irho_target < irhoarr[i] && irho_target > irhoarr[i+1] ) 
			{
			result = muarr[i] + (irho_target-irhoarr[i]) * ( (muarr[i+1]-muarr[i]) / ( irhoarr[i+1] - irhoarr[i] ) )    ; 

			// obviously (j+2)-(j+1)=1, but this factor is left 
			// in there to make the interpolation more transparent.
			found=true;
			break;
       			}
		
		}
	}

  
if(!found)
	{
	cout << "\n CRITICAL ERROR!, target irho value = " << irho_target << " is beyond the range of the interpolation file. exiting. \n\n";
	exit(1);
	}
else
	return result;

}

//*****************************************************************
double bren_array_dotprod(const int L, vector<double>&  A, vector<double>& B )
{
int i;
double NA=0.0, NB=0.0;
double result=0.0;


double normcheckA=0.0;
double normcheckB=0.0;

for(i=0;i<L;i++)	//---first get the vector norms for each.
	{
	NA+= pow(  A.at(i), 2);
	NB+= pow(  B.at(i), 2);
	}
NA=sqrt(NA);
NB=sqrt(NB);

//! cout << "\n norm A= " << NA << " norm B=" << NB << endl;

//----------------now multiply them out:-------------------------

for(i=0;i<L;i++)	//---first get the vector norms for each.
	{
	result += (A.at(i)/NA)*(B.at(i)/NB);

	normcheckA += (A.at(i)/NA)*(A.at(i)/NA);
	normcheckB += (B.at(i)/NB)*(B.at(i)/NB);
	}

if( fabs(normcheckA-1.0) > 1E-10  ||  fabs(normcheckB-1.0) > 1E-10 )
	{
	cout << "\n ERROR: normchecks don't work out. exiting.\n\n ";
	exit(1);
	}



return result;
}

//****************************************************************************

int space_tpoints_linear(const double tmin, const double dt_obs, const double tmax, const int nbins, double * tx)
{
int i=0;
double temp=tmin;
int n=0;

i=0;
while( i < floor((tmax -tmin)/dt_obs)    )
	{
	temp = tmin + (double(i)*dt_obs);
	tx[i]=temp;

	i++;
	}
	
if (i != nbins )
	{
	cout << "\n ERROR in space tpoints_linear : nbins =" << nbins << ", i= " << i << endl;
	exit(1);
	}


return i;

}
//****************************************************************************
int space_tpoints_logarithmic(const double t0, const double tf, const int Np10, const int nbins, double * tx)
{
int i;
double l10ti; //----log base 10 of the time point 'i' -raise 10 to it to get the time point.

double l10t0 = log10(t0);
if (isnan(l10t0) || isinf(l10t0) )
	{
	cout << "\n ERROR: in space_tpoints_logarithmic, log of initial point is NaN or inf, exiting. \n";
	exit(1);
	}
for(i=0;i<nbins;i++)
	{
	l10ti   = l10t0 + (double(i)/double(Np10));
	tx[i] = pow (10.0, l10ti );
	}
	
//---- nbins is defined as = ceil(   Np10* log10( tf/t0)  );
//---- in the gillespie function that calls this.

return 1;
}

//****************************************************************************
double get_local_array_avg(const double *x, const int centerpos, const int range)
{
double result =0.0;
int i;

for( i=centerpos - range; i <= centerpos+range; i++)
	{
	result +=x[i];
	}
result = result/(2*range+1);

return result;
}

//****************************************************************************
double get_array_max( const double *x, const double *y, const int size)
{
double result =0.0;
int i;

//---not bothering with regression for the moment. double c0=0.0, c1=0.0, c2=0.0;

for( i=0; i < size; i++)
	{
	if(y[i] > result)
		{
		result = y[i];
		}
	}

return result;
}

//*****************************************************************
int VNN_SNG_calc(double * potential,const int w, const double E0)
{ // potential has allocated size [NNRANGE_full+1], ( NNRANGE_full =2*w), hence 2*w+1


int p = 2*w+1;
double Omega;
int n,i,j,x,bound;

double y[p];

//--- THE INDEX STARTS AT ZERO: i.e. WE ENUMERATE THE EMPTY LATTICE SITES THAT EXIST
//--- IN BETWEEN THE POINTS OF INTEREST: VNN[0] IS THE INTERACTION POTENTIAL OF 
//--- TWO NUCLEOSOMES THAT ARE DIRECTLY ADJACENT.

for(n=0;n < p;n ++) { potential[n]=0.0;  y[n]=0.0; }  // INITIALIZE	

for(n=0 ; n<p ; n++) 
	{ 
	if(n <= w)	// n, and not n-1 here since '0' is the starting index
	   {
	   for(i=0;i<=n;i++)
	      {
	      for(j=0;j<=(n-i);j++)
		{	
		y[n] += gsl_sf_exp(E0*double(i+j)); 
		}
	      }
	   } 

	else 
	   for(i=0;i<=w;i++)
	      {
	      bound=min(w, n-i);
	      for(j=0;j<=bound;j++)
		{	
		y[n] +=  gsl_sf_exp(E0*double(i+j))  ; 
	        }
	      }
	}

Omega=0.0;
for(i=0;i<=w;i++)
	{
	for(j=0;j<=w;j++)
	   {	
	   Omega += gsl_sf_exp(E0*double(i+j));
	   }
	} 

for(n=0;n<p;n++)
	{ 
	potential[n] = -1.0*gsl_sf_log(y[n]/Omega); 
	}



return 1;
}
//*****************************************************************
int VNN_LNG_calc(double * potential,const int w, const double E0)
{// potential has allocated size [NNRANGE_full+1], ( NNRANGE_full =2*w), hence 2*w+1

int a = 2*w+1;
double Omega;
int n=0;

for(n=0 ; n<a ; n++) //the max value is (a-1), and there the potential becomes zero.
	{
	potential[n] = (a-(n+1))*E0;
	}

return 1;
}
//*****************************************************************
int VNN_SNG_calc_smallp(double * potential, const int NNRANGE, const int a, const double rm, const double E0)
{// potential has allocated size [NNRANGE_full+1], ( NNRANGE_full =2*w), hence 2*w+1
 //-- HERE, 'a' is the footprint size, and 'rm' is the min-potential distance IN UNITS OF A
 //-- (i.e. a*rm) is the point of minimum interaction, which is generally not an integer


double r;
double Omega;
int n,i,j,x,bound;
int w=NNRANGE/2;

bool LJ = false;
double y[a];

if (LJ)
	{
	cout << "invoking the Lennard-Jones potential.\n";

	for(n=0 ; n<NNRANGE ; n++) // a is the particle minimum, NNRANGE is how far out we calculated it (_almost_ everywhere).
		{
		//---	potential[n] = E0*( pow( ( rm*double(a)/double(n+1) ),12.0) - 2.0 * pow((rm*double(a)/double(n+1) ), 	6.0) );
		//--- this was the old Lennard-Jones potential (included attractive minimum... here's the new one:
	
		r=(double(n+1)/double(NNRANGE));
		potential[n] = E0*(gsl_sf_pow_int ((1.0/r), 12) - 2.0*gsl_sf_pow_int ((1.0/r), 6) +1);
		}
	}

else	//-------------------- JUST THE REGULAR SNG POTENTIAL ---------------------------
	{


	//--- THE INDEX STARTS AT ZERO: i.e. WE ENUMERATE THE EMPTY LATTICE SITES THAT EXIST
	//--- IN BETWEEN THE POINTS OF INTEREST: VNN[0] IS THE INTERACTION POTENTIAL OF 
	//--- TWO NUCLEOSOMES THAT ARE DIRECTLY ADJACENT.

	for(n=0; n < a;n ++) { potential[n]=0.0;  y[n]=0.0; }  // INITIALIZE	

	for(n=0 ; n<a ; n++) 
		{ 
		if(n <= w)	// n, and not n-1 here since '0' is the starting index
		   {
		   for(i=0;i<=n;i++)
		      {
		      for(j=0;j<=(n-i);j++)
			{	
			y[n] += gsl_sf_exp(E0*double(i+j)); 
			}
		      }
		   } 

		else 
		   for(i=0;i<=w;i++)
		      {
		      bound=min(w, n-i);
		      for(j=0;j<=bound;j++)
			{	
			y[n] +=  gsl_sf_exp(E0*double(i+j))  ; 
		        }
		      }
		}

	Omega=0.0;
	for(i=0;i<=w;i++)
		{
		for(j=0;j<=w;j++)
		   {	
		   Omega += gsl_sf_exp(E0*double(i+j));
		   }
		} 

	for(n=0;n<a;n++)
		{ 
		potential[n] = -1.0*gsl_sf_log(y[n]/Omega); 
		}


	}

//---------- JUST DEFAULT BACK TO THE SNG POTENTIAL


return 1;
}
//*****************************************************************
int VNN_LNG_calc_smallp(double * potential, const int a, const double E0)
{// potential has allocated size [NNRANGE_full+1], ( NNRANGE_full =2*w), hence 2*w+1
 // the index refers to the number of empty intervening lattice sites., hence v[0] is for adjacent particles.

int n=0;
for(n=0 ; n<a ; n++) //the max value is (a-1), and there the potential becomes zero.
	{
	potential[n] = (a-(n+1))*E0;
	}

return 1;
}

//****************************************************************
int coarse_grain(const double * Vin_full, double * xcoarse, double * Vout_coarse, const int L, const int p, const double CGF) //--coarse-grain the 2-body interaction potential into a smaller system.
{
// -----here, 'p' is the length of the original interaction. 
// -----'L' is the length that you want to calculate the coarsened interaction out to (if it is longer than CGF*p, then the extra positions will just be filled in with zeros)
// -----xcoarse is the array of positions along the x-axis after coarsening.

//-- WARNING -Vin MUST be spaced by exactly 1, and Vin[0] must correspond to the potential energy of direct neighbours (i.e. gap of zero, distance of 1), hence V[n] means 'energy inherent in a gap of n, or equivalently a distance of n+1'

int i,j,n;

for(j=0;j<L;j++) //-----the array is only of size 'L' -the L'th position should map back to zero.
	{
	xcoarse[j]=(j+1)*CGF;
	}
//-------------------------------------now coarse-grain the y-potential----------------------------

bool found;
for(i=0;i<L;i++) //----cycle through the 'i' positions on the CG'ed lattice looking for matches 
	{
	found=false;
   
	if(i*CGF>=p)	//--- we just automatically set the potential to zero and skip this case to avoid array read-over.
		{
		Vout_coarse[i]=0.0;
		found = true;
		}
	else
	   {
	   for (j=0;j<p;j++) // ---- cycle through the full lattice sites 'j', and check if we can find that x-value EXACTLY,
		{	     // ---- if so, then just set the coarse-value equal ---

		if( xcoarse[i] == (j+1) )
			{
			Vout_coarse[i] = Vin_full[j];
			found=true;
			break;	
			}
		}

	   if(!found) // ----IF NOT THEN LOOK FOR A BOUNDING REGION AND INTERPOLATE.
		{	
		for (j=0;j<p;j++) 
			{
			if( xcoarse[i] > (j+1) && (xcoarse[i] < (j+1+1)) ) 
				{
				Vout_coarse[i] = Vin_full[j] + (xcoarse[i]-(j+1))*((Vin_full[j+1]) - Vin_full[j]/((j+2)-(j+1)) )    ; // obviously (j+2)-(j+1)=1, but this factor is left 
				// in there to make the interpolation more transparent.
				found=true;
				break;
           			}
			}
		}


   //---- the only remaining possibility is that it's the last one and it's
   //---- beyond the interaction range -BUT THIS SHOULD ONLY HAPPEN ONCE! -thats
   //---- why we check for if i==xsizez_coarse
   
	if( xcoarse[i] > p )
		{
		found=true;
		Vout_coarse[i]=0.0;
		}
  
	if(!found)
		{
		cout << "\n error, didnt find a good interpolation\n\n";
		}
	  }//---end the "else" on whether i > p
	}//end the for-loop through i.
	

return 1;
}


//******************************************************************************
string bren_itoa( const int x )
{
string result;

std::stringstream ss; 
ss << x;
result = ss.str();

return result;
}

//******************************************************************************
double max(const double a, const double b)
{
if(a>b)
	return a;
else
	return b;
}
//=================================
double min(const double a, const double b)
{
if(a<b)
	return a;
else
	return b;
}
//=================================
int min(const int a, const int b)
{
if(a<b)
	return a;
else
	return b;
}
//=================================
int max(const int a, const int b)
{
if(a>b)
	return a;
else
	return b;
}


//******************************************************************************

double interpolate_general( const double * x1, const double * y1, const int size_in, const double * x2, double * y2_out, const int size_out)
{ 
// --- takes x1,y1, arrays of size size_in, and assumes linear dependance to 
// --- interpolate y2_out according to x2 points, each of the latter of which are of size size_out

int i, j; //---counter indices
int num_x1_points = size_in;
int num_y1_points = size_in;
int num_x2_points = size_out;
double result=0.0;

double x1_min, x1_nexttomin, y1_min, y1_nexttomin, x1_max, x1_nexttomax, y1_max, y1_nexttomax;

for ( i=0; i<=(num_x1_points-2); i++ )
	{	
	if( (x1[i+1]-x1[i])*(x1[i+2]-x1[i+1]) <0  ) 
		{
		cout << "\n ERROR : x1 input values must be monotonic";
		}
	}

int direction =1;

if( x1[2]-x1[1] <0)
	{
	//------------- y1 - DECREASING ----------------------
	 cout << "\n WARNING: in general interpolation file. order is inverted.\n"; 
 	direction = -1;
	x1_min         = x1[num_x1_points-1];
	x1_nexttomin   = x1[num_x1_points-2];
	y1_min         = y1[num_x1_points-1];
	y1_nexttomin   = y1[num_x1_points-2];

	x1_max         = x1[0];
	x1_nexttomax   = x1[1];
	y1_max         = y1[0];	
	y1_nexttomax   = y1[1];

	}
else
	{

	//------------- y1 - INCREASING ----------------------
	x1_min         = x1[0];
	x1_nexttomin   = x1[1];
	y1_min         = y1[0];
	y1_nexttomin   = y1[1];

	x1_max         = x1[num_x1_points-1];
	x1_nexttomax   = x1[num_x1_points-2];
	y1_max         = y1[num_x1_points-1];	
	y1_nexttomax   = y1[num_x1_points-2];
	}


//-------- GOT THE X-AXIS, NOW GET THE CG'ed potential

bool found=false;

for (i=0; i < num_x2_points; i++)
	{ 
	found=false;

//-----BOUNDARY CASES--------   
	if( x2[i] <= x1_min )	//----- we just automatically set the potential to zero and skip this case to avoid array read-over.
		{
		found = true;
	        y2_out[i] = y1_min + (x2[i]-x1_min)*( (y1_min - y1_nexttomin)/(x1_min - x1_nexttomin ) )    ;
        
		}
	else if( x2[i] >= x1_max)
		{
		found = true;
        	y2_out[i] = y1_max + (x2[i]-x1_max)*( (y1_max - y1_nexttomax)/(x1_max-x1_nexttomax) )  ;
		}
        else
        	{
	//------SOMEWHERE WITHIN THE BOUNDARIES------------
        	for (j=0; j<(num_x1_points-1);j++)
			{
           		if(found)
				{ 
				cout << "\n ERROR in interpolating_general. repeating after found\n";
				 exit(1);
				}
			if( x2[i] >= x1[j] && x2[i] < x1[j+1] )
				{
				y2_out[i] = y1[j] + (x2[i]-x1[j])*((y1[j+1] - y1[j])/(x1[j+1]-x1[j]) ) ;
				found=true;
				break;
 				}
           		}

        	}
            
	if(!found)
		{
		cout << "\n ERROR : couldnt find an interpolation point at x2[ " << i << "] = " << x2[i] << endl;
		exit(1); 
		}
  	}

result = 0.0;
for (i=0; i < num_x2_points; i++)
	{ 
	result += y2_out[i];
	}
// ------- end the for-loop through i.

return result;	    
}
//******************************************************************************
int pick_index_from_rand_norm_dist( double rand01, double * distr, int L )
{
//  takes a normalized probability distribution distr, of length L, 
//  and selects an index from a uniformly distributed random input variable rand01

double sum = 0.0;
bool found = false;
int result;
int i;

for(i=0;i<L;i++)
	{
	sum += distr[i];

	if(  (sum > rand01) && (found == false)  )
		{
		result = i;
		found  = true;
		}
	}

if ( fabs(sum-1.0) > 0.000001)
	{
	cout << "\n ERROR in pick_index_from_rand_dist: distribution is not normalized.\n";
	exit(1);	
	}

if (!found)
	{
	cout << "\n ERROR in pick_index_from_rand_dist. Failed to make a selection found=false.\n";
	exit(1);
	}

else 
	{
	return result;
	}
}

//********************************************************************************

int voidsize_corr_inc_states_from_string(string config ,int repnum, gsl_matrix * M2, gsl_vector * v1, const int L )
{
int N=0,i=0; 
double prev_density=0.0;
vector<int> voidvec;
string delim1 = "_";
string delim2 = "-";


//--- number of particles (equivalently, voids).
//--- example taken from : http://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c 

size_t pos = 0;
string token;

pos   = config.find(delim1);
if ( pos == std::string::npos )
	{  //---- empty string with just 0
	if(config!="0")
		{
		cout << "\n ERRORin voidsize_corr_inc_states_from_string: string not 0 but still not finding delim1 exiting.\n ";
		exit(1);
		}
	voidvec.push_back(L);
	}
else
	{
	token = config.substr(0, pos);
	N     = atoi(token.c_str());//--- number of particles
	config.erase(0, pos + delim1.length());

	while ((pos = config.find(delim2)) != std::string::npos) 
		{
		token = config.substr(0, pos);
		voidvec.push_back(atoi(token.c_str()));
		config.erase(0, pos + delim2.length());
		}

	//----- THEN GRAB THE LAST ONE ------
	token = config.substr(0, config.size());
	voidvec.push_back(atoi(token.c_str()));
	// config.erase(0, delim2.length());
	}

//-------------- NOW THE STRING CHARACTERS ARE CONVERTED INTO A VECTOR IN ORDER ------------

int sum=0;
for(i=0; i < voidvec.size(); i++)
	{
	sum+=voidvec.at(i);
	}
if (sum != L)
	{
	cout << "\n ERROR in voidsize_corr_inc_states_from_string. Voids don't sum to L.";
	exit(1);
	}




int Vi;
int Vip1;

for(i=0; i < (voidvec.size()-1); i++)
	{
	Vi   = voidvec.at(i);
	Vip1 = voidvec.at(i+1);

	prev_density = gsl_vector_get(v1,Vi);
	gsl_vector_set(v1,Vi, (prev_density+repnum) );

	prev_density = gsl_matrix_get(M2,Vi,Vip1);
	gsl_matrix_set(M2,Vi,Vip1,(prev_density+repnum) );

	}

Vi   = voidvec.at( voidvec.size()-1);
Vip1 = voidvec.at(0);		//----periodic boundaries

prev_density = gsl_vector_get(v1,Vi);
gsl_vector_set(v1,Vi,prev_density+(repnum) );

prev_density = gsl_matrix_get(M2,Vi,Vip1);
gsl_matrix_set(M2,Vi,Vip1,prev_density+(repnum) );


//-------------- NOW THE OCCURANCE VALUES IN THE MATRIX/VECTORS HAVE BEEN INCREMENTED -------


}
//****************************************************************************888
int bren_print_gsl_matrix_to_stream(const gsl_matrix * M, ostream &out)
{
int i,j,nrow,ncol;

nrow=M->size1;
ncol=M->size2;

out << endl;
out << endl;

for(i=0;i<ncol;i++)
	{
	for(j=0;j<nrow;j++)
		{
		out << gsl_matrix_get(M,i,j) << " \t ";
		}
	out << endl;
	}

return 1;
}



//**********************************************************************************************

int voidsize_corr_get_configs_from_file (std::string filepath, std::string prefix, const int time_index, gsl_matrix * M2, gsl_vector * v1, const int numruns, const int L)
{//---- matrix of neighbour voids M, and vector of single voids v should be allocated and initialized in main.

const int charlength =400;	//----the length of the string kept in memory for file-path manipulation.
string tempstring;
int    temprep=0,i=0;

char cpath[charlength];
clear_charray(cpath, charlength );
sprintf(cpath, "%s%s_%d.txt",filepath.c_str(),prefix.c_str(),time_index); //----input file should always be in the working directory.

ifstream fin( cpath );
if (fin.fail())
	{
	cout << "\n ERROR in voidsize_corr_get_configs_from_file, failed to open file. exiting. \n";
	exit(1);
	}

vector<string> configs; //---- vector of independent configurations that are possible.
vector<int>  rep;	//---- number of times such a configuration was repeated. 

int Numconfigs=0;
fin >> tempstring >> temprep; 

while(!fin.eof())
	{
	Numconfigs++;
	configs.push_back (tempstring);
	rep.push_back     (temprep);

	fin >> tempstring >> temprep; 
	}

//--------- CROSSCHECK ---------
int sum=0;
for (i=0; i<Numconfigs; i++ )
	{
	sum+=rep.at(i);
	}
if (sum != numruns)
	{
	cout << "\n ERROR in voidsize_corr_get_configs_from_file. reps don't sum to numruns.";
	exit(1);
	}
//------------------------------

for (i=0; i<Numconfigs; i++ )
	{
	tempstring = configs.back();
	configs.pop_back();

	temprep    = rep.back();
	rep.pop_back();

	voidsize_corr_inc_states_from_string( tempstring , temprep, M2, v1, L );
	}


// This function multiplies the elements of vector a by the constant factor x. The result a_i \leftarrow x a_i is stored in a. 

}
//*************************************************************************


#endif //--this ends the clause as to whether or not this symbol (i.e. this file) has been defined already.
