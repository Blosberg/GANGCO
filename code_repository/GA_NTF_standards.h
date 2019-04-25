//--------GA_NTF_standards.h -contains a few functions useful for the GC implementation of the nucl GA.-----
// ---last updated on  Fri Jan 17 18:27:38 CET 2014  by  Brendan.Osberg  at location  th-ws-e537

//  changes from  Fri Jan 17 18:27:38 CET 2014 : added conditional in check_states so that n is incremented only if krm_b==true --controls number of remodelling candidates.

//  changes from  Fri Nov 29 10:52:02 CET 2013 : added remod_cand. function and int min_dist to GAdata; this way remodellers know when they've hit neighbours for the case of HNG.

//----------------------------------------------------------------------------


#ifndef GA_STANDARDS  //---check whether it's been defined already so we don't do it twice. 
#define GA_STANDARDS

#include <fstream>   
#include <iostream>  
#include <iomanip>  
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <queue>
#include <math.h>
#include <gsl/gsl_rng.h>


//*****************************************************************
int VNTF_calc(double * potential,const int w, const double E0)
{	//-- REMEMBER, THE POTENTIAL NOW REFERS TO THE NUMBER OF EMPTY LATTICE SITES
	//-- IN BETWEEN, NOT THE NUMBER OF INTERSTITIALS.
double Omega;
int n,i,x;

double y[w];

for(n=0;n < w;n ++) { potential[n]=0.0;  y[n]=0.0; }  // INITIALIZE	

for(n=0 ; n<w ; n++) 
	{ 
	   for(i=0;i<=n;i++)
	      {
		y[n] += gsl_sf_exp(E0*double(i)); 
	      }

	}

Omega=0.0;
for(i=0;i<w;i++)
	{
	Omega += gsl_sf_exp(E0*double(i));
	} 

for(n=0;n<w;n++)
	{ 
	potential[n] = -1.0*gsl_sf_log(y[n]/Omega); 
	}

return 1;
}

//*****************************************************************
double GAdata::muTF(int x )
{
/*-----------------------
if((x<F0) || (x >= (F0+F1)))
	 return muTF0;
else if((x>=F0)  && (x<(F0+F1)) )
	 return muTF1;
else
	{
	*(log) << "\n\n unable to assign muTF value, exiting. \n\n ";
	(*log).close();
	exit(1);
	}
if(x==F0)
	{
   return muTF0;	
	}
else if(x==F1)
   {
   return muTF1;
   }
else
	{
	return muTF_nonspec;
	}
------------------------*/

if(x==F0 || x==F1 )
   return muTF1;	
else
   return muTF0;	

}

//***************************************************************************************
int react(int x,int Rtype, GAdata &P)
{
switch(Rtype)
	{
	case 0: P.remove_Nuc( x ); 	break;
	case 1: P.add_Nuc(x); 		break;
	case 2: P.slide_Nuc_left( x); 	break;
	case 3: P.slide_Nuc_right( x); 	break;

	case 4: P.remove_TF( x); 	break;
	case 5: P.add_TF( x); 		break;
	case 6: P.slide_TF_left( x); 	break;
	case 7: P.slide_TF_right( x); 	break;

	default:

	P.flag=102;
	P.recentRx=-10;
	P.process_error(x);	

	}
}

//******************************************************************************
int GAdata::check_states( void )
{
int i = 0;
int x = 0;
int first_part, last_part;
int pL,pR;

first_part = pos[Llim-1].part_right;
last_part  = pos[0].part_left;  // either way.

bool wrapped=false;


// *log << "counter = " << counter << endl;

x=0;

if(pos[0].state == 0)
	{
	pR=first_part; 
	}
else
	{
	x  = first_part; 
	pR = x;
	}

pL=last_part;

if(partnum ==0)
	{
	for(i=0;i<Llim;i++)	
		{
		if(  pos[i].state == 2 || pos[i].state ==1 || pos[i].part_left != i || pos[i].part_right != i )
		    {
		    flag = 2022;
		    process_error(i);
		    }
		}
	
	}

else
  {
  while(1)
	{
	if(pos[x].state ==0)//------------------------
	   {
	   pL = pos[x].part_left;
	   pR = pos[x].part_right;
	
	   i=right(pL);
	   while(1)
		{
		if(i == pR )
			break;
		else if(  pos[i].state == 2 || pos[i].state ==1 || pos[i].part_left !=pL || pos[i].part_right != pR )
		    {
		    flag = 2000;
		    process_error(i);
		    }
		else
		    i=right(i);
		}

	   if(i<x)//wrapped around
		{
		wrapped=true;
		break;
		}
	   else
		x=i; 
	   } //---end of 'if(pos[x].state=0)'
	if(pos[x].state ==1)//-----------------------------------
	   {
	   if( (x==first_part) && wrapped  )
		break;

	   if( pos[x].part_left != pL   || x != pR) 
		    {
		    flag = 2001;
		    process_error(x);
		    }
	
	   pR = pos[x].part_right; // now change pR
	   pL = x;
	
	   x=right(x);
	   if(x==0)	
		wrapped=true;
	   }//---end of 'if(pos[x].state=1)'
	if(pos[x].state ==2)//----------------------------------
	   {
	   if( (x==first_part) && wrapped  )
		break;

	    if( pos[x].part_left != pL   || x != pR) 
		    {
		    flag = 2002;
		    process_error(x);
		    }
	
		
	   pR = pos[x].part_right; // now change pR
	   pL = x;
	
           x=right(x);	
	   if(x==0)
	       wrapped=true;	  

	   for(i=1;i<m;i++)
		{
		if( (pos[x].state != 2+i) || (pos[x].part_left !=pL ) || (pos[x].part_right != pR) )
			{
			flag = 2003;
			process_error(x);
			}

		x=right(x);	
		if(x==0)
			wrapped=true;
		}

	   }//---end of 'if(pos[x].state=2)'
		if((pos[x].state >2)  &&  (pos[x].state != (2+distance(pos[x].part_left,x)) ) ) 
		{
	   	flag = 2004;
		process_error(x);
		}	
	}
 } //---end of 'else' from 'if(partnum=0)' -i.e. the contingency that there are particles.


int n=0;//--counter for the number of valid Nucl. pairs.
x=first_part;
int N_remod_test=0;

while( krm_b  &&  Nucnum >= 2  )
  {
   pR = pos[x].part_right;
   if( pos[x].state ==1  && pos[pR].state ==1 && remod_cand(x,pR) )	
	{ //---ARE THEY BOTH NUCL'S AND ARE THEY CANDIDATES FOR REMODELLING?
	n++;
	}
   if( pR == first_part) 
	{
	break;
	}
   else
	{
	x=pR; //keep going.
	}

  }

if(n!=N_remod)	//----  having counted all the way around, n should now be identical to the number of 
   {		//----  pairs that are eligible for remodelling.
   
   flag=102356;
   recentRx=-20;
   process_error(x);	
   exit(1);
   }

return 1;
} //----end of function

//*****************************************************************
int GAdata::cleanup_occ( void)
{
int x=0,j=0,k=0;

for (x=0;x<Llim;x++)
	{
	if(pos[x].state == 1)
		{
		j++;
		if( (pos[x].t_bind_N <0) || (pos[x].t_bind_N > t) )
			{ 
			flag = 4011; 
			process_error(x);
			}
		else
			{
			pos[x].t_occ_N += (tf - pos[x].t_bind_N);
			pos[x].t_bind_N = tf;
			}
		}
	else if( pos[x].state == 2 )
		{
		k++;
		if( (pos[x].t_bind_TF <0) || (pos[x].t_bind_TF > t) )
			{ 
			flag = 4012; 
			process_error(x);
			}
		else
			{
			pos[x].t_occ_TF += (tf - pos[x].t_bind_TF);
			pos[x].t_bind_TF = tf;
			}
		}

	}
if((j!= Nucnum) || (k!=TFnum))
	{
	flag=4013;
	process_error(1);
	}

return (j+k);
}
//*****************************************************************
int GAdata::reset( void)
{

// this will yield a duplicate function definition next time you try to compite with this.
// you copied it into the GA_absolute_standards file to avoid the TF def.

int x;
   Nucnum  = 0;
   TFnum   = 0;
   partnum = 0;
   a0=0.0;
   t=0.0;
   counter=0;

//------HERE WE SET THE INITIAL REACTION RATES----------------
for(x=0;x<Llim;x++)	
	{
	
	pos[x].state=0;		// initialize each lattice pos to 0 (nothing)
	pos[x].part_right = x;
	pos[x].part_left  = x;


	*pos[x].a_removeN = 0.0;
	*pos[x].a_sNL = 0.0;
	*pos[x].a_sNR = 0.0;

	*pos[x].a_removeTF = 0.0;
	*pos[x].a_sTFL = 0.0;
	*pos[x].a_sTFR = 0.0;

	*pos[x].a_addN  = ka_N*gsl_sf_exp(muN(x)); 	//-- ADDING NUCLEOSOMES  -------------
	*pos[x].a_addTF = ka_TF*gsl_sf_exp(muTF(x));	//---ADDING TRANSCRIPTION FACTORS -------	


	a0+=*pos[x].a_addN;
	a0+=*pos[x].a_addTF;

	//-----------------AND SET AVERAGE OCCUPANCY TO ZERO INITIALLY------------

	pos[x].t_occ_TF = 0.0;
	pos[x].t_occ_N  = 0.0;

	pos[x].t_bind_TF = -1.0;
	pos[x].t_bind_N  = -1.0;
	} //---end for x=0..Llim	

}

//**********************************************************************
int GAdata::plot_snapshot(void)
{
int i,j,x;
char cpath[200];
ofstream fout;

int N_there_now=0;
int TF_there_now=0;

for(j=0;j<200;j++)
	{
	cpath[j]='\0';
	}

//----------------------------------------------------------------
if(plotnum < 10)
	{
	sprintf(cpath, "%sdens_snap_index_0000%d.txt", path.c_str(),plotnum);
	}
else if(plotnum < 100)
	{
	sprintf(cpath, "%sdens_snap_index_000%d.txt", path.c_str(),plotnum);
	}
else if(plotnum < 1000)
	{
	sprintf(cpath, "%sdens_snap_index_00%d.txt", path.c_str(),plotnum);
	}
else if(plotnum < 10000)
	{
	sprintf(cpath, "%sdens_snap_index_0%d.txt", path.c_str(),plotnum);
	}
else 
	{
	sprintf(cpath, "%sdens_snap_index_%d.txt", path.c_str(),plotnum);
	}

fout.open(cpath);
double current_t_occN =0.0;
double current_t_occTF=0.0;

//--------------------------------------------------------------------------------------------
for(x=0;x<Llim;x++)
	{
	current_t_occN =0.0;	
	current_t_occTF=0.0;	

	if(pos[x].state==1)
		{
		N_there_now=1;
		current_t_occN += (t - pos[x].t_bind_N);
		}
	else
		N_there_now=0;

	if(pos[x].state==2)
		{
		TF_there_now=1;
		current_t_occTF += (t - pos[x].t_bind_TF);
		}
	else
		TF_there_now=0;


	current_t_occN  += pos[x].t_occ_N;
	current_t_occTF += pos[x].t_occ_TF;

	fout << x << " \t " << (1.0/t)*current_t_occN << " \t " << (1.0/t)*current_t_occTF << " \t " << N_there_now << " \t " << TF_there_now << endl;
	}

fout.close();
*timestamps  << t << endl; //---print out the time of this plot


plotnum ++;

return 1;
}
//**********************************************************************
int GAdata::plot_snapshot_kymo( ofstream * foutN, ofstream * foutTF )
{
int i,j,x;
char cpath[200];

int N_there_now=0;

for(j=0;j<200;j++)
	{
	cpath[j]='\0';
	}
bool foundN=false;
bool foundTF=false;

//--------------------------------------------------------------------------------------------
for(x=0;x<Llim;x++)
	{

	if(pos[x].state==1)
		{
		*foutN << t << " \t " << x << endl;
		foundN=true;
		}

	if(pos[x].state==2)
		{
		*foutTF << t << " \t " << x << endl;
		foundTF=true;
		}
	}

if(foundN)
	{
	*foutN << endl;
	}

if(foundTF)
	{
	*foutTF << endl;
	}


plotnum_kymo ++;

return 1;
}

//****************************************************************************************


#endif //--this ends the clause as to whether or not this symbol (i.e. this file) has been defined already.
