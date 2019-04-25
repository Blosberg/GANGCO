//----------------------------------------------------------------------------
//---LAST UPDATED Friday SEPTEMBER 30th 2011 11:55PM

//---THESE ARE SOME GENERAL FUNCTIONS USED BY ALL IMPLEMENTATIONS OF THE GA SIMULATION
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
double	GAdata::dEsL(int type, int x, int * neighbours )
{
double iL=0.0;	// old interaction energy with the left particle
double iR=0.0;	// old interaction energy with the right particle
double iLp=0.0;	// new interaction energy with the left particle
double iRp=0.0;	// new interaction energy with the right particle

double mux=0.0, muxp=0.0;


mux  = muN(x);
muxp = muN(left(x));

if(type != 1)      //-----IN THE CANONICAL ENSEMBLE IT SHOULD BE A NUC.
	{	
	flag=6;
	interaction_error(type, x, neighbours  );
	}
	//--------------------------------------------------------------------
if(Nucnum==1)
	{
	iLp = 0.0;
	iRp = 0.0;
	iL  = 0.0;
	iR  = 0.0;
	}
else 
  {
  if ( type == 1)	//----*THIS* PARTICLE IS A NUCL----
	{

		iL  = interaction_NN( distance(neighbours[0], x ) );
		iLp = interaction_NN( distance(neighbours[0], left(x) )  );

	if(pos[neighbours[0]].state != 1)
		{	
		flag=1;	
		interaction_error(type, x, neighbours );
		}

		iR  = interaction_NN(  distance( x, neighbours[1] ) );
		iRp = interaction_NN(  distance(  left(x), neighbours[1] ) );
	if(pos[neighbours[1]].state != 1)
		{
		flag=2;	
		interaction_error(type, x, neighbours ); 
		}
	}

   //--------------------------------------------------------------------
  }//----END OF 'if(partnum>1)'-------

return (iLp + iRp - iL -iR   + mux - muxp);

}
//*****************************************************************
double	GAdata::dEsR(int type, int x, int * neighbours )
{
double iL=0.0;	// old interaction energy with the left particle
double iR=0.0;	// old interaction energy with the right particle
double iLp=0.0;	// new interaction energy with the left particle
double iRp=0.0;	// new interaction energy with the right particle

double mux=0.0, muxp=0.0;




mux  = muN(x);
muxp = muN(right(x));


if ( type != 1)	//----*THIS* PARTICLE SHOULD ALWAYS BE A NUCL----
	{	
	flag=12;
	interaction_error(type, x, neighbours );
	}
//--------------------------------------------------------------------

if(Nucnum==1)
	{
	iLp = 0.0;
	iRp = 0.0;
	iL  = 0.0;
	iR  = 0.0;
	}
else 
   {
   if ( type == 1)	//----*THIS* PARTICLE IS A NUCL----
	{


	iL  = interaction_NN(  distance(neighbours[0], x ) );
	iLp = interaction_NN(  distance(neighbours[0], right(x) ) );

	iR  = interaction_NN(  distance(x, neighbours[1] ) );
	iRp = interaction_NN(  distance( right(x), neighbours[1] ) );

	}
	//--------------------------------------------------------------------

   }//----END OF 'if(partnum>1)'

return (iLp + iRp - iL -iR  + mux - muxp );

}

//***************************************************************************************
int react(int x,int Rtype, GAdata &P)
{
switch(Rtype)
	{

//---------------   DELETE THIS -----------------------------

	if ( ! (P.t >=0 && P.t <= P.tf ))
		{
		(*P.log) << "\n ERROR in react: P.t=" << P.t << endl;
		}

//---------------- DOWN TO HERE -----------------------------

	case 0: P.slide_Nuc_left( x); 	break;
	case 1: P.slide_Nuc_right( x); 	break;

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

if(Nucnum ==0)
	{
	for(i=0;i<Llim;i++)	
		{
		if(  pos[i].state ==1 || pos[i].part_left != i || pos[i].part_right != i )
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

return 1;
} //----end of function
//*****************************************************************
double GAdata::cleanup_occ( void)
{
int x=0,j=0;
double sum=0.0;

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
			pos[x].t_occ_N += (t - pos[x].t_bind_N); //--increment the occupation time.
			pos[x].t_bind_N = t;		// --- reset the binding time to "now"
			}
		}

	}

//---check that we have Nucnum particles.
if( j!= Nucnum )
	{
	flag=4013;
	process_error(1);
	exit(1);
	}

return sum;
}

//*****************************************************************
int GAdata::reset( void)
{
int x;
   Nucnum  = 0;
   a0=0.0;
   t=0.0;
   counter=0;

//------HERE WE SET THE INITIAL REACTION RATES----------------
for(x=0;x<Llim;x++)	
	{
	
	pos[x].state=0;		// initialize each lattice pos to 0 (nothing)
	pos[x].part_right = x;
	pos[x].part_left  = x;


	*pos[x].a_sNL = 0.0;
	*pos[x].a_sNR = 0.0;

	pos[x].t_occ_N  = 0.0;
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

//--------------------------------------------------------------------------------------------
for(x=0;x<Llim;x++)
	{
	current_t_occN =0.0;	

	if(pos[x].state==1)
		{
		N_there_now=1;
		current_t_occN += (t - pos[x].t_bind_N);
		}
	else
		N_there_now=0;


	current_t_occN  += pos[x].t_occ_N;

	fout << x << " \t " << (1.0/t)*current_t_occN << " \t " << N_there_now << " \t " <<  endl;
	}

fout.close();
*timestamps  << t << endl; //---print out the time of this plot


plotnum ++;

return 1;
}
//**********************************************************************
int GAdata::plot_snapshot_kymo( ofstream * fout)
{
int i,j,x;
char cpath[200];

int positions[Nucnum];


int N_there_now     = 0;
int N_found_already = 0;

for(j=0;j<200;j++)
	{
	cpath[j]='\0';
	}

//--------------------------------------------------------------------------------------------
for(x=0;x<Llim;x++)
	{

	if(pos[x].state==1)
		{
		positions[N_found_already]=x;
		N_found_already++;
		}

	}
if (N_found_already != Nucnum)
	{
	flag  = 17777;
	process_error(x);
	} 

int temp_wrap=wraparound_counter;

while(temp_wrap != 0)
	{
	if(temp_wrap <0 )
		{
		wrap_array(positions,Nucnum,1);
		temp_wrap++;
		}
	else if(temp_wrap >0)
		{
		wrap_array(positions,Nucnum,-1);
		temp_wrap--;
		}
	
	}

//--------------------------------------------------------------------------------------------
for(i=0;i<Nucnum;i++)
	{
	*fout << t << " \t " << positions[i] << endl;
	}
//--------------------------------------------------------------------------------------------

*fout << endl;

plotnum_kymo ++;

return 1;
}
//********************************************************************************************
int wrap_array( int * positions, const int size, const int direction)
{
int i=0;
int temp1;
int temp2;

if(direction == 1)
	{
	temp1 = positions[0];

	for(i=1;i<size;i++)
		{
		temp2        = positions[i];

		positions[i] = temp1;
		temp1        = temp2;
		}
	positions[0] = temp1;
	}
else if(direction == -1)
	{

	temp1 = positions[size-1];

	for(i=(size-2);i>=0;i--)
		{
		temp2        = positions[i];

		positions[i] = temp1;
		temp1        = temp2;
		}

	positions[size-1] = temp1;

	}
return 1;

}

#endif //--this ends the clause as to whether or not this symbol (i.e. this file) has been defined already.
