//---THESE ARE the absolute basic functions *ALL* (really all this time) IMPLEMENTATIONS OF THE GA SIMULATION
// ---last updated on  Fri Nov 29 10:52:02 CET 2013  by  Brendan.Osberg  at location  th-ws-e537

//  changes from  Fri Nov 29 10:52:02 CET 2013 : added remod_cand. function and int min_dist to GAdata; this way remodellers know when they've hit neighbours for the case of HNG.

//----next to last update accounted for the additional "boltzmann_on_removal" option in k_E
//----------------------------------------------------------------------------

#ifndef GA_ABSOLUTE_STANDARDS  //---check whether it's been defined already so we don't do it twice. 
#define GA_ABSOLUTE_STANDARDS

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
#include <gsl/gsl_blas.h>
#include <string.h>

//*****************************************************************
double GAdata::interaction_NN(int x)
{
if(x<=0)
	{
	flag=3001;
	process_error(x);
	}

if(x>NNRANGE)
	 return 0.0;
else
	return VNN[x-1];
}
//*****************************************************************
double GAdata::muN(int x )
{
return muNarray[x];
//return 0.0;
}
//*****************************************************************
double GAdata::interaction_NTF( int x)
{
if(x<=0)
	{
	flag=3002;
	process_error(x);
	}

if(x>NTFRANGE)
	 return 0.0;
else
	return VNTF[x-1];
}
//*****************************************************************
double GAdata::interaction_TFTF( int x)
{
if(x<m)
	{
	flag = 1002;
	process_error(x);
	}
else
	return	0.0;
}
//*****************************************************************
void GAdata::interaction_error(int type, int x, int *neighbours )
{

(*log) << "\n error! in interaction, type =" << type << ", x= " << x << endl;
(*log) << "\n flag = " << flag << endl;

(*log) << "neighbours[0] = " << neighbours[0] << endl;
(*log) << "neighbours[1] = " << neighbours[1] << endl;

(*log) << "\n\n unable to calculate interaction, exiting. \n\n ";
(*log).close();
exit(1);

}
//*****************************************************************
void GAdata::process_error( int x )
{

(*log) << "\n error! in some process, here x=" << x << endl;

(*log) << " pos[" << x << "] state = " << pos[x].state << endl;
(*log) << " recent reaction =" << recentRx << endl;
(*log) << " counter=" << counter << endl;
(*log) << " t=" << t  << endl;

(*log) << " Nucnum=" << Nucnum << endl;
(*log) << " flag = " << flag << endl;

(*log).close();
exit(1);

}

//*****************************************************************
double GAdata::k_E(const double dE)
{
//------ use just this line if you want to 
//------ have the "always evaluate energetics upon addition" option.

if ( boltzmann_on_uphill )
	{
	if (dE > 0.0)
		return gsl_sf_exp(-1.0*dE);
	else
 		return 1.0;
	}
else if (boltzmann_on_add || boltzmann_on_removal || boltzmann_on_addrem_intermed ) 
	{				
	//-- here we always return the energy in any of these circumstances, 
	//-- this function wouldn't get called unless it was appropriate to 
	//-- return the actual energy change. For the intermediate case, this 
	//-- function only gets fed the BZalpha fraction of the total interaction energy

	if(dE > 100)
		{//-----underflow catch.
		return 0.0;
		}
	else	
		{
		return gsl_sf_exp(-1.0*dE);
		}
	}
else
	{
	cout << "\n ERROR: not sure when to take the boltzmann factor.\n exiting.\n\n";
	exit(1); 
	}

}

//*****************************************************************
int GAdata::left(const int x)
{
int result;

if( x-1 == -1)
	result = Llim-1;
else if( (x-1 >= 0) && (x-1 < Llim-1))
	result = x-1;
else
	{
	cout << "\n error! incorrect bounds in left function!\n";
	exit(1);
	}

return result;
}
//*****************************************************************
int GAdata::right(const int x)
{
int result;

if( x == Llim-1)
	result = 0;
else if( (x+1 >= 1 ) && (x+1 < Llim))
	result = x+1;
else
	{
	cout << "\n error! incorrect bounds in right function !\n";
	exit(1);
	}

return result;
}
//*****************************************************************
int GAdata::distance(const int a, const int b)	// ORDER MATTERS!!!
{	//----b is assumed to be always to the right of a (with wrap-around)
int temp = b-a;

if( temp > 0)
	return temp;
else
	return (Llim + temp );	//result is ALWAYS > 0, hence interation with SELF is zer
}
//*****************************************************************
int GAdata::remod_cand(const int a, const int b)	
{
//--- ASSUMING THERE IS A NUCL TO THE LEFT AT a, AND TO THE RIGHT AT b, 
//--- IS THIS PAIR A CANDIDATE FOR BEING REMODELLED CLOSER TOGETHER
//--- MUST BE WITHIN RANGE, BUT FURTHER AWAY THAN min_dist. 
//----b is assumed to be always to the right of a (with wrap-around)
int result;

int d = distance(a,b);

if ( d <= RMrange  &&  d > min_dist  )
	{	//---within range, but further away than min_dist.
	result = 1;
	}
else
	{	//---either of the above two conditions not true.
	result = 0;
	}

return result;
}
//********************************************************************************
int GAdata::printout_states( ofstream *fout)
{
int n = 0;
int x  = Llim-1;
int xp = pos[x].part_right;

if(Nucnum == 0)
	return 0;
else
	{
	while(1)
		{
		  if(pos[xp].state == 1)
			{
			*fout << t << " \t " << xp  << " \t " << "-10" << endl;
			n++;
			}
		  else if(pos[xp].state == 2)
			{
			*fout << t  << " \t " << "-10" << " \t " << xp << endl;
			n++;
			}
		  else
			{
			flag=115;
			process_error(-100);
			}
		  x=xp;
		  xp = pos[x].part_right;

		  if(  (xp < x) || (Nucnum==1))
			break;
		}
	}

return n;
}
//********************************************************************************
int GAdata::printout_avgs( ofstream *fout)
{
int x=0;

for(x=0;x<Llim;x++)
	{
	*fout << x << " \t " << pos[x].t_occ_N << endl;
	}

return 1;
}
//********************************************************************************
int GAdata::reset_prevnext_pointers()
{
int x;
for(x=0;x<Llim;x++)
	{
	pos[x].part_left  = x;
	pos[x].part_right = x;
	}

return 1;
}
//*****************************************************************
/*
bool GAdata::printtime(void)
{
int i=0,j=0;
bool result=false;
double temp=0.0;

if(t < short_print_timescale)
  {
  temp = (double(plotnum)/(double(Nplots2makeshort)))*short_print_timescale; // see what fraction 
								//of tf we've 	completed so far

  if(t > temp)
	{
	result	=true;
	}
  else	
	result = false;
  }
else
  {
  temp = (double(plotnum-Nplots2makeshort)/(double(Nplots2makelong)))*(tf-short_print_timescale); // see what fraction of tf we've 	
  if((t-short_print_timescale) > temp)
	{
	result	=true;
	}
  else	
	result = false;
  }
return result;
}
*/

//*****************************************************************
bool GAdata::printtime_kymo(void)
{
int i=0,j=0;
bool result=false;
double temp=0.0;

temp = (double(plotnum_kymo)/(double(Nplots2make_kymo) ) )  *(tf); // see what fraction of tf we've 	

  if(t > temp)
	{
	result	=true;
	}
  else	
	result = false;
return result;
}
//**********************************************************************

double GAdata::interaction_dEadd(int type, int x, int *neighbours )
{
double iL=0.0;	// interaction energy with the left particle
double iR=0.0;	// interaction energy with the right particle
double iO=0.0;	// "Outer" interaction energy between the left and right particle
		// that would exist if this intervening particle weren't there.

if( partnum  == 0)
	return 0;

int locL = pos[x].part_left;
int locR = pos[x].part_right; 

// neighbours[0]= type left
// neighbours[1]= location left
// neighbours[2]= type right
// neighbours[3]= location right

	//--------------------------------------------------------------------
if ( type == 1)	//----*THIS* PARTICLE IS A NUCL----
	{
	if(neighbours[0] ==1)	//left neighbour is a nucl
		iL = interaction_NN(  distance(neighbours[1], x ) );
	else if(neighbours[0] ==2) // left neighbour is a TF
		iL = interaction_NTF( (distance(neighbours[1], x )-(m-1)) );
	else
		{	
		flag=13;
		interaction_error(type, x, neighbours );
		}

	if(neighbours[2] ==1)	// right neighbour is a nucl
		iR = interaction_NN(  distance(x, neighbours[3] ) );
	else if(neighbours[2] ==2) 
		iR = interaction_NTF(  distance(x, neighbours[3] ) );	else
		{	flag=14;	interaction_error(type, x, neighbours );}
	}
	//--------------------------------------------------------------------
else if ( type == 2) //----*THIS* PARTICLE IS A TF----
	{
	if(neighbours[0] ==1)
		iL = interaction_NTF(  distance(neighbours[1], x ) );
	else if(neighbours[0] ==2)
		{//-----error catch-----
		if(distance(neighbours[1], x ) < m )
			{
			(*log) << "\n error! TF overlap!\n\n";
			flag=15;	
			interaction_error(type, x, neighbours );
			}
		else
			iL=0.0;
		}
	else
		{	flag=16;	interaction_error(type, x, neighbours );}

	if(neighbours[2] ==1)
		iR = interaction_NTF( (distance(x,neighbours[3] ) - (m -1)  ) );
	else if(neighbours[2] ==2)
		{//-----error catch-----
		if(distance(x,neighbours[3] ) < m )
			{
			(*log) << "\n error! TF overlap!\n\n";
			flag=17;	
			interaction_error(type, x, neighbours );
			}
		else
			iR=0.0;
		}

	}
	//--------------------------------------------------------------------
else
	{	
	flag=17;
	interaction_error(type, x, neighbours );
	}


if(neighbours[1] == neighbours[3]) // left and right neighbours are the same particle:
	iO=0.0;			   // it can't interact with itself
else 
	{
	if ( (neighbours[0] ==1) &&  (neighbours[2] ==1) ) 
		iO= interaction_NN(  distance(neighbours[1], neighbours[3] ) );
	else if ( (neighbours[0] ==1) &&  (neighbours[2] ==2) ) 
		iO= interaction_NTF( distance(neighbours[1], neighbours[3] ) );
	else if ( (neighbours[0] ==2) &&  (neighbours[2] ==1) ) 
		iO= interaction_NTF( distance(neighbours[1], neighbours[3] ) - (m-1) );
	else if ( (neighbours[0] ==2) &&  (neighbours[2] ==2) ) 
		iO= interaction_TFTF( distance(neighbours[1], neighbours[3] )  );


	else
		{	
		flag=18;	interaction_error(type, x, neighbours );}
		}


	return iL + iR -iO; //----- ALWAYS returns the change in energy produced by addition.
			    //----- selection of uphill/downhill addition/removal criteria is left to higher functions.
}


//***************************************************************************************
int react(int x,int Rtype, GAdata &P)
{
switch(Rtype) ///------ the version in this file is base-4
	{
	case 0: P.remove_Nuc( x ); 	break;
	case 1: P.add_Nuc(x); 		break;
	case 2: P.slide_Nuc_left( x); 	break;
	case 3: P.slide_Nuc_right( x); 	break;

	default:

	P.flag=102;
	P.recentRx=-10;
	P.process_error(x);	

	}
}

//*****************************************************************
double	GAdata::dEsL(int type, int x, int * neighbours )
{
double iL=0.0;	// old interaction energy with the left particle
double iR=0.0;	// old interaction energy with the right particle
double iLp=0.0;	// new interaction energy with the left particle
double iRp=0.0;	// new interaction energy with the right particle

double mux=0.0, muxp=0.0;

if ( type == 1)	//----*THIS* PARTICLE IS A NUCL----
	{
	mux  = muN(x);
	muxp = muN(left(x));
	}
else if ( type == 2) //----*THIS* PARTICLE IS A TF----
	{
	cout << "\n ERROR: TFs should never slide\n";
	exit(1);
	/*
	mux  = muTF(x);
	muxp = muTF(left(x));
	*/
	}
else		     //-----DON'T KNOW WHAT THIS IS.
	{	
	flag=6;
	interaction_error(type, x, neighbours  );
	}
	//--------------------------------------------------------------------
if(partnum==1)
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

	if(neighbours[0] ==1)	//left neighbour is a nucl
		{
		iL  = interaction_NN( distance(neighbours[1], x ) );
		iLp = interaction_NN( distance(neighbours[1], left(x) )  );
		}
	else if(neighbours[0] ==2) // left neighbour is a TF
		{
		iL  = interaction_NTF(  distance(neighbours[1], x )- (m-1)  );
		iLp = interaction_NTF(  distance(neighbours[1], left(x) )-(m-1)  );
		}
	else
		{	
		flag=1;	
		interaction_error(type, x, neighbours );
		}

	if(neighbours[2] ==1)	// right neighbour is a nucl
		{
		iR  = interaction_NN(  distance( x, neighbours[3] ) );
		iRp = interaction_NN(  distance(  left(x), neighbours[3] ) );
		}
	else if(neighbours[2] ==2) // right neighbour is a TF
		{
		iR  = interaction_NTF(  distance( x, neighbours[3]) );	
		iRp = interaction_NTF(  distance( left(x), neighbours[3] ) );	
		}
	else
		{
		flag=2;	
		interaction_error(type, x, neighbours ); 
		}
	}
	//--------------------------------------------------------------------
else if ( type == 2) //----*THIS* PARTICLE IS A TF----
	{

	if(neighbours[0] ==1)
		{
		iL  = interaction_NTF( distance(neighbours[1], x ) );
		iLp = interaction_NTF( distance(neighbours[1], left(x) ) );
		}
	else if(neighbours[0] ==2)
		{//-----error catch-----
		if( ( distance(neighbours[1], x ) < m) || ( distance(neighbours[1], left(x) ) < m ) )
			{	flag=3;		interaction_error(type, x, neighbours );}
		else
			{
			iL=0.0;
			iLp=0.0;
			}
		}
	else
		{	
		flag=4;	
		interaction_error(type, x, neighbours );
		}

	if(neighbours[2] ==1)
		{
		iR  = interaction_NTF(  (distance( x,neighbours[1] ) - (m-1) ) );
		iRp = interaction_NTF(  (distance( left(x),neighbours[1] ) - (m-1) ) );
		}
	else if(neighbours[2] ==1)
		{//-----error catch-----
		if( distance(x,neighbours[3]) < m )
			{
			flag=5;
			interaction_error(type, x, neighbours ); 
			}
		else
			iR=0.0;
		}

	}
	//--------------------------------------------------------------------
  }//----END OF 'if(partnum>1)'-------

return (iLp + iRp - iL -iR   + mux - muxp); //-------ALWAYS RETURNS THE CHANGE IN ENERGY 
					    //---- SELECTION AGAINST 0 (UPHILL/DOWNHILL) IS DONE BY HIGHER-ORDER FUNCTIONS.

}
//*****************************************************************
double	GAdata::dEsR(int type, int x, int * neighbours )
{
double iL=0.0;	// old interaction energy with the left particle
double iR=0.0;	// old interaction energy with the right particle
double iLp=0.0;	// new interaction energy with the left particle
double iRp=0.0;	// new interaction energy with the right particle

double mux=0.0, muxp=0.0;



if ( type == 1)	//----*THIS* PARTICLE IS A NUCL----
	{
	mux  = muN(x);
	muxp = muN(right(x));
	}
else if ( type == 2) //----*THIS* PARTICLE IS A TF----
	{

	cout  << "\n ERROR: TFs should never slide.\n";
	exit(1);
	/*
	mux  = muTF(x);
	muxp = muTF(right(x));
	*/
	}
else		     //-----DON'T KNOW WHAT THIS IS
	{	
	flag=12;
	interaction_error(type, x, neighbours );
	}
//--------------------------------------------------------------------

if(partnum==1)
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

	if(neighbours[0] ==1)	//left neighbour is a nucl
		{
		iL  = interaction_NN(  distance(neighbours[1], x ) );
		iLp = interaction_NN(  distance(neighbours[1], right(x) ) );
		}
	else if(neighbours[0] ==2) // left neighbour is a TF
		{
		iL  = interaction_NTF(  ( distance(neighbours[1], x ) - (m-1) ) );
		iLp = interaction_NTF(  ( distance(neighbours[1], right(x) )-  (m-1) ) );
		}
	else
		{	
		flag=7;	interaction_error(type, x, neighbours );
		}

   if(neighbours[2] ==1)	// right neighbour is a nucl
		{
		iR  = interaction_NN(  distance(x, neighbours[3] ) );
		iRp = interaction_NN(  distance( right(x), neighbours[3] ) );
		}
	else if(neighbours[2] ==2) 
		{
		iR  = interaction_NTF(  distance(x, neighbours[3] ) );	
		iRp = interaction_NTF(  distance(right(x), neighbours[3] ) );	
		}
	else
		{	flag=8;	interaction_error(type, x, neighbours );}
	}
	//--------------------------------------------------------------------
   else if ( type == 2) //----*THIS* PARTICLE IS A TF----
	{

	if( distance(x,neighbours[3]) < m )//--error flag
		{	//---things shouldn't be closer than 'm' to the right.
		flag=11;
		interaction_error(type, x, neighbours );
		}//---this contingency will hopefull always be ignored.

	if(neighbours[0] ==1)
		{
		iL  = interaction_NTF( distance(neighbours[1], x ) );
		iLp = interaction_NTF( distance(neighbours[1], right(x) ) );
		}
	else if(neighbours[0] ==2)
		{//-----error catch-----
		if(  distance(x ,neighbours[3] ) <= m )
		   {	
		   flag=9;
		   interaction_error(type, x, neighbours );
		   }
		else
			{
			iL=0.0;
			iLp=0.0;
			}
		}
	else
		{
		flag=10;
		interaction_error(type, x, neighbours  );
		}

	if(neighbours[2] ==1)
		{
		iR  = interaction_NTF(  ( distance(x,neighbours[3]) - (m-1)) );
		iRp = interaction_NTF(  ( distance( right(x),neighbours[3]) - (m-1) ) );
		}
	else if(neighbours[2] ==2)
		{//-----error catch-----
			iR=0.0;
		}

	}
	//--------------------------------------------------------------------
   }//----END OF 'if(partnum>1)'

return (iLp + iRp - iL -iR  + mux - muxp ); //-------ALWAYS RETURNS THE CHANGE IN ENERGY 
					    //---- SELECTION AGAINST 0 (UPHILL/DOWNHILL) IS DONE BY HIGHER-ORDER FUNCTIONS.

}

//*****************************************************************
int GAdata::reset( void)
{
int x;
   Nucnum  = 0;
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


	*pos[x].a_addN  = ka_N*gsl_sf_exp(muN(x)); 	//-- ADDING NUCLEOSOMES  -------------


	a0+=*pos[x].a_addN;

	//-----------------AND SET AVERAGE OCCUPANCY TO ZERO INITIALLY------------

	pos[x].t_occ_N  = 0.0;

	pos[x].t_bind_N  = -1.0;
	} //---end for x=0..Llim	

}

//**********************************************************************
int GAdata::plot_snapshot_kymo( ofstream * foutN )
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

	}

if(foundN)
	{
	*foutN << endl;
	}


plotnum_kymo ++;

return 1;
}
//***************************************************************************
int  assign_transfermatrix_elements(gsl_matrix *T, const string NGtype, const int a, const double mu, const double * VNN)
{//---- "a" is the footprint size
int i,j;
int result=0;

if(T->size1 != T->size2)
	{
	cout << "\nERROR: size1/2 not equal in T matrix where we're setting elements. Exiting.\n";
	exit(1);
	}

if(T->size1 != a+1)
	{
	cout << "\nERROR: size1 not equal to a+1 in T matrix where we're setting elements. Exiting.\n";
	exit(1);
	}

gsl_matrix_set (T, 0  , 0, 1.0);
gsl_matrix_set (T, a  , 0, 1.0);

gsl_matrix_set (T, 0  , 1, gsl_sf_exp(mu) );
gsl_matrix_set (T, a  , 1, gsl_sf_exp(mu) );

for(i=1; i<(a); i++)
	{
	gsl_matrix_set(T,i,i+1,1.0);
	}

if( strcmp(NGtype.c_str(),"SNG")==0 ||  strcmp(NGtype.c_str(),"LNG")==0  )
        {
 
	for(i=1; i<(a); i++)
		{
		gsl_matrix_set(T,i,1,gsl_sf_exp(mu-VNN[i-1]));
		}
       result=1;
        }
else if(  strcmp(NGtype.c_str(),"HNG")!=0 )
	{
	cout << "\nERROR: NGtype not matched in transfer matrix assignment. Exiting \n";
	exit(1);
	}

return result;
}
//*****************************************************************************
int bren_matrix_pow(const gsl_matrix *T, const int n, gsl_matrix * output)
{

int i;

if(T->size1 != T->size2)
	{
	cout << "\n ERROR:matrix sent to bren_matrix_pow is not square. Exiting.\n";
	exit(1);
	}
int size=T->size1;

gsl_matrix * temp  = gsl_matrix_alloc (size, size);
gsl_matrix_set_zero(temp);


gsl_matrix_set_identity (output);

for(i=0;i<n;i++)
	{

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, T, output, 0.0, temp);
	
	gsl_matrix_memcpy ( output, temp);	
	
	}

return 1;
}
//*******************************************************************************
double bren_get_matrix_trace(const gsl_matrix * T)
{
int i;
double result=0.0;

if(T->size1 != T->size2)
        {
        cout << "\n ERROR:matrix sent to bren_matrix_pow is not square. Exiting.\n";
        exit(1);
        }
int size=T->size1;

for (i=0;i<size;i++)
	{
	result += gsl_matrix_get(T,i,i);
	}

return result;

}



#endif //--this ends the clause as to whether or not this symbol (i.e. this file) has been defined already.
