/*Copyright by Binbin Lu*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>

void nodeExisted(double *nodeXlist,double *nodeYlist,int *n,double *X,double *Y, int *tag);
void edgelength(double *nodeXlist,double *nodeYlist,int *n, double *edgelength, int *longlat);
void gc_el(double *lon1, double *lon2, double *lat1, double *lat2, double *gel);
void dist_p2e(double *nodeXlist,double *nodeYlist,int *n,int *idx,double *X,double *Y, double *dist, double *fx, double *fy);
void footxy(double *nodeXlist,double *nodeYlist,double *X,double *Y, double *dist, double *fx, double *fy);
int maxindex(double *list, int n);
int minindex(double *list, int n);
void addDegree(int *NIDs, int *n,int *nid, int *DL);
void minusDegree(int *NIDs, int *n,int *nid, int *DL);

void nodeExisted(double *nodeXlist,double *nodeYlist,int *n,double *X,double *Y, int *tag)
{
	int N=*n, i;
	double nx,ny;
	for (i=0; i<N; i++)
	   {
		 nx=nodeXlist[i];
		 if(nx==X[0])
		 {
			 ny=nodeYlist[i];
			 if (ny==Y[0])
			    {tag[0]=i+1;
			     break;}
			 }
		 else
		    {continue;}
		   }
	}

void edgelength(double *nodeXlist,double *nodeYlist,int *n, double *edgelength, int *longlat)
{
	int N=*n, i;
	double el[1],gel[1];
	el[0]=(double)0;
	if (longlat[0]==0)
	{
		for(i=0; i<N-1; i++)
		{
		    el[0]=el[0]+hypot((nodeXlist[i+1]-nodeXlist[i]),(nodeYlist[i+1]-nodeYlist[i]));
			} 
		}
	else
	{
		for(i=0; i<N-1; i++)
		{
			gc_el(nodeXlist+i+1,nodeXlist+i,nodeYlist+i+1,nodeYlist+i+1, gel);
		    el[0]=el[0]+gel[0];
			} 
		}
		edgelength[0]=el[0];			
	}

void gc_el(double *lon1, double *lon2, double *lat1, double *lat2, double *gel)
{
	double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C;
    double w, R, a, f, D, H1, H2;
    double lat1R, lat2R, lon1R, lon2R, DE2RA;
    
    DE2RA = M_PI/180;
    a = 6378.137;              /* WGS-84 equatorial radius in km */
    f = 1.0/298.257223563;     /* WGS-84 ellipsoid flattening factor */
    
    lat1R = lat1[0]*DE2RA;
    lat2R = lat2[0]*DE2RA;
    lon1R = lon1[0]*DE2RA;
    lon2R = lon2[0]*DE2RA;
    
    F = ( lat1R + lat2R )/2.0;
    G = ( lat1R - lat2R )/2.0;
    L = ( lon1R - lon2R )/2.0;

    sinG2 = R_pow_di( sin( G ), 2 );
    cosG2 = R_pow_di( cos( G ), 2 );
    sinF2 = R_pow_di( sin( F ), 2 );
    cosF2 = R_pow_di( cos( F ), 2 );
    sinL2 = R_pow_di( sin( L ), 2 );
    cosL2 = R_pow_di( cos( L ), 2 );

    S = sinG2*cosL2 + cosF2*sinL2;
    C = cosG2*cosL2 + sinF2*sinL2;

    w = atan( sqrt( S/C ) );
    R = sqrt( S*C )/w;

    D = 2*w*a;
    H1 = ( 3*R - 1 )/( 2*C );
    H2 = ( 3*R + 2 )/( 2*S );

    gel[0]= D*( 1 + f*H1*sinF2*cosG2 - f*H2*cosF2*sinG2 );
	}
	
void dist_p2e(double *nodeXlist,double *nodeYlist,int *n,int *idx,double *X,double *Y, double *dist, double *fx, double *fy)
{
	int N=*n;
	double d=100000000000000;
	double FX=0,FY=0;
	for (int i=0; i<N-1;i++)
	{
		footxy(nodeXlist+i,nodeYlist+i,X,Y, dist, fx, fy);
		if (dist[0]<d)
		{
			d=dist[0];
			FX=fx[0];
			FY=fy[0];
			idx[0]=i+1;
			}
		}
		dist[0]=d;
		fx[0]=FX;
		fy[0]=FY; 
	}
	/*footxy is to compute the coordinate of the foot point from a given point to a line (defined by two points)*/
	/*The lengths of nodeXlist and nodeYlist are all supposed to be 2*/
void footxy(double *nodeXlist,double *nodeYlist, double *X,double *Y, double *dist, double *fx, double *fy)
{
	double A,B,C;
	int Xmax, Xmin, Ymax, Ymin;
	A=nodeYlist[1]-nodeYlist[0];
	B=nodeXlist[0]-nodeXlist[1];
	C=nodeYlist[0]*nodeXlist[1]-nodeYlist[1]*nodeXlist[0];
	dist[0]=fabs(A*X[0]+B*Y[0]+C)/sqrt(A*A+B*B);
	if (A==(double)0)
	{
		if (B==(double)0)
		{
			fx[0]=nodeXlist[0];
			fy[0]=nodeYlist[0];
			}
		else
		{
			fy[0]=nodeYlist[0];
			fx[0]=X[0];
			}
		}
	else
	{
		if (B==(double)0)
		{
			fx[0]=nodeXlist[0];
			fy[0]=Y[0];
			}
		else
		{
			fx[0]=B*B*X[0]/(A*A+B*B)-A*B*Y[0]/(A*A+B*B)-A*C/(A*A+B*B);
			fy[0]=-A*B*X[0]/(A*A+B*B)+A*A*Y[0]/(A*A+B*B)-B*C/(A*A+B*B);
			}
		}
	
	if (A==(double)0)
	{
		Xmax=maxindex(nodeXlist, (int)2);
		Xmin=minindex(nodeXlist, (int)2);
		if (fx[0]<nodeXlist[Xmin])
		{
			dist[0]=hypot((nodeXlist[Xmin]-X[0]),(nodeYlist[Xmin]-Y[0]));
			fx[0]=nodeXlist[Xmin];
			}
		if (fx[0]>nodeXlist[Xmax])
		{
			dist[0]=hypot((nodeXlist[Xmax]-X[0]),(nodeYlist[Xmax]-Y[0]));
			fx[0]=nodeXlist[Xmin];
			}
		}
	else
	{
		Ymax=maxindex(nodeYlist, (int)2);
		Ymin=minindex(nodeYlist, (int)2);
		if (fy[0]<nodeYlist[Ymin])
		{
			dist[0]=hypot((nodeXlist[Ymin]-X[0]),(nodeYlist[Ymin]-Y[0]));
			fx[0]=nodeXlist[Ymin];
			fy[0]=nodeYlist[Ymin];
			}
		if (fy[0]>nodeYlist[Ymax])
		{
			dist[0]=hypot((nodeXlist[Ymax]-X[0]),(nodeYlist[Ymax]-Y[0]));
			fx[0]=nodeXlist[Ymax];
			fy[0]=nodeYlist[Ymax];
			}
		}
	}

int maxindex(double *list, int n)
{
	int tag=0, N=n;
	double e=-100000000000000;
	for (int i=0;i<N;i++)
	{
		if(e<list[i])
		{
			tag=i;
			e=list[i];
			}
		
		}
		return(tag);
	}

int minindex(double *list, int n)
{
	int tag=0, N=n;
	double e=100000000000000;
	for (int i=0;i<N;i++)
	{
		if(e>list[i])
		{
			tag=i;
			e=list[i];
			}
		
		}
		return(tag);
	}

void addDegree(int *NIDs, int *n,int *nid, int *DL)
{
	int N=*n;
	for (int i=0; i<N; i++)
	{
		if (NIDs[i]==nid[0])
		{
			DL[i]=DL[i]+1;
			break;
			}
		}
	}
void minusDegree(int *NIDs, int *n,int *nid, int *DL)
{
	int N=*n;
	for (int i=0; i<N; i++)
	{
		if (NIDs[i]==nid[0])
		{
			DL[i]=DL[i]-1;
			break;
			}
		}
	}
