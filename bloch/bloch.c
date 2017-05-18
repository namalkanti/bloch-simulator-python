#include <stdio.h>
#include <math.h>

#define GAMMA   26753.0
#define TWOPI	6.283185

#define DEBUG



void multmatvec(double *mat, double *vec, double *matvec)

	/* Multiply 3x3 matrix by 3x1 vector. */

{
*matvec++ = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2];
*matvec++ = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2];
*matvec++ = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2];
}



void addvecs(double *vec1, double *vec2, double *vecsum)

	/* Add two 3x1 Vectors */

{
*vecsum++ = *vec1++ + *vec2++;
*vecsum++ = *vec1++ + *vec2++;
*vecsum++ = *vec1++ + *vec2++;
}




void adjmat(double *mat, double *adj)

/* ======== Adjoint of a 3x3 matrix ========= */

{
*adj++ = (mat[4]*mat[8]-mat[7]*mat[5]);
*adj++ =-(mat[1]*mat[8]-mat[7]*mat[2]);
*adj++ = (mat[1]*mat[5]-mat[4]*mat[2]);
*adj++ =-(mat[3]*mat[8]-mat[6]*mat[5]);
*adj++ = (mat[0]*mat[8]-mat[6]*mat[2]);
*adj++ =-(mat[0]*mat[5]-mat[3]*mat[2]);
*adj++ = (mat[3]*mat[7]-mat[6]*mat[4]);
*adj++ =-(mat[0]*mat[7]-mat[6]*mat[1]);
*adj++ = (mat[0]*mat[4]-mat[3]*mat[1]);
}


void zeromat(double *mat)

/* ====== Set a 3x3 matrix to all zeros	======= */

{
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
}


void eyemat(double *mat)

/* ======== Return 3x3 Identity Matrix  ========= */

{
zeromat(mat);
mat[0]=1;
mat[4]=1;
mat[8]=1;

}

double detmat(double *mat)

/* ======== Determinant of a 3x3 matrix ======== */

{
double det;

det = mat[0]*mat[4]*mat[8];
det+= mat[3]*mat[7]*mat[2];
det+= mat[6]*mat[1]*mat[5];
det-= mat[0]*mat[7]*mat[5];
det-= mat[3]*mat[1]*mat[8];
det-= mat[6]*mat[4]*mat[2];

return det;
}


void scalemat(double *mat, double scalar)

/* ======== multiply a matrix by a scalar ========= */

{
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
}


void invmat(double *mat, double *imat)

/* ======== Inverse of a 3x3 matrix ========= */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */

{
double det;
int count;

det = detmat(mat);	/* Determinant */
adjmat(mat, imat);	/* Adjoint */

for (count=0; count<9; count++)
	*imat++ /= det;
}


void addmats(double *mat1, double *mat2, double *matsum)

/* ====== Add two 3x3 matrices.	====== */

{
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
}


void multmats(double *mat1, double *mat2, double *matproduct)

/* ======= Multiply two 3x3 matrices. ====== */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */

{
*matproduct++ = mat1[0]*mat2[0] + mat1[3]*mat2[1] + mat1[6]*mat2[2];
*matproduct++ = mat1[1]*mat2[0] + mat1[4]*mat2[1] + mat1[7]*mat2[2];
*matproduct++ = mat1[2]*mat2[0] + mat1[5]*mat2[1] + mat1[8]*mat2[2];
*matproduct++ = mat1[0]*mat2[3] + mat1[3]*mat2[4] + mat1[6]*mat2[5];
*matproduct++ = mat1[1]*mat2[3] + mat1[4]*mat2[4] + mat1[7]*mat2[5];
*matproduct++ = mat1[2]*mat2[3] + mat1[5]*mat2[4] + mat1[8]*mat2[5];
*matproduct++ = mat1[0]*mat2[6] + mat1[3]*mat2[7] + mat1[6]*mat2[8];
*matproduct++ = mat1[1]*mat2[6] + mat1[4]*mat2[7] + mat1[7]*mat2[8];
*matproduct++ = mat1[2]*mat2[6] + mat1[5]*mat2[7] + mat1[8]*mat2[8];
}


void calcrotmat(double nx, double ny, double nz, double *rmat)

	/* Find the rotation matrix that rotates |n| radians about
		the vector given by nx,ny,nz				*/
{
double ar, ai, br, bi, hp, cp, sp;
double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
double phi;

phi = sqrt(nx*nx+ny*ny+nz*nz);

if (phi == 0.0)
	{
	*rmat++ = 1;
	*rmat++	= 0;
	*rmat++ = 0;
	*rmat++ = 0;
	*rmat++ = 1;
	*rmat++	= 0;
	*rmat++ = 0;
	*rmat++ = 0;
	*rmat++ = 1;
	}

/*printf("calcrotmat(%6.3f,%6.3f,%6.3f) -> phi = %6.3f\n",nx,ny,nz,phi);*/

else
	{
	/* First define Cayley-Klein parameters 	*/
	hp = phi/2;
	cp = cos(hp);
	sp = sin(hp)/phi;	/* /phi because n is unit length in defs. */
	ar = cp;
	ai = -nz*sp;
	br = ny*sp;
	bi = -nx*sp;

 	/* Make auxiliary variables to speed this up	*/

	arar = ar*ar;
	aiai = ai*ai;
	arai2 = 2*ar*ai;
	brbr = br*br;
	bibi = bi*bi;
	brbi2 = 2*br*bi;
	arbi2 = 2*ar*bi;
	aibr2 = 2*ai*br;
	arbr2 = 2*ar*br;
	aibi2 = 2*ai*bi;


	/* Make rotation matrix.	*/

	*rmat++ = arar-aiai-brbr+bibi;
	*rmat++ = -arai2-brbi2;
	*rmat++ = -arbr2+aibi2;
	*rmat++ =  arai2-brbi2;
	*rmat++ = arar-aiai+brbr-bibi;
	*rmat++ = -aibr2-arbi2;
	*rmat++ =  arbr2+aibi2;
	*rmat++ =  arbi2-aibr2;
	*rmat++ = arar+aiai-brbr-bibi;
	}
}



void zerovec(double *vec)

/*	Set a 3x1 vector to all zeros	*/

{
*vec++=0;
*vec++=0;
*vec++=0;
}


int times2intervals( double *endtimes, double *intervals, long n)
/* ------------------------------------------------------------
	Function takes the given endtimes of intervals, and
	returns the interval lengths in an array, assuming that
	the first interval starts at 0.

	If the intervals are all greater than 0, then this
	returns 1, otherwise it returns 0.
   ------------------------------------------------------------ */

{
int allpos;
int count;
double lasttime;

allpos=1;
lasttime = 0.0;

for (count = 0; count < n; count++)
	{
	intervals[count] = endtimes[count]-lasttime;
	lasttime = endtimes[count];
	if (intervals[count] <= 0)
		allpos =0;
	}

return (allpos);
}






void blochsim(double *b1real, double *b1imag,
		double *xgrad, double *ygrad, double *zgrad, double *tsteps,
		int ntime, double *e1, double *e2, double df,
		double dx, double dy, double dz,
		double *mx, double *my, double *mz, int mode)

	/* Go through time for one df and one dx,dy,dz.		*/

{
int tcount;
double gammadx;
double gammady;
double gammadz;
double rotmat[9];
double amat[9], bvec[3];	/* A and B propagation matrix and vector */
double arot[9], brot[3];	/* A and B after rotation step. */
double decmat[9];		/* Decay matrix for each time step. */
double decvec[3];		/* Recovery vector for each time step. */
double rotx,roty,rotz;		/* Rotation axis coordinates. */
double imat[9], mvec[3];
double mcurr0[3];		/* Current magnetization before rotation. */
double mcurr1[3];		/* Current magnetization before decay. */

eyemat(amat); 		/* A is the identity matrix.	*/
eyemat(imat); 		/* I is the identity matrix.	*/

zerovec(bvec);
zerovec(decvec);
zeromat(decmat);

gammadx = dx*GAMMA;	/* Convert to Hz/cm */
gammady = dy*GAMMA;	/* Convert to Hz/cm */
gammadz = dz*GAMMA;	/* Convert to Hz/cm */


mcurr0[0] = *mx;		/* Set starting x magnetization */
mcurr0[1] = *my;		/* Set starting y magnetization */
mcurr0[2] = *mz;		/* Set starting z magnetization */


for (tcount = 0; tcount < ntime; tcount++)
	{
		/*	Rotation 	*/

	rotz = -(*xgrad++ * gammadx + *ygrad++ * gammady + *zgrad++ * gammadz +
								df*TWOPI ) * *tsteps;
	rotx = (- *b1real++ * GAMMA * *tsteps);
	roty = (+ *b1imag++ * GAMMA * *tsteps++);
	calcrotmat(rotx, roty, rotz, rotmat);

	if (mode == 1)
		{
		multmats(rotmat,amat,arot);
		multmatvec(rotmat,bvec,brot);
		}
	else
		multmatvec(rotmat,mcurr0,mcurr1);


		/* 	Decay	*/

	decvec[2]= 1- *e1;
	decmat[0]= *e2;
	decmat[4]= *e2++;
	decmat[8]= *e1++;

	if (mode == 1)
		{
		multmats(decmat,arot,amat);
		multmatvec(decmat,brot,bvec);
		addvecs(bvec,decvec,bvec);
		}
	else
		{
		multmatvec(decmat,mcurr1,mcurr0);
		addvecs(mcurr0,decvec,mcurr0);
		}

	/*
	printf("rotmat = [%6.3f  %6.3f  %6.3f ] \n",rotmat[0],rotmat[3],
	  			rotmat[6]);
	printf("         [%6.3f  %6.3f  %6.3f ] \n",rotmat[1],rotmat[4],
				rotmat[7]);
	printf("         [%6.3f  %6.3f  %6.3f ] \n",rotmat[2],rotmat[5],
				rotmat[8]);
	printf("A = [%6.3f  %6.3f  %6.3f ] \n",amat[0],amat[3],amat[6]);
	printf("    [%6.3f  %6.3f  %6.3f ] \n",amat[1],amat[4],amat[7]);
	printf("    [%6.3f  %6.3f  %6.3f ] \n",amat[2],amat[5],amat[8]);
	printf(" B = <%6.3f,%6.3f,%6.3f> \n",bvec[0],bvec[1],bvec[2]);
	printf("<mx,my,mz> = <%6.3f,%6.3f,%6.3f> \n",
		amat[6] + bvec[0], amat[7] + bvec[1], amat[8] + bvec[2]);

	printf("\n");
	*/

	if (mode == 2)		/* Sample output at times.  */
					/* Only do this if transient! */
		{
		*mx = mcurr0[0];
		*my = mcurr0[1];
		*mz = mcurr0[2];

		mx++;
		my++;
		mz++;
		}
	}



	/* If only recording the endpoint, either store the last
		point, or calculate the steady-state endpoint. */

if (mode==0)		/* Indicates start at given m, or m0. */
	{
	*mx = mcurr0[0];
	*my = mcurr0[1];
	*mz = mcurr0[2];
	}

else if (mode==1)	/* Indicates to find steady-state magnetization */
	{
	scalemat(amat,-1.0);		/* Negate A matrix 	*/
	addmats(amat,imat,amat);	/* Now amat = (I-A)		*/
	invmat(amat,imat);		/* Reuse imat as inv(I-A) 	*/
	multmatvec(imat,bvec,mvec);	/* Now M = inv(I-A)*B		*/
	*mx = mvec[0];
	*my = mvec[1];
	*mz = mvec[2];
	}


}



void blochsimfz(double *b1real, double *b1imag, double *xgrad, double *ygrad, double *zgrad,
		double *tsteps,
		int ntime, double *t1, double *t2, double *dfreq, int nfreq,
		double *dxpos, double *dypos, double *dzpos, int npos,
		double *mx, double *my, double *mz, int mode)


{
int count;
int poscount;
int fcount;
int totpoints;
int totcount = 0;

int ntout;

double t1val;
double t2val;

double *e1;
double *e2;
double *e1ptr;
double *e2ptr;
double *tstepsptr;
double *dxptr, *dyptr, *dzptr;
double *t1ptr, *t2ptr;


if (mode & 2)
	ntout = ntime;
else
	ntout = 1;

// 	/* First calculate the E1 and E2 values at each time step. */
//
// e1 = (double *) malloc(ntime * sizeof(double));
// e2 = (double *) malloc(ntime * sizeof(double));
//
// e1ptr = e1;
// e2ptr = e2;
// tstepsptr = tsteps;
//
// for (count=0; count < ntime; count++)
// 	{
// 	*e1ptr++ = exp(- *tstepsptr / t1);
// 	*e2ptr++ = exp(- *tstepsptr++ / t2);
// 	}

totpoints = npos*nfreq;

for (fcount=0; fcount < nfreq; fcount++) {
    dxptr = dxpos;
    dyptr = dypos;
    dzptr = dzpos;
		t1ptr = t1;
		t2ptr = t2;
    for (poscount=0; poscount < npos; poscount++) {
				/* First calculate the E1 and E2 values at each time step. */

				e1 = (double *) malloc(ntime * sizeof(double));
				e2 = (double *) malloc(ntime * sizeof(double));

				e1ptr = e1;
				e2ptr = e2;
				tstepsptr = tsteps;

				t1val = *t1ptr++;
				t2val = *t2ptr++;
				for (count=0; count < ntime; count++) {
					*e1ptr++ = exp(- *tstepsptr / t1val);
					*e2ptr++ = exp(- *tstepsptr++ / t2val);
				}

				if (mode == 3)	/* Steady state AND record all time points. */

					{	/* First go through and find steady state, then
							repeat as if transient starting at steady st.*/

					blochsim(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, ntime,
						e1, e2, *dfreq, *dxptr, *dyptr,
						*dzptr, mx, my, mz, 1);

					blochsim(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, ntime,
						e1, e2, *dfreq, *dxptr++, *dyptr++,
						*dzptr++, mx, my, mz, 2);
					}
				else
					{
					blochsim(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, ntime,
						e1, e2, *dfreq, *dxptr++, *dyptr++,
						*dzptr++, mx, my, mz, mode);
					}

				mx += ntout;
				my += ntout;
				mz += ntout;

				totcount++;
				if ((totpoints > 40000) && ( ((10*totcount)/totpoints)> (10*(totcount-1)/totpoints) ))
					printf("%d%% Complete.\n",(100*totcount/totpoints));

				free(e1);
				free(e2);
		}
    dfreq++;
}

// free(e1);
// free(e2);

}
