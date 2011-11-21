#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "sepp_vect.h"

void Inverse2x2(double **mat, double **res) {

	double mat_11, mat_12, mat_21, mat_22, hh;

	mat_11 = mat[0][0];

	mat_12 = mat[0][1];

	mat_21 = mat[1][0];

	mat_22 = mat[1][1];

	hh = mat_11*mat_22-mat_12*mat_21;

	res[0][0] = mat_22  / hh;

	res[0][1] = -mat_21 / hh;

	res[1][0] = -mat_12 / hh;

	res[1][1] =  mat_11 / hh;

}


void Inverse3x3(double mat[3][3], double res[3][3]) {

	double mat_11,mat_12,mat_13,mat_21,mat_22,mat_23,mat_31,mat_32,mat_33,det;

	mat_11 = mat[0][0];

	mat_12 = mat[0][1];

	mat_13 = mat[0][2];

	mat_21 = mat[1][0];

	mat_22 = mat[1][1];

	mat_23 = mat[1][2];

	mat_31 = mat[2][0];

	mat_32 = mat[2][1];

	mat_33 = mat[2][2];



	res[0][0] =   mat_33*mat_22-mat_32*mat_23;

	res[0][1] = -(mat_33*mat_21-mat_31*mat_23);

	res[0][2] =   mat_32*mat_21-mat_31*mat_22;

	res[1][0] = -(mat_33*mat_12-mat_32*mat_13);

	res[1][1] =   mat_33*mat_11-mat_31*mat_13;

	res[1][2] = -(mat_32*mat_11-mat_31*mat_12);

	res[2][0] =   mat_23*mat_12-mat_22*mat_13;

	res[2][1] = -(mat_23*mat_11-mat_21*mat_13);

	res[2][2] =   mat_22*mat_11-mat_21*mat_12;



	det  =  mat_11*(res[0][0])-mat_21*(-res[0][1])+mat_31*(res[0][2]);



	res[0][0] /= det;

	res[0][1] /= det;

	res[0][2] /= det;

	res[1][0] /= det;

	res[1][1] /= det;

	res[1][2] /= det;

	res[2][0] /= det;

	res[2][1] /= det;

	res[2][2] /= det;

}


void Inverse3x3diag(double mat[3][3],double vec[3], double res[3][3]) {

	double mat_11, mat_12, mat_13, mat_21, mat_22, mat_23, mat_31, mat_32,
	mat_33, det;

	mat_11 = mat[0][0] + vec[0];

	mat_12 = mat[0][1];

	mat_13 = mat[0][2];

	mat_21 = mat[1][0];

	mat_22 = mat[1][1] + vec[1];

	mat_23 = mat[1][2];

	mat_31 = mat[2][0];

	mat_32 = mat[2][1];

	mat_33 = mat[2][2] + vec[2];



	res[0][0] = mat_33*mat_22-mat_32*mat_23;

	res[0][1] = -(mat_33*mat_21-mat_31*mat_23);

	res[0][2] = mat_32*mat_21-mat_31*mat_22;

	res[1][0] = -(mat_33*mat_12-mat_32*mat_13);

	res[1][1] = mat_33*mat_11-mat_31*mat_13;

	res[1][2] = -(mat_32*mat_11-mat_31*mat_12);

	res[2][0] = mat_23*mat_12-mat_22*mat_13;

	res[2][1] = -(mat_23*mat_11-mat_21*mat_13);

	res[2][2] = mat_22*mat_11-mat_21*mat_12;



	det  =  mat_11*(res[0][0])-mat_21*(-res[0][1])+mat_31*(res[0][2]);



	res[0][0] /= det;

	res[0][1] /= det;

	res[0][2] /= det;

	res[1][0] /= det;

	res[1][1] /= det;

	res[1][2] /= det;

	res[2][0] /= det;

	res[2][1] /= det;

	res[2][2] /= det;

}



SEXP normData(SEXP xS, SEXP E_SX, SEXP nnS){

	double *x=REAL(xS);

	double *E_SX_Data=REAL(E_SX);

	double euclid_Scalar_x,euclid_Scalar_E_SX, xt[2], E_SX_t[3], epsv;

	int nn=(INTEGER(nnS))[0];

	int nn2;

	void** ret=Calloc(1, void*);

	double **E_SX_n;

	E_SX_n	=(double**)Calloc(nn, double*);

	int j,l;


	int *nna;

	epsv = 0.0;

	nna=(int*)Calloc(1, int);

	*nna=nn;

	for(l=0; l<nn; l++) {

		E_SX_n[l]=(double*)Calloc(3, double);

	}

	nn2 = nn + nn;

	for (j=0;j<nn;j++) {

		xt[0] = x[j];

		xt[1] = x[nn+j];

		euclid_Scalar_x = sqrt(pow(xt[0],2) + pow(xt[1],2));

		E_SX_t[0] = E_SX_Data[j] + epsv;

		E_SX_t[1] = E_SX_Data[nn+j] + epsv;

		E_SX_t[2] = E_SX_Data[nn2+j] + epsv;

		euclid_Scalar_E_SX = sqrt(pow(E_SX_t[0],2) + pow(E_SX_t[1],2) +
				pow(E_SX_t[2],2));

		E_SX_n[j][0] = E_SX_t[0] / euclid_Scalar_E_SX * euclid_Scalar_x;

		E_SX_n[j][1] = E_SX_t[1] / euclid_Scalar_E_SX * euclid_Scalar_x;

		E_SX_n[j][2] = E_SX_t[2] / euclid_Scalar_E_SX * euclid_Scalar_x;


	}

	ret[0]=(void*)nna;

	ret[1]=(void*)E_SX_n;


	return R_MakeExternalPtr(ret,R_NilValue,R_NilValue);
}



SEXP sparseFarmsC(SEXP xS, SEXP cycS, SEXP XXS, SEXP nnS) {

	double *x=REAL(xS);

	double *XX=REAL(XXS);



	int cyc=(INTEGER(cycS))[0];

	int nn=(INTEGER(nnS))[0];


	double fac1,f1,f2,epsv,ess;

	double LPsiL[3][3], LPsi[3][2], Psi[2][2], sum1[2][3], sum2[3][3], tt[3],
	tt1[3], xt[2], inva[3][3], invp[2][2], tmplapla[3];

	void** ret=Calloc(4, void*);	

	double **lapla,**E_SX_n, **L;

	int *nna;

	int i,j,k,l;

	nna=(int*)Calloc(1, int);

	*nna=nn;


	L=(double**)Calloc(2, double*);

	lapla=(double**)Calloc(nn, double*);

	E_SX_n	=(double**)Calloc(nn, double*);



	for(l=0; l<2; l++) {

		L[l]=(double*)Calloc(3, double);

	}



	for(l=0; l<nn; l++) {

		lapla[l]=(double*)Calloc(3, double);

		E_SX_n[l]=(double*)Calloc(3, double);

	}



	fac1=0.9;

	f1= sqrt(fac1*XX[0]);

	f2= sqrt(fac1*XX[3]);

	epsv = 1e-10;

	L[0][0] = epsv * f1;

	L[0][1] = 1 / sqrt(2) * f1;

	L[0][2] = f1;

	L[1][0] = f2;

	L[1][1] = 1 / sqrt(2) * f2;

	L[1][2] = epsv * f2;

	Psi[0][0] = (1-fac1)*XX[0];

	Psi[0][1] = 0.0;

	Psi[1][0] = 0.0;

	Psi[1][1] = (1-fac1)*XX[3];

	for(l=0; l<nn; l++) {

		for(k=0; k<3; k++) {

			lapla[l][k] = 1.0;

		}

	}



	invp[0][0] = 1.0/Psi[0][0];

	invp[0][1] = 0.0;

	invp[1][0] = 0.0;

	invp[1][1] = 1.0/Psi[1][1];



	for (i=0;i<cyc;i++){

		MULMat_t(LPsi,L,invp,3,2,2);

		MULMat(LPsiL,LPsi,L,3,2,3);

		for(l=0; l<2; l++) {

			for(k=0; k<3; k++) {

				sum1[l][k] = 0.0;

			}
		}

		for(l=0; l<3; l++) {
			for(k=0; k<3; k++) {

				sum2[l][k] = 0.0;

			}
		}


		for (j=0;j<nn;j++) {

			tmplapla[0] = lapla[j][0];

			tmplapla[1] = lapla[j][1];

			tmplapla[2] = lapla[j][2];

			Inverse3x3diag(LPsiL,tmplapla,inva);

			xt[0] = x[j];

			xt[1] = x[nn+j];

			MULMVr(tt,LPsi,xt,3,2);

			MULMVr(tt1,inva,tt,3,3);

			if (tt1[0]<0.0) tt1[0] = 0.0;

			if (tt1[1]<0.0) tt1[1] = 0.0;

			if (tt1[2]<0.0) tt1[2] = 0.0;

			for(l=0; l<3; l++) {

				for(k=0; k<3; k++) {

					ess = inva[l][k] + tt1[l]*tt1[k];

					sum2[l][k] += ess;

					if (l==k) lapla[j][k]= 1/sqrt(epsv+ess);

				}

				for(k=0; k<2; k++) {

					sum1[k][l] += xt[k]*tt1[l];

				}

			}

		}

		Inverse3x3(sum2,inva);

		MULMat(L,sum1,inva,2,3,3);

		for(l=0; l<2; l++) {

			for(k=0; k<3; k++) {

				if (L[l][k]<0.0) L[l][k] = 0.0;

			}

		}

	}

	MULMat_t(LPsi,L,invp,3,2,2);

	MULMat(LPsiL,LPsi,L,3,2,3);

	for (j=0;j<nn;j++) {

		tmplapla[0] = lapla[j][0];

		tmplapla[1] = lapla[j][1];

		tmplapla[2] = lapla[j][2];

		Inverse3x3diag(LPsiL,tmplapla,inva);


		xt[0] = x[j];
		xt[1] = x[nn+j];

		MULMVr(tt,LPsi,xt,3,2);

		MULMVr(tt1,inva,tt,3,3);


		if (tt1[0]<0.0) tt1[0] = 0.0;

		if (tt1[1]<0.0) tt1[1] = 0.0;

		if (tt1[2]<0.0) tt1[2] = 0.0;



		for(k=0; k<3; k++) {

			E_SX_n[j][k]=tt1[k];

		}

	}



	ret[0]=(void*)nna;

	ret[1]=(void*)L;

	ret[2]=(void*)E_SX_n;

	ret[3]=(void*)lapla;

	return R_MakeExternalPtr(ret,R_NilValue,R_NilValue);

}





SEXP getL(SEXP pointer) {

	void** ret=R_ExternalPtrAddr(pointer);

	double **L=(double**)ret[1];

	int i,j;

	SEXP rueckgabe;

	PROTECT(rueckgabe = allocVector(REALSXP, 6));

	for(j=0; j<3; j++) {
		for(i=0; i<2; i++) {

			REAL(rueckgabe)[i+j*2]=L[i][j];

		}

	}

	UNPROTECT(1);

	return rueckgabe;

}



SEXP getEss(SEXP pointer) {

	void** ret=R_ExternalPtrAddr(pointer);

	double **E_SX_n=(double**)ret[2];

	int nn=*((int*)ret[0]);

	int i,j;

	SEXP rueckgabe;

	PROTECT(rueckgabe = allocVector(REALSXP, 3*nn));

	for(j=0; j<3; j++) {
		for(i=0; i<nn; i++) {

			REAL(rueckgabe)[i+j*nn]=E_SX_n[i][j];

		}

	}

	UNPROTECT(1);

	return rueckgabe;

}



SEXP getLap(SEXP pointer) {

	void** ret=R_ExternalPtrAddr(pointer);

	double **lapla=(double**)ret[3];

	int nn=*((int*)ret[0]);

	int i,j;

	SEXP rueckgabe;

	PROTECT(rueckgabe = allocVector(REALSXP, 3*nn));

	for(i=0; i<nn; i++) {

		for(j=0; j<3; j++) {

			REAL(rueckgabe)[i*3+j]=lapla[i][j];

		}

	}

	UNPROTECT(1);

	return rueckgabe;

}






SEXP getE_SX_norm(SEXP pointer) {

	void** ret=R_ExternalPtrAddr(pointer);

	double **E_SX_n=(double**)ret[1];

	int nn=*((int*)ret[0]);

	int i,j;

	SEXP rueckgabe;

	PROTECT(rueckgabe = allocVector(REALSXP, 3*nn));

	for(j=0; j<3; j++) {
		for(i=0; i<nn; i++) {

			REAL(rueckgabe)[i+j*nn]=E_SX_n[i][j];

		}

	}

	UNPROTECT(1);

	return rueckgabe;

}






SEXP deinit(SEXP pointer) {

	void** ret=R_ExternalPtrAddr(pointer);	

	int nn=*((int*)ret[0]);

	double **L=(double**)ret[1];

	double **E_SX_n=(double**)ret[2];

	double **lapla=(double**)ret[3];

	int l;


	for(l=0; l<2; l++) {

		Free(L[l]);

	}

	for(l=0; l<nn; l++) {

		Free(E_SX_n[l]);

		Free(lapla[l]);

	}

	Free(ret[3]);

	Free(ret[2]);

	Free(ret[1]);

	Free(ret[0]);

	Free(ret);

	return pointer;

}




SEXP deinit_ESX(SEXP pointer) {

	void** ret=R_ExternalPtrAddr(pointer);	

	int nn=*((int*)ret[0]);

	double **E_SX_n=(double**)ret[1];

	int l;

	for(l=0; l<nn; l++) {

		Free(E_SX_n[l]);

	}

	Free(ret[1]);

	Free(ret[0]);

	Free(ret);

	return pointer;

}

