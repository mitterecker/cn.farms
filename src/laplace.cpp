#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <Rmath.h>
#include <complex>
#include <errno.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <string>
#include <time.h>

using namespace std;

const int steps=30;
const double borderFactor=0.005;
const double integrationFactor=1e-50;
const double integrationFailureFactor=5.0;
const double sigmaFactor=0.1;
double debug=false;

const double Pi=M_PI;
const double E=M_E;
const double sqrt3=1.7320508075688772935274463415058723669428052538104;
const double sqrt2Pi=2.5066282746310005024157652848110452530069867406099;
const double logSqrt2Pi=0.91893853320467274178032973640561763986139747363778;

double eps1=1e-10;
double eps2=1e-10;


double logSigmaFactor;
double log1MSigmaFactor;
double sqrtMLogSigmaFactor;
double sqrtM2LogSigmaFactor;
double logBorderFactor;
double logIntegrationFactor;
double logMultiplyFactor=log(0.01);
ofstream debugOutput;


void quadratic(double a, double b, double c, double& D, double& z1,
		double& z2) {
	double maxCoeff=fmax(fmax(fabs(a), fabs(b)), fabs(c));
	a=a/maxCoeff;
	b=b/maxCoeff;
	c=c/maxCoeff;
	D=b*b-4.*a*c;
	if(D<0) {
		z1=-b/(2.*a);
		z2=sqrt(-D)/(2.*a);
	} else {
		double sqrtD=sqrt(D);
		z1=(-b-sqrtD)/(2.*a);
		z2=(-b+sqrtD)/(2.*a);
		if(z2<z1) {
			double zsave=z1;
			z1=z2;
			z2=zsave;
		}
	}
}

double evalLogGauss(double z, double my, double sigma) {
	return -((-my + z)*(-my + z))/(2.*sigma*sigma) - log(2*Pi)/2. - log(sigma);
}

double evalLogUnnormalizedPosterior(double a, double b, double c, double sigmaZ,
		double logNormfact, double z) {
	double z2=z*z;
	return (a*z2+b*z+c-fabs(z)/sigmaZ)+logNormfact;
}

double evalApproximation(double leftWeight, double rightWeight, double leftMy,
		double rightMy, double leftSigma, double rightSigma, double z) {
	if(z<0) {
		return log(leftWeight)+evalLogGauss(z, leftMy, leftSigma);
	}
	else {
		return log(rightWeight)+evalLogGauss(z, rightMy, rightSigma);
	}
}


void computeParameters(double a, double b, double c, double sigmaZ,
		double logNormfact, int& Case, double& maxZ, double& logMaxValue,
		double &alpha, double& leftWeight, double& rightWeight, double&leftMy,
		double& rightMy, double& leftSigma, double& rightSigma,double& moment1,
		double& moment2, double& entropy, double& crossentropy) {

	//"sigmaZ" is considered to be "b" in a Laplace distribution

	//Cases:
	//1: Laplace
	//2: Gauss
	//3: Mixt. of Truncated Gauss

	if(fabs(a)<=eps1) {
		Case=1;
		moment1=0.0;
		moment2=2.0*sigmaZ*sigmaZ;
		entropy=log(2.0*exp(1.0)*sigmaZ);
		crossentropy=log(2.0*exp(1.0)*sigmaZ);
		logMaxValue=evalLogUnnormalizedPosterior(a, b, c, sigmaZ,
				logNormfact, 0.0);
		alpha=exp(logMaxValue+log(2.0*sigmaZ));
		maxZ=0.0;
		return;
	}

	leftSigma=sqrt(-1.0/(2.0*a));
	rightSigma=sqrt(-1.0/(2.0*a));
	leftMy=(-b-1.0/sigmaZ)/(2.0*a);
	rightMy=(-b+1.0/sigmaZ)/(2.0*a);

	double areaLeft=0.5*erfc(-(-leftMy/leftSigma)/sqrt(2.0));
	double areaRight=1.0-0.5*erfc(-(-rightMy/rightSigma)/sqrt(2.0));

	//Special Case: Zero == Laplace
	if(areaLeft<=eps2&&areaRight<=eps2) {
		Case=1;
		moment1=0.0;
		moment2=2.0*sigmaZ*sigmaZ;
		entropy=log(2.0*exp(1.0)*sigmaZ);
		crossentropy=log(2.0*exp(1.0)*sigmaZ);
		logMaxValue=evalLogUnnormalizedPosterior(a, b, c, sigmaZ,
				logNormfact, 0.0);
		alpha=exp(logMaxValue+log(2.0*sigmaZ));
		maxZ=0.0;
		return;
	}

	//Special Case: Outlier == Gauss
	if(areaLeft<=eps2) {
		//if(areaLeft==0.0) {
		Case=2;
		moment1=rightMy;
		moment2=rightMy*rightMy+rightSigma*rightSigma;
		entropy=0.5*log(2.0*Pi*exp(1.0)*rightSigma*rightSigma);
		crossentropy=-(-log(2.0*sigmaZ)-moment1/sigmaZ);
		leftWeight=0.0;
		rightWeight=1.0;
		logMaxValue=evalLogUnnormalizedPosterior(a, b, c, sigmaZ, logNormfact,
				rightMy);
		alpha=exp(logMaxValue+log(rightSigma*sqrt(2*Pi)));
		maxZ=rightMy;
		return;
	}

	//Special Case: Outlier == Gauss
	if(areaRight<=eps2) {
		//if(areaRight==0.0) {
		Case=2;
		moment1=leftMy;
		moment2=leftMy*leftMy+leftSigma*leftSigma;
		entropy=0.5*log(2.0*Pi*exp(1.0)*leftSigma*leftSigma);
		crossentropy=-(-log(2.0*sigmaZ)+moment1/sigmaZ);
		leftWeight=1.0;
		rightWeight=0.0;
		logMaxValue=evalLogUnnormalizedPosterior(a, b, c, sigmaZ, logNormfact,
				leftMy);
		alpha=exp(logMaxValue+log(leftSigma*sqrt(2*Pi)));
		maxZ=leftMy;
		return;
	}

	Case=3;

	double weight=exp(evalLogGauss(0, rightMy, rightSigma)-
			evalLogGauss(0, leftMy, leftSigma)+log(areaLeft)-log(areaRight));
	leftWeight=weight/(weight+1.0);
	rightWeight=1.0/(weight+1.0);

	//Computation of Moments and Entropy
	double I0left=1.0;
	double I0right=1.0;
	double I1left=-exp(evalLogGauss(-leftMy/leftSigma, 0.0, 1.0))/areaLeft;
	double I1right=exp(evalLogGauss(-rightMy/rightSigma, 0.0, 1.0))/areaRight;
	double I2left=(-leftMy/leftSigma)*I1left+I0left;
	double I2right=(-rightMy/rightSigma)*I1right+I0right;

	double leftMoment1=1.0*leftMy*I0left+1.0*leftSigma*I1left;
	double rightMoment1=1.0*rightMy*I0right+1.0*rightSigma*I1right;
	double leftMoment2=1.0*leftMy*leftMy*I0left+2.0*leftMy*leftSigma*I1left+
			1.0*leftSigma*leftSigma*I2left;
	double rightMoment2=1.0*rightMy*rightMy*I0right+
			2.0*rightMy*rightSigma*I1right+1.0*rightSigma*rightSigma*I2right;

	moment1=leftWeight*(1.0*leftMy*I0left+1.0*leftSigma*I1left)+rightWeight*
			(1.0*rightMy*I0right+1.0*rightSigma*I1right);
	moment2=leftWeight*(1.0*leftMy*leftMy*I0left+2.0*leftMy*leftSigma*I1left+
			1.0*leftSigma*leftSigma*I2left)+rightWeight*
			(1.0*rightMy*rightMy*I0right+2.0*rightMy*rightSigma*I1right+1.0*
					rightSigma*rightSigma*I2right);

	double leftEntropy=-(-log(sqrt(2.0*Pi)*leftSigma*areaLeft)+
			(-leftMoment2/2.0+leftMoment1*leftMy-leftMy*leftMy/2.0)/
			(leftSigma*leftSigma));
	double rightEntropy=-(-log(sqrt(2.0*Pi)*rightSigma*areaRight)+
			(-rightMoment2/2.0+rightMoment1*rightMy-rightMy*rightMy/2.0)/
			(rightSigma*rightSigma));
	entropy=-leftWeight*log(leftWeight)-rightWeight*log(rightWeight)+
			leftWeight*leftEntropy+rightWeight*rightEntropy;

	crossentropy=-(-log(2.0*sigmaZ)+
			(leftWeight*leftMoment1-rightWeight*rightMoment1)/sigmaZ);

	//Computation of MaxZ
	double max1Z;
	if(leftMy<0.0)
		max1Z=leftMy;
	else
		max1Z=0.0;
	double max2Z;
	if(rightMy>0.0)
		max2Z=rightMy;
	else
		max2Z=0.0;
	if(log(leftWeight/areaLeft)+evalLogGauss(max1Z, leftMy, leftSigma)>
	log(rightWeight/areaRight)+evalLogGauss(max2Z, rightMy, rightSigma))
		maxZ=max1Z;
	else
		maxZ=max2Z;
	logMaxValue=evalLogUnnormalizedPosterior(a, b, c, sigmaZ, logNormfact, maxZ);
	alpha=exp(logMaxValue-evalApproximation(leftWeight/areaLeft,
			rightWeight/areaRight, leftMy, rightMy, leftSigma, rightSigma, maxZ));
}


extern "C" SEXP momentsGauss(SEXP it, SEXP eps1S, SEXP eps2S, SEXP aS, SEXP bS,
		SEXP cS, SEXP sigmaZS, SEXP normfactS, SEXP evalProb, SEXP methodS) {
	logBorderFactor=log(borderFactor);
	logIntegrationFactor=log(integrationFactor);

	int elements=LENGTH(aS);
	int method=0;

	SEXP RET;
	SEXP normConstRET, moment1RET, moment2RET, maxRET, entropyRET,
	crossentropyRET, caseRET;

	PROTECT(RET = allocVector(VECSXP, 7));
	PROTECT(normConstRET = allocVector(REALSXP, elements));
	PROTECT(moment1RET = allocVector(REALSXP, elements));
	PROTECT(moment2RET = allocVector(REALSXP, elements));
	PROTECT(maxRET = allocVector(REALSXP, elements));
	PROTECT(entropyRET = allocVector(REALSXP, elements));
	PROTECT(crossentropyRET = allocVector(REALSXP, elements));
	PROTECT(caseRET = allocVector(REALSXP, elements));

	eps1=(double)(REAL(eps1S)[0]);
	eps2=(double)(REAL(eps2S)[0]);

	for(int j=0; j<elements; j++) {

		double a=(double)(REAL(aS)[j]);
		double b=(double)(REAL(bS)[j]);
		double c=(double)(REAL(cS)[j]);
		double sigmaZ=(double) ((REAL(sigmaZS))[0]);
		double logNormfact=log((double)(REAL(normfactS)[j]));

		int Case;
		double maxZ;
		double logMaxValue;
		double alpha;
		double leftWeight, rightWeight;
		double leftMy, rightMy;
		double leftSigma, rightSigma;
		double moment1;
		double moment2;
		double entropy;
		double crossentropy;


		computeParameters(a, b, c, sigmaZ, logNormfact, 
				Case,
				maxZ, logMaxValue,
				alpha, leftWeight, rightWeight, leftMy, rightMy, leftSigma,
				rightSigma,	moment1, moment2, entropy, crossentropy);


		if((double)(REAL(evalProb)[0])==1.0)
			REAL(normConstRET)[j]=alpha;
		else
			REAL(normConstRET)[j]=1.0;
		REAL(moment1RET)[j]=moment1;
		REAL(moment2RET)[j]=moment2;
		REAL(maxRET)[j]=maxZ;
		REAL(entropyRET)[j]=entropy;
		REAL(crossentropyRET)[j]=crossentropy;
		REAL(caseRET)[j]=(double)Case;
	}

	SET_VECTOR_ELT(RET, 0, normConstRET);
	SET_VECTOR_ELT(RET, 1, moment1RET);
	SET_VECTOR_ELT(RET, 2, moment2RET);
	SET_VECTOR_ELT(RET, 3, maxRET);
	SET_VECTOR_ELT(RET, 4, entropyRET);
	SET_VECTOR_ELT(RET, 5, crossentropyRET);
	SET_VECTOR_ELT(RET, 6, caseRET);

	SEXP namesRET;
	PROTECT(namesRET = allocVector(STRSXP, 7));
	SET_STRING_ELT(namesRET, 0, mkChar("normConst"));
	SET_STRING_ELT(namesRET, 1, mkChar("moment1"));
	SET_STRING_ELT(namesRET, 2, mkChar("moment2"));
	SET_STRING_ELT(namesRET, 3, mkChar("max"));
	SET_STRING_ELT(namesRET, 4, mkChar("Entropy"));
	SET_STRING_ELT(namesRET, 5, mkChar("CrossEntropy"));
	SET_STRING_ELT(namesRET, 6, mkChar("Case"));

	setAttrib(RET, R_NamesSymbol, namesRET);
	UNPROTECT(9);
	return RET;
}

void evalIntegration(double z, double a, double b, double c, double sigmaZ,
		double logNormfact, double logFmax, double& newLogFmax, double& value,
		double& entropy, double& crossentropy) {
	double eval=evalLogUnnormalizedPosterior(a, b, c, sigmaZ, logNormfact, z);
	if(eval>newLogFmax)
		newLogFmax=eval;
	value=exp(eval-logFmax);
	entropy=value*(eval-logFmax);
	crossentropy=value*(-log(sigmaZ)-fabs(z)/sigmaZ);
}

void sumTrapez(double a, double b, double eps, int n, double* integral,
		double* err, bool& fail,
		void func(double val, double apar, double bpar, double cpar,
				double sigmapar, double normfactpar, double logFmax,
				double& newLogFmax, double& value, double& entropy,
				double& crossentropy),
				double ap, double bp, double cp, double sigmap,
				double normfactp, double logFmax, double& newLogFmax) {
	integral[0]=0.0;
	integral[1]=0.0;
	integral[2]=0.0;
	integral[3]=0.0;
	integral[4]=0.0;


	double h=(b-a)/((double)n);
	for(int i=0; i<n; i++) {
		double xvalue=a+h*((2.*i-1.)/2.);
		double addmoment;
		double addentropy;
		double addcrossentropy;
		func(xvalue, ap, bp, cp, sigmap, normfactp, logFmax, newLogFmax,
				addmoment, addentropy, addcrossentropy);
		integral[0]=integral[0]+addmoment;
		addmoment=addmoment*xvalue;
		integral[1]=integral[1]+addmoment;
		addmoment=addmoment*xvalue;
		integral[2]=integral[2]+addmoment;
		integral[3]=integral[3]+addentropy;
		integral[4]=integral[4]+addcrossentropy;
	}
	integral[0]=integral[0]*h;
	integral[1]=integral[1]*h;
	integral[2]=integral[2]*h;
	integral[3]=integral[3]*h;
	integral[4]=integral[4]*h;
}

extern "C" SEXP momentsIntegrationTrapez(SEXP it, SEXP eps1S, SEXP eps2S,
		SEXP aS, SEXP bS, SEXP cS, SEXP sigmaZS, SEXP normfactS,
		SEXP evalProb, SEXP methodS) {
	logBorderFactor=log(borderFactor);
	logIntegrationFactor=log(integrationFactor);

	logBorderFactor=logIntegrationFactor;

	int elements=LENGTH(aS);
	int method=0;

	SEXP RET;
	SEXP normConstRET, moment1RET, moment2RET, maxRET, entropyRET,
	crossentropyRET, caseRET;

	PROTECT(RET = allocVector(VECSXP, 7));
	PROTECT(normConstRET = allocVector(REALSXP, elements));
	PROTECT(moment1RET = allocVector(REALSXP, elements));
	PROTECT(moment2RET = allocVector(REALSXP, elements));
	PROTECT(maxRET = allocVector(REALSXP, elements));
	PROTECT(entropyRET = allocVector(REALSXP, elements));
	PROTECT(crossentropyRET = allocVector(REALSXP, elements));
	PROTECT(caseRET = allocVector(REALSXP, elements));

	eps1=(double)(REAL(eps1S)[0]);
	eps2=(double)(REAL(eps2S)[0]);

	for(int j=0; j<elements; j++) {

		double a=(double)(REAL(aS)[j]);
		double b=(double)(REAL(bS)[j]);
		double c=(double)(REAL(cS)[j]);
		double sigmaZ=(double) ((REAL(sigmaZS))[0]);
		double logNormfact=log((double)(REAL(normfactS)[j]));

		int Case;
		double maxZ;
		double logMaxValue;
		double alpha;
		double leftWeight, rightWeight;
		double leftMy, rightMy;
		double leftSigma, rightSigma;
		double moment1;
		double moment2;
		double entropy;
		double crossentropy;

		computeParameters(a, b, c, sigmaZ, logNormfact,
				Case,
				maxZ, logMaxValue,
				alpha, leftWeight, rightWeight, leftMy, rightMy, leftSigma,
				rightSigma,	moment1, moment2, entropy, crossentropy);

		if(method==1||method==2) {
			logNormfact=logNormfact+log(2.0);
		}
		double left, right, D, variable;
		quadratic(a, b+1.0/sigmaZ, c-logBorderFactor-logMaxValue+logNormfact,
				D, left, variable);
		quadratic(a, b-1.0/sigmaZ, c-logBorderFactor-logMaxValue+logNormfact,
				D, variable, right);

		if(method==1) {
			left=0.0;
		}
		else if(method==2) {
			right=0.0;
		}


		double logFmax=evalLogUnnormalizedPosterior(a, b, c, sigmaZ,
				logNormfact, maxZ);
		double newLogFmax=logFmax;


		double integral[5];
		double err[5];
		bool fail;
		sumTrapez(left, right, numeric_limits<double>::epsilon(), 15000,
				integral, err, fail, evalIntegration, a, b, c, sigmaZ,
				logNormfact, logFmax, newLogFmax);

		bool repeatQuadrature=false;
		for(int i=0; i<3; i++) {
			if(isnan(integral[i]))
				repeatQuadrature=true;
		}
		if(repeatQuadrature) {
			sumTrapez(left, right, numeric_limits<double>::epsilon(), 15000,
					integral, err, fail, evalIntegration, a, b, c, sigmaZ,
					logNormfact, newLogFmax, newLogFmax);
			logFmax=newLogFmax;
		}

		double normConst=integral[0];
		if((double)(REAL(evalProb)[0])==1.0)
			REAL(normConstRET)[j]=exp(log(normConst)+logFmax);
		else
			REAL(normConstRET)[j]=1.0;
		if(normConst==0.0)
			normConst=1.0;
		REAL(moment1RET)[j]=integral[1]/normConst;
		REAL(moment2RET)[j]=integral[2]/normConst;
		REAL(entropyRET)[j]=-(integral[3]/normConst-log(normConst));
		REAL(crossentropyRET)[j]=-(integral[4]/normConst);
		REAL(caseRET)[j]=(double)Case;
	}

	SET_VECTOR_ELT(RET, 0, normConstRET);
	SET_VECTOR_ELT(RET, 1, moment1RET);
	SET_VECTOR_ELT(RET, 2, moment2RET);
	SET_VECTOR_ELT(RET, 3, maxRET);
	SET_VECTOR_ELT(RET, 4, entropyRET);
	SET_VECTOR_ELT(RET, 5, crossentropyRET);
	SET_VECTOR_ELT(RET, 6, caseRET);

	SEXP namesRET;
	PROTECT(namesRET = allocVector(STRSXP, 7));
	SET_STRING_ELT(namesRET, 0, mkChar("normConst"));
	SET_STRING_ELT(namesRET, 1, mkChar("moment1"));
	SET_STRING_ELT(namesRET, 2, mkChar("moment2"));
	SET_STRING_ELT(namesRET, 3, mkChar("max"));
	SET_STRING_ELT(namesRET, 4, mkChar("Entropy"));
	SET_STRING_ELT(namesRET, 5, mkChar("CrossEntropy"));
	SET_STRING_ELT(namesRET, 6, mkChar("Case"));

	setAttrib(RET, R_NamesSymbol, namesRET);
	UNPROTECT(9);
	return RET;
}
