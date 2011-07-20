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
const double sqrtPi=1.7724538509055160272981674833411451827975494561224;
const double sqrt2=1.4142135623730950488016887242096980785696718753769;
const double sqrt3=1.7320508075688772935274463415058723669428052538104;
const double sqrt2Pi=2.5066282746310005024157652848110452530069867406099;
const double sqrt2DPi=0.79788456080286535587989211986876373695171726232987;
const double logarithm2=0.69314718055994530941723212145817656807550013436026;
const double log2Pi=1.8378770664093454835606594728112352797227949472756;
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


void quadratic(double a, double b, double c, double& D, double& z1, double& z2) {
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

double evalLogGauss(double z, double my, double sigma, double logSigma) {
	return -((-my + z)*(-my + z))/(2.*sigma*sigma) - log2Pi/2. - logSigma;
}

double evalLogUnnormalizedPosterior(double a, double b, double c, double sigmaZ, double logNormfact, double z) {
	double z2=z*z;
	return (a*z2+b*z+c-fabs(z)/sigmaZ)+logNormfact;
}

double evalApproximation(double logLeftWeight, double logRightWeight, double leftMy, double rightMy, double leftSigma, double logLeftSigma, double rightSigma, double logRightSigma, double z) {
	if(z<0) {
		return logLeftWeight+evalLogGauss(z, leftMy, leftSigma, logLeftSigma);
	}
	else {
		return logRightWeight+evalLogGauss(z, rightMy, rightSigma, logRightSigma);
	}
}


/*
double asymp(double x) {
  double sumx=1.0;
  double x2=x*x;
  double zaehler=1.0;
  double nenner=2*x2;
  for(int i=1; i<100; i++) {
    if(i%2==0) {
      sumx=sumx+zaehler/nenner;
    }
    else {
      sumx=sumx-zaehler/nenner;
    }
    nenner=nenner*2*x2;
    zaehler=zaehler*(2*i+1);
  }
  sumx=sumx/(x*sqrt(Pi))
  return sumx;
}

*/

double chainErfc(double x)
{
  double eps=1E-15;
  double bound=1E-30;

  double a;
  double delta;

  double b=x;
  double f=b;
  if(fabs(b)<bound)
    f=bound;
  double c=f;
  double d=0.0;
  int i;
  for(i=1; i<=20; i++) {
    a=i*0.5;
    d=a*d+b;
    if(fabs(d)<bound)
      d=bound;
    c=b+a/c;
    if(fabs(c)<bound)
      c=bound;
    d=1.0/d;
    delta=d*c;
    f *= delta;
    if(fabs(delta-1.0)<=eps) break;
  }
  //printf("%d Iterations\n", i);
  return 1/(f*sqrtPi);
}






void computeParameters(double a, double b, double c, double sigmaZ, double logSigmaZ,
		double logNormfact, int& Case, double& maxZ, double& logMaxValue,
		double &alpha, double& leftWeight, double& rightWeight, double&leftMy,
		double& rightMy, double& leftSigma, double& rightSigma,double& moment1,
		double& moment2, double& entropy, double& crossentropy) {

	//"sigmaZ" is considered to be "b" in a Laplace distribution
	double thres=5.0;

	if(fabs(a) < 1E-4 && fabs(b) < 1E-2) {
		Case=0;
		//Rprintf("CaseSpecial\n");
		maxZ=0.0;
		logMaxValue=evalLogUnnormalizedPosterior(a, b, c, sigmaZ, logNormfact, 0.0);
		alpha=exp(logMaxValue)*2.0*sigmaZ;
		moment1=0.0;
		moment2=2.0*sigmaZ*sigmaZ;
		//entropy=log(2.0*exp(1.0)*sigmaZ);
		entropy=logarithm2+1.0+logSigmaZ;
		//crossentropy=log(2.0*exp(1.0)*sigmaZ);
		crossentropy=logarithm2+1.0+logSigmaZ;
		return;
	}

	double sqrtMa=sqrt(-a);
	//leftSigma=sqrt(-1.0/(2.0*a));
	leftSigma=1.0/(sqrt2*sqrtMa);
	//rightSigma=sqrt(-1.0/(2.0*a));
	rightSigma=leftSigma;
	double logLeftSigma=log(leftSigma);
	double logRightSigma=logLeftSigma;
	leftMy=(-b-1.0/sigmaZ)/(2.0*a);
	rightMy=(-b+1.0/sigmaZ)/(2.0*a);

	double argLeft=(b+(1.0/sigmaZ))/(2.0*sqrtMa);
	double argRight=(b-(1.0/sigmaZ))/(2.0*sqrtMa);
	
	double areaLeft=0.5*erfc(argLeft);
	//double areaRight=1.0-0.5*erfc(argRight);
	double areaRight=0.5*erfc(-argRight);
	double logAreaLeft=log(areaLeft);
	double logAreaRight=log(areaRight);

	double I1left;
	double I1right;
	
	if(argLeft >= 0.0 && argRight <= 0.0) {
		Case=100;
		if(argLeft > thres) {
			Case=Case+10;
			double chainLeft=chainErfc(argLeft);
			//I1left=-(sqrt(2.0))/(sqrt(Pi)*chainLeft);
			I1left=-sqrt2DPi/chainLeft;
			//Rprintf("1..logAreaLeft OLD: %.50lf\n", logAreaLeft);
			logAreaLeft=-argLeft*argLeft+log(0.5*chainLeft);
			//Rprintf("1..logAreaLeft NEW: %.50lf\n", logAreaLeft);
		}
		else {
			//I1left=-exp(evalLogGauss(-leftMy/leftSigma, 0.0, 1.0, 0.0))/areaLeft;
			//I1left=-exp(-argLeft*argLeft-log(2*Pi)/2.)/areaLeft;
			I1left=-exp(-argLeft*argLeft-log2Pi/2.)/areaLeft;
		}
		
		if(argRight < -thres) {
			Case=Case+1;
			double chainRight=chainErfc(-argRight);
			//I1right=(sqrt(2.0))/(sqrt(Pi)*chainRight);
			I1right=sqrt2DPi/chainRight;
			//Rprintf("1..logAreaRight OLD: %.50lf\n", logAreaRight);
			logAreaRight=-argRight*argRight+log(0.5*chainRight);
			//Rprintf("1..logAreaRight NEW: %.50lf\n", logAreaRight);
		}
		else {
			//I1right=exp(evalLogGauss(-rightMy/rightSigma, 0.0, 1.0, 0.0))/areaRight;
			//I1right=exp(-argRight*argRight-log(2*Pi)/2.)/areaRight;
			I1right=exp(-argRight*argRight-log2Pi/2.)/areaRight;
		}
	}
	else if(argLeft <= 0.0 && argRight <= 0.0) {
		if(areaLeft >= 0.99) {
			//Rprintf("Case21\n");
			Case=210;
			maxZ=leftMy;
			logMaxValue=evalLogUnnormalizedPosterior(a, b, c, sigmaZ, logNormfact, leftMy);
			//alpha=exp(logMaxValue+log(leftSigma*sqrt(2*Pi)));
			alpha=exp(logMaxValue+logLeftSigma+logSqrt2Pi);
			leftWeight=1.0;
			rightWeight=0.0;
			moment1=leftMy;
			moment2=leftMy*leftMy+leftSigma*leftSigma;
			//entropy=0.5*log(2.0*Pi*exp(1.0)*leftSigma*leftSigma);
			entropy=0.5*log2Pi+0.5+logLeftSigma;
			//crossentropy=-(-log(2.0*sigmaZ)+moment1/sigmaZ);
			crossentropy=logarithm2+logSigmaZ-moment1/sigmaZ;
			return;
		}
		else {
			//Rprintf("Case22\n");
			Case=220;
			//I1left=-exp(evalLogGauss(-leftMy/leftSigma, 0.0, 1.0, 0.0))/areaLeft;
			//I1left=-exp(-argLeft*argLeft-log(2*Pi)/2.)/areaLeft;
			I1left=-exp(-argLeft*argLeft-log2Pi/2.)/areaLeft;
			
			if(argRight < -thres) {
				Case=Case+1;
				double chainRight=chainErfc(-argRight);
				//I1right=(sqrt(2.0))/(sqrt(Pi)*chainRight);
				I1right=sqrt2DPi/chainRight;
				//Rprintf("22..logAreaRight OLD: %.50lf\n", logAreaRight);
				logAreaRight=-argRight*argRight+log(0.5*chainRight);
				//Rprintf("22..logAreaRight NEW: %.50lf\n", logAreaRight);
			}
			else {
				//I1right=exp(evalLogGauss(-rightMy/rightSigma, 0.0, 1.0, 0.0))/areaRight;
				//I1right=exp(-argRight*argRight-log(2*Pi)/2.)/areaRight;
				I1right=exp(-argRight*argRight-log2Pi/2.)/areaRight;
			}
		}
	}
	else if(argLeft >= 0.0 && argRight >= 0.0) {
		if(areaRight >= 0.99) {
			//Rprintf("Case31\n");
			Case=310;
			maxZ=rightMy;
			logMaxValue=evalLogUnnormalizedPosterior(a, b, c, sigmaZ, logNormfact, rightMy);
			alpha=exp(logMaxValue+logRightSigma+logSqrt2Pi);
			leftWeight=0.0;
			rightWeight=1.0;
			moment1=rightMy;
			moment2=rightMy*rightMy+rightSigma*rightSigma;
			//entropy=0.5*log(2.0*Pi*exp(1.0)*rightSigma*rightSigma);
			entropy=0.5*log2Pi+0.5+logRightSigma;
			//crossentropy=-(-log(2.0*sigmaZ)-moment1/sigmaZ);
			crossentropy=logarithm2+logSigmaZ+moment1/sigmaZ;
			maxZ=rightMy;
			return;
		}
		else {
			//Rprintf("Case32\n");
			Case=320;
			if(argLeft>thres) {
				Case=Case+1;
				double chainLeft=chainErfc(argLeft);
				//I1left=-(sqrt(2.0))/(sqrt(Pi)*chainLeft);
				I1left=-sqrt2DPi/chainLeft;
				//Rprintf("32..logAreaLeft OLD: %.50lf\n", logAreaLeft);
				logAreaLeft=-argLeft*argLeft+log(0.5*chainLeft);
				//Rprintf("32..logAreaLeft NEW: %.50lf\n", logAreaLeft);
			}
			else {
				//I1left=-exp(evalLogGauss(-leftMy/leftSigma, 0.0, 1.0, 0.0))/areaLeft;
				//I1left=-exp(-argLeft*argLeft-log(2*Pi)/2.)/areaLeft;
				I1left=-exp(-argLeft*argLeft-log2Pi/2.)/areaLeft;
			}
			//I1right=exp(evalLogGauss(-rightMy/rightSigma, 0.0, 1.0, 0.0))/areaRight;
			//I1right=exp(-argRight*argRight-log(2*Pi)/2.)/areaRight;
			I1right=exp(-argRight*argRight-log2Pi/2.)/areaRight;
		}
	}
	
	//Rprintf("Case%d\n", Case);
	
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

	if(max1Z!=0.0)
		maxZ=max1Z;
	if(max2Z!=0.0)
		maxZ=max2Z;
	if(max1Z==0.0&&max2Z==0.0)
		maxZ=0.0;
	logMaxValue=evalLogUnnormalizedPosterior(a, b, c, sigmaZ, logNormfact, maxZ);
	
	
	//double weight=-I1right/I1left;
	//leftWeight=weight/(weight+1.0);
	//rightWeight=1.0/(weight+1.0);

	double I0left=1.0;
	double I0right=1.0;
	
	double I2left=(-leftMy/leftSigma)*I1left+I0left;
	double I2right=(-rightMy/rightSigma)*I1right+I0right;

	double leftMoment1=1.0*leftMy*I0left+1.0*leftSigma*I1left;
	double rightMoment1=1.0*rightMy*I0right+1.0*rightSigma*I1right;
	double leftMoment2=1.0*leftMy*leftMy*I0left+2.0*leftMy*leftSigma*I1left+1.0*leftSigma*leftSigma*I2left;
	double rightMoment2=1.0*rightMy*rightMy*I0right+2.0*rightMy*rightSigma*I1right+1.0*rightSigma*rightSigma*I2right;
	
	
	leftWeight=I1right/(I1right-I1left);
	rightWeight=I1left/(I1left-I1right);
	
	double logLeftWeight=log(leftWeight);
	double logRightWeight=log(rightWeight);
	alpha=exp(logMaxValue-evalApproximation(logLeftWeight-logAreaLeft, logRightWeight-logAreaRight, leftMy, rightMy, leftSigma, logLeftSigma, rightSigma, logRightSigma, maxZ));
	
	
	//moment1=leftWeight*(1.0*leftMy*I0left+1.0*leftSigma*I1left)+rightWeight*(1.0*rightMy*I0right+1.0*rightSigma*I1right);
	//Rprintf("moment1 OLD: %.50lf\n", moment1);
	//moment1=(leftMy*I1right-rightMy*I1left)/(I1right-I1left);
	moment1=rightMy-leftWeight/(a*sigmaZ);
        //Rprintf("moment1 NEW: %.50lf\n", moment1);
	
	
	//moment2=leftWeight*(1.0*leftMy*leftMy*I0left+2.0*leftMy*leftSigma*I1left+1.0*leftSigma*leftSigma*I2left)+rightWeight*(1.0*rightMy*rightMy*I0right+2.0*rightMy*rightSigma*I1right+1.0*rightSigma*rightSigma*I2right);
	//Rprintf("moment2 OLD: %.50lf\n", moment2);
	//moment2 = (I1right*leftMy*leftMy+I1right*I1left*leftMy*leftSigma-I1left*rightMy*rightMy-I1left*I1right*rightMy*rightSigma)/(I1right-I1left)+leftSigma*leftSigma; 
	//moment2 = rightMy*rightMy+(b*leftWeight)/(a*a*sigmaZ)+leftSigma*leftSigma+(leftWeight*I1left)/(sqrtMa*sqrtMa*sqrtMa*sqrt(2.0)*sigmaZ);
        moment2 = (-0.5*a+0.25*(b-1/sigmaZ)*(b-1/sigmaZ)+(b*leftWeight)/(sigmaZ)+(leftWeight*I1left*sqrtMa)/(sqrt(2.0)*sigmaZ))/(a*a);
	//Rprintf("moment2 NEW: %.50lf\n", moment2);
	
	//double leftEntropy=-(-log(sqrt(2.0*Pi)*leftSigma)-logAreaLeft+(-leftMoment2/2.0+leftMoment1*leftMy-leftMy*leftMy/2.0)/(leftSigma*leftSigma));
	double leftEntropy=logSqrt2Pi+logLeftSigma+logAreaLeft-(-leftMoment2/2.0+leftMoment1*leftMy-leftMy*leftMy/2.0)/(leftSigma*leftSigma);
	//double rightEntropy=-(-log(sqrt(2.0*Pi)*rightSigma)-logAreaRight+(-rightMoment2/2.0+rightMoment1*rightMy-rightMy*rightMy/2.0)/(rightSigma*rightSigma));
	double rightEntropy=logSqrt2Pi+logRightSigma+logAreaRight-(-rightMoment2/2.0+rightMoment1*rightMy-rightMy*rightMy/2.0)/(rightSigma*rightSigma);
	
	//entropy=-leftWeight*log(leftWeight)-rightWeight*log(rightWeight)+leftWeight*leftEntropy+rightWeight*rightEntropy;
	entropy=-leftWeight*logLeftWeight-rightWeight*logRightWeight+leftWeight*leftEntropy+rightWeight*rightEntropy;
	//crossentropy=-(-log(2.0*sigmaZ)+(leftWeight*leftMoment1-rightWeight*rightMoment1)/sigmaZ);
	crossentropy=logarithm2+logSigmaZ-(leftWeight*leftMoment1-rightWeight*rightMoment1)/sigmaZ;
	
	/*
	Rprintf("a:          %50.25lf\n", a);
	Rprintf("b:          %50.25lf\n", b);
	Rprintf("I1right:    %50.25lf\n", I1right);
	Rprintf("I1left:     %50.25lf\n", I1left);
	Rprintf("leftMy:     %50.25lf\n", leftMy);
	Rprintf("leftSigma:  %50.25lf\n", leftSigma);
	Rprintf("rightMy:    %50.25lf\n", rightMy);
	Rprintf("rightSigma: %50.25lf\n", rightSigma);
	*/
	
	
	/*Rprintf("argLeft:      %29.25lf   argRight:     %29.25lf\n", argLeft, argRight);
	Rprintf("areaLeft:     %29.25lf   areaRight:    %29.25lf\n", areaLeft, areaRight);
	Rprintf("leftMy:       %29.25lf   rightMy:      %29.25lf\n", leftMy, rightMy);
	Rprintf("leftSigma:    %29.25lf   rightSigma:   %29.25lf\n", leftSigma, rightSigma);
	Rprintf("leftweight:   %29.25lf   rightweight:  %29.25lf\n", leftWeight, rightWeight);
	Rprintf("I1left:       %29.25lf   I1right:      %29.25lf\n", I1left, I1right);
	Rprintf("I2left:       %29.25lf   I2right:      %29.25lf\n", I2left, I2right);
	Rprintf("leftMoment1:  %29.25lf   rightMoment1: %29.25lf\n", leftMoment1, rightMoment1);
	Rprintf("leftMoment2:  %29.25lf   rightMoment2: %29.25lf\n", leftMoment2, rightMoment2);
	Rprintf("moment1:      %29.25lf\n", moment1);
	Rprintf("moment2:      %29.25lf\n", moment2);
	Rprintf("a:            %29.25lf\n", a);
	Rprintf("b:            %29.25lf\n", b);
	Rprintf("1/sigmaz:     %29.25lf\n", 1/sigmaZ);*/
	
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
	double sigmaZ=(double) ((REAL(sigmaZS))[0]);
	double logSigmaZ=log(sigmaZ);


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
		//Rprintf("%d:\n", j);


		computeParameters(a, b, c, sigmaZ, logSigmaZ, logNormfact, 
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
	crossentropy=value*(-log(2*sigmaZ)-fabs(z)/sigmaZ);
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
	double sigmaZ=(double) ((REAL(sigmaZS))[0]);
	double logSigmaZ=log(sigmaZ);

	for(int j=0; j<elements; j++) {

		double a=(double)(REAL(aS)[j]);
		double b=(double)(REAL(bS)[j]);
		double c=(double)(REAL(cS)[j]);
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

		computeParameters(a, b, c, sigmaZ, logSigmaZ, logNormfact,
				Case,
				maxZ, logMaxValue,
				alpha, leftWeight, rightWeight, leftMy, rightMy, leftSigma,
				rightSigma,	moment1, moment2, entropy, crossentropy);

		if(method==1||method==2) {
			logNormfact=logNormfact+log(2.0);
		}
		double left, right, D, variable;
		quadratic(a, b+1.0/sigmaZ, c-logBorderFactor-logMaxValue+logNormfact, D, left, variable);
		quadratic(a, b-1.0/sigmaZ, c-logBorderFactor-logMaxValue+logNormfact, D, variable, right);

		if(method==1) {
			left=0.0;
		}
		else if(method==2) {
			right=0.0;
		}


		double logFmax=evalLogUnnormalizedPosterior(a, b, c, sigmaZ, logNormfact, maxZ);
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
