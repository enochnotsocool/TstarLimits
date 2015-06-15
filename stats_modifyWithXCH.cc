#include <iostream>
#include <cmath>
#include <cassert>
#include <sstream>

#include <TGraph.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TRandom3.h>

#include "binneddata.hh"
#include "fit.hh"
#include "statistics.hh"

using namespace std;

//##################################################################################################################################
// User Section 1
//
// (Change the code outside the User Sections only if you know what you are doing)
//

////////////////////////////////////////////////////////////////////////////////
// magic numbers
////////////////////////////////////////////////////////////////////////////////
// use Markov chain Monte Carlo (MCMC) to marginalize nuisance parameters
int useMCMC = 0;
// IMPORTANT: With useMCMC = 1, the systematic uncertanties are included in the limit calculation by default. Use the PAR_NUIS[] array below to control what uncertainties are included

// set the factor that defines the upper bound for the signal xs used by the MCMC as xsUpperBoundFactor*stat-only_limit
const double xsUpperBoundFactor=3.0;

// number of pseudoexperiments (when greater than 0, expected limit with +/- 1 and 2 sigma bands is calculated)
//const int NPES=200; // 200 (the more pseudo-experiments, the better. However, 200 is a reasonable choice)
const int NPES=5; 
//const int NPES=300; 
//const int NPES=400; 

// number of samples of nuisance parameters for Bayesian MC integration (when greater than 0, systematic uncertanties are included in the limit calculation)
//const int NSAMPLES=5000; // 4000 (larger value is better but it also slows down the code. 4000 is a reasonable compromise between the speed and precision)
int NSAMPLES=1000;
//const int NSAMPLES=10000;
//const int NSAMPLES=800000;

// alpha (1-alpha=confidence interval)
const double ALPHA=0.05;

// left side tail
const double LEFTSIDETAIL=0.0;
//const double LEFTSIDETAIL=0.1;

// output file name
//const string OUTPUTFILE="results/.root";
const string OUTPUTFILE=".root";

// center-of-mass energy
//const double sqrtS = 7000.;

// histogram binning
const int NBINS=9+3;    // Tgluon
const int NBINStr=8+2;    // Tgamma

double BOUNDARIES[NBINS+1] = {  350, 450, 550, 650, 750, 850, 950, 1050, 1150, 1250, 1350, 1450, 1550}; // default
double BOUNDARIEStr[NBINStr+1] = { 50, 200,  350, 500, 650, 800, 950, 1100, 1250, 1400, 1550}; // try to combine with Tr

// parameters
double SIGMASS=0;
const int NTgammaBkgShapeUnc = 8-1-1;
const int NPARS=18+1+NTgammaBkgShapeUnc;
const int NBKGPARS=6;
const int POIINDEX=0; // which parameter is "of interest"

//input the event yield after selection directly
double yields_af_selection = 0;
double yields_af_selection2 = 0;
double yields_af_selection3 = 0;
double BKGyields_af_selection3 = 0;
double Nsigma_tg = 15.;
double Nsigma_bkg_shape = 1.;

//const char* PAR_NAMES[NPARS]    = { "xs", "lumi", "sys", "stat", "sys2", "stat2", "p0", "p1", "p2", "p0_2", "p1_2", "p2_2", "n0", "n1", "n2", "n0_2", "n1_2", "n2_2","sys3","err_bkg3"};
const char* PAR_NAMES[NPARS]    = { "xs", "lumi", "sys", "stat", "sys2", "stat2", "p0", "p1", "p2", "p0_2", "p1_2", "p2_2", "n0", "n1", "n2", "n0_2", "n1_2", "n2_2","sys3","nbkg3_0","nbkg3_1","nbkg3_2","nbkg3_3","nbkg3_4","nbkg3_5"};//,"nbkg3_6"};

      double PAR_GUESSES[NPARS] = {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0,   1.00376e+04, -2.51828e+02, 1.12311e+02, 5.93010e+02, 8.04015e+01, 1.06339e+02,  0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0};//, 0};

      double PAR_MIN[NPARS]     = { 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0, -9999, -9999, 0, -9999, -9999, -1.*Nsigma_tg, -1.*Nsigma_tg, -1.*Nsigma_tg, -1.*Nsigma_tg, -1.*Nsigma_tg, -1.*Nsigma_tg,0.5, -1.*Nsigma_bkg_shape, -1.*Nsigma_bkg_shape, -1.*Nsigma_bkg_shape, -1.*Nsigma_bkg_shape, -1.*Nsigma_bkg_shape, -1.*Nsigma_bkg_shape};//, -1.*Nsigma_bkg_shape};

      double PAR_MAX[NPARS]     = {  300,  2.0, 2.0, 2.0, 2.0, 2.0, 999999, 9999, 9999, 999999, 9999, 9999, Nsigma_tg, Nsigma_tg, Nsigma_tg, Nsigma_tg, Nsigma_tg, Nsigma_tg,2.0, Nsigma_bkg_shape, Nsigma_bkg_shape, Nsigma_bkg_shape, Nsigma_bkg_shape, Nsigma_bkg_shape, Nsigma_bkg_shape};//, Nsigma_bkg_shape};

// Error is very improtant to TMinute, and they should be non-zero. 
// Otherwise, the number of parameters returned by TMinute will be wrong.
// (http://root.cern.ch/root/html/TMinuit.html#TMinuit:GetNumPars)
      double PAR_ERR[NPARS]     = { 1.0e-05,   0.044, 0.10,  0.10, 0.10, 0.10,   1.55022e+04, 1.05142e+01, 1.62065e+00, 3.56977e+02, 2.09964e+01, 1.53236e+00, 1, 1, 1, 1, 1, 1, 0.11, 1, 1, 1, 1, 1, 1};//, 1};

// 1,2 = signal (2 not used in the fit); 0,3 = background (3 not used in the fit)
      int PAR_TYPE[NPARS]       = {    1,      2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3};//, 3}; 

// 0 = not varied, >=1 = nuisance parameters with different priors (1 = Lognormal, 2 = Gaussian, 3 = Gamma, >=4 = Uniform)
      int PAR_NUIS[NPARS]       = {    0,      0,  1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 1, 2, 2, 2, 2, 2, 2};//, 4}; // default
      //int PAR_NUIS[NPARS]       = {    0,      0,  1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4};//, 4}; // default

// [Signal Unc for EXO limit]
int Nmass = 16;
double SignalTotalUnc[48][3] = {
    {500,          10.8, 5574.091797}, // Tgamma500GeV
    {550,          10.7, 2573.889404}, // Tgamma550GeV
    {600,          10.8, 1226.346802}, // Tgamma600GeV
    {650,          10.7, 632.710022}, // Tgamma650GeV
    {700,          10.8, 317.526978}, // Tgamma700GeV
    {750,          10.7, 174.833603}, // Tgamma750GeV
    {800,          10.8, 94.925888}, // Tgamma800GeV
    {850,          11.1, 53.797935}, // Tgamma850GeV
    {900,          11.0, 30.495079}, // Tgamma900GeV
    {950,          10.8, 17.935661}, // Tgamma950GeV
    {1000,          10.9, 9.918372}, // Tgamma1000GeV
    {1100,          11.1, 3.615197}, // Tgamma1100GeV
    {1200,          11.1, 1.293794}, // Tgamma1200GeV
    {1300,          11.1, 0.448085}, // Tgamma1300GeV
    {1400,          11.0, 0.162459}, // Tgamma1400GeV
    {1500,          11.3, 0.052235}, // Tgamma1500GeV
    {500,           0.0, 0.000000}, // TGluon500GeV
    {550,          79.5, 0.891844}, // TGluon550GeV
    {600,           0.0, 0.000000}, // TGluon600GeV
    {650,          72.4, 0.286988}, // TGluon650GeV
    {700,           0.0, 0.000000}, // TGluon700GeV
    {750,           0.0, 0.000000}, // TGluon750GeV
    {800,         103.3, 0.027793}, // TGluon800GeV
    {850,          59.8, 0.037518}, // TGluon850GeV
    {900,          74.2, 0.013288}, // TGluon900GeV
    {950,         102.2, 0.004644}, // TGluon950GeV
    {1000,           0.0, 0.000000}, // TGluon1000GeV
    {1100,          73.7, 0.001885}, // TGluon1100GeV
    {1200,          74.2, 0.000788}, // TGluon1200GeV
    {1300,           0.0, 0.000000}, // TGluon1300GeV
    {1400,           0.0, 0.000000}, // TGluon1400GeV
    {1500,          74.7, 0.000063}, // TGluon1500GeV
    {500,          10.4, 118.563499}, // TGluonTgamma500GeV
    {550,          11.0, 48.464127}, // TGluonTgamma550GeV
    {600,          11.8, 28.309416}, // TGluonTgamma600GeV
    {650,          11.1, 13.491109}, // TGluonTgamma650GeV
    {700,          10.4, 8.132792}, // TGluonTgamma700GeV
    {750,           9.7, 4.369513}, // TGluonTgamma750GeV
    {800,           9.7, 2.551290}, // TGluonTgamma800GeV
    {850,           9.1, 1.448032}, // TGluonTgamma850GeV
    {900,           7.7, 0.874790}, // TGluonTgamma900GeV
    {950,           9.8, 0.502066}, // TGluonTgamma950GeV
    {1000,           9.4, 0.294752}, // TGluonTgamma1000GeV
    {1100,           9.3, 0.120106}, // TGluonTgamma1100GeV
    {1200,           9.1, 0.041918}, // TGluonTgamma1200GeV
    {1300,           9.5, 0.013884}, // TGluonTgamma1300GeV
    {1400,          10.2, 0.005703}, // TGluonTgamma1400GeV
    {1500,           9.8, 0.002582}}; // TGluonTgamma1500GeV


string TgammaBkgShapeUnc[NTgammaBkgShapeUnc]={
    "Mergecate00",
    "Mergecate01",
    "Mergecate02",
    "Mergecate03",
    "SignalAndControlRegionscate0",
    //"SignalAndControlRegionscate3",
    "SignalAndControlRegionscate4"//,
    //"DataMCComparisoncate0",
    //"DataMCComparisoncate3",
    //"DataMCComparisoncate4"
};
TRandom3 *rnd3;
RandomPrior *MorphingParameter;
TH1D *histTgammaBkg;
TH1D *histTgammaBkgShapeUncPlus[NTgammaBkgShapeUnc];
TH1D *histTgammaBkgShapeUncMinus[NTgammaBkgShapeUnc];
bool UsingHiggsInsteadOfTheat = true;
//
// End of User Section 1
//##################################################################################################################################

// input files vector
vector<string> INPUTFILES;
vector<string> INPUTFILES2;
vector<string> INPUTFILES3;

// covariance matrix
double COV_MATRIX[NPARS][NPARS];
TMatrixDSym covMatrix = TMatrixDSym(NBKGPARS);
TVectorD eigenValues = TVectorD(NBKGPARS);
TMatrixD eigenVectors = TMatrixD(NBKGPARS,NBKGPARS);

// constrain S to be positive in the S+B fit
//const bool posS = 0;
const bool posS = 1;

// use B-only fit with fixed but non-zero signal when calculating the covariance matrix used for background systematics
const bool BonlyFitForSyst = 1;

// shift in the counter used to extract the covariance matrix
int shift = 2;

// branching fraction 
double BR = 1.;
double BranchRatioForTgluon = 1.;

// resonance shape type
string ResShapeType = "Combined";

// T+gluon
TH1D* HISTCDF=0; // signal CDF
TH1D* HISTCDF2=0; // signal CDF
// T+gamma
TH1D* HISTCDF3=0; // signal CDF
TH1D* HISTCDFBKG3=0; // BKG CDF

////////////////////////////////////////////////////////////////////////////////
// function integral
////////////////////////////////////////////////////////////////////////////////
vector <double> INTEGRAL(double *x0, double *xf, double *par)
//std::pair <double,double> INTEGRAL(double *x0, double *xf, double *par)
//double INTEGRAL(double *x0, double *xf, double *par)
{
  if(DEBUGMODE)
      printf("[INTEGRAL] xs : %f p(2~20:%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)\n",par[0], 
              par[1], par[2], par[3], 
              par[4], par[5], par[6], par[7],
              par[8], par[9], par[10], par[11],
              par[12], par[13], par[14], par[15],
              par[16], par[17], par[18], par[19]);
  double xs=par[0];
  double lumi=par[1];
  double sys=par[2];
  double stat=par[3];
  double sys2=par[4];
  double stat2=par[5];
  double p0=par[6];
  double p1=par[7];
  double p2=par[8];
  double p0_2=par[9];
  double p1_2=par[10];
  double p2_2=par[11];
  double n[NBKGPARS] = {0.};
  n[0]=par[12];
  n[1]=par[13];
  n[2]=par[14];
  n[3]=par[15];
  n[4]=par[16];
  n[5]=par[17];
  double sys3 = par[18];
  //double err_bkg3 = par[19];    // considered for background shape uncertainties for Tgamma channel
  double nbkg3[NTgammaBkgShapeUnc];
  for(int i=0;i<NTgammaBkgShapeUnc;i++)
      nbkg3[i] = par[19+i];

  if( COV_MATRIX[0+shift][0+shift]>0. && (n[0]!=0. || n[1]!=0. || n[2]!=0. || n[3]!=0. || n[4]!=0. || n[5]!=0.) )
  {
    double g[NBKGPARS] = {0.};
    for(int v=0; v<NBKGPARS; ++v)
    {
      for(int k=0; k<NBKGPARS; ++k) g[k]=n[v]*eigenValues(v)*eigenVectors[k][v];
      p0 += g[0];
      p1 += g[1];
      p2 += g[2];
      p0_2 += g[3];
      p1_2 += g[4];
      p2_2 += g[5];
    }
  }

  // uses Simpson's 3/8th rule to compute the background integral over a short interval
  // also use a power series expansion to determine the intermediate intervals since the pow() call is expensive
  double dx=(xf[0]-x0[0])/3.;
  double x=x0[0];

  double a=p0/(1+exp((x-p1)/p2));
  double b=dx*(-p0)/p2*(exp((p1+x)/p2))/(exp(p1/p2)+exp(x/p2))/(exp(p1/p2)+exp(x/p2));
  double c=0.5*dx*dx*p0/p2/p2*exp((p1+x)/p2)*(exp(x/p2)-exp(p1/p2))/(exp(p1/p2)+exp(x/p2))/(exp(p1/p2)+exp(x/p2))/(exp(p1/p2)+exp(x/p2));
  double d=0.166666667*dx*dx*dx*(-1)*p0/p2/p2/p2*(exp((p1+x)/p2))*((-4)*(exp((p1+x)/p2))+exp(2*p1/p2)+exp(2*x/p2))/(exp(p1/p2)+exp(x/p2))/(exp(p1/p2)+exp(x/p2))/(exp(p1/p2)+exp(x/p2))/(exp(p1/p2)+exp(x/p2));

  double bkg=(xf[0]-x0[0])*(a+0.375*(b+c+d)+0.375*(2*b+4*c+8*d)+0.125*(3*b+9*c+27*d));//up to 3rd order
  //double bkg=(xf[0]-x0[0])*(a+0.375*(b+c)+0.375*(2*b+4*c)+0.125*(3*b+9*c));//up to 2nd order
  //double bkg=(xf[0]-x0[0])*(a+0.375*(b)+0.375*(2*b)+0.125*(3*b));//up to 1st order 
  //double x1 = (xf[0]+x0[0])/2;
  //double a1 = p0/(1+exp((x1-p1)/p2));
  //printf("x/a/a1/b/c/d/bkg = %f, %f, %f, %f, %f, %f, %f\n", x, a*50, a1*50, b, c, d, bkg);

  if(bkg<0.) bkg=0.0000001;

  int bin1=HISTCDF->GetXaxis()->FindBin(xf[0]);
  int bin2=HISTCDF->GetXaxis()->FindBin(x0[0]);
  if(bin1<1) bin1=1;
  if(bin1>HISTCDF->GetNbinsX()) bin1=HISTCDF->GetNbinsX();
  if(bin2<1) bin1=1;
  if(bin2>HISTCDF->GetNbinsX()) bin2=HISTCDF->GetNbinsX();
  //double sig=xs*lumi*(HISTCDF->GetBinContent(bin1)-HISTCDF->GetBinContent(bin2));

  //input signal event direct
  //printf("sys : %f, stat : %f\n", sys, stat);
  //double sig=sys*stat*yields_af_selection*xs*lumi*(HISTCDF->GetBinContent(bin1)-HISTCDF->GetBinContent(bin2));
  double sig=sys*yields_af_selection*xs*lumi*(HISTCDF->GetBinContent(bin1)-HISTCDF->GetBinContent(bin2));
  //printf("yields %f\n", yields_af_selection);

  //second data
  a=p0_2/(1+exp((x-p1_2)/p2_2));
  b=dx*(-p0_2)/p2_2*(exp((p1_2+x)/p2_2))/(exp(p1_2/p2_2)+exp(x/p2_2))/(exp(p1_2/p2_2)+exp(x/p2_2));
  c=0.5*dx*dx*p0_2/p2_2/p2_2*exp((p1_2+x)/p2_2)*(exp(x/p2_2)-exp(p1_2/p2_2))/(exp(p1_2/p2_2)+exp(x/p2_2))/(exp(p1_2/p2_2)+exp(x/p2_2))/(exp(p1_2/p2_2)+exp(x/p2_2));
  d=0.166666667*dx*dx*dx*(-1)*p0_2/p2_2/p2_2/p2_2*(exp((p1_2+x)/p2_2))*((-4)*(exp((p1_2+x)/p2_2))+exp(2*p1_2/p2_2)+exp(2*x/p2_2))/(exp(p1_2/p2_2)+exp(x/p2_2))/(exp(p1_2/p2_2)+exp(x/p2_2))/(exp(p1_2/p2_2)+exp(x/p2_2))/(exp(p1_2/p2_2)+exp(x/p2_2));
  double bkg2=(xf[0]-x0[0])*(a+0.375*(b+c+d)+0.375*(2*b+4*c+8*d)+0.125*(3*b+9*c+27*d));//up to 3rd order

  bin1=HISTCDF2->GetXaxis()->FindBin(xf[0]);
  bin2=HISTCDF2->GetXaxis()->FindBin(x0[0]);
  if(bin1<1) bin1=1;
  if(bin1>HISTCDF2->GetNbinsX()) bin1=HISTCDF2->GetNbinsX();
  if(bin2<1) bin1=1;
  if(bin2>HISTCDF2->GetNbinsX()) bin2=HISTCDF2->GetNbinsX();
  double sig2=sys2*yields_af_selection2*xs*lumi*(HISTCDF2->GetBinContent(bin1)-HISTCDF2->GetBinContent(bin2));
  if(bkg2<0.) bkg2=0.0000001;

  bin1=HISTCDFBKG3->GetXaxis()->FindBin(xf[0]);
  bin2=HISTCDFBKG3->GetXaxis()->FindBin(x0[0]);
  if(bin1<1) bin1=1;
  if(bin1>HISTCDFBKG3->GetNbinsX()) bin1=HISTCDFBKG3->GetNbinsX();
  if(bin2<1) bin1=1;
  if(bin2>HISTCDFBKG3->GetNbinsX()) bin2=HISTCDFBKG3->GetNbinsX();

  double err_bkg3 = 1.;
  double IsNormal = 0.;
  for(int i=0;i<NTgammaBkgShapeUnc;i++){
      IsNormal += nbkg3[i];
  }
  //if(par[19]==1. /*Normal*/){
  //    err_bkg3 = 1.;
  //}else{
  if(IsNormal!=0)
      if(UsingHiggsInsteadOfTheat){
          err_bkg3 = -1.;
              double lamda_= 0.;
              int bin_ = histTgammaBkg->GetXaxis()->FindBin((x0[0]+xf[0])/2);
              err_bkg3 = histTgammaBkg->GetBinContent(bin_);
              double a_lamda = 0.;
              double b_lamda = 0.;
              double c_lamda = 0.;
              for(int i=0;i<NTgammaBkgShapeUnc;i++){
                  //lamda_ = rnd3->Gaus(0, 1);
                  //lamda_ = exp(rnd3->Gaus(0, 1));   // log-normal (corresponding to 1.0 width and 1.0 mean)
                  //lamda_ = MorphingParameter->getRandom();   // log-normal using statistics
                  //lamda_ = rnd3->Uniform(-10, 10);
                  lamda_ = nbkg3[i];

                  a_lamda = 0.;
                  b_lamda = 0.;
                  c_lamda = 0.;
                  if(fabs(lamda_)>=1){      // fixed a bug on 110914
                      a_lamda =  max(lamda_, 0.);
                      b_lamda = -1.*fabs(lamda_);
                      c_lamda =  max(-1.*lamda_, 0.);
                  }else{
                      a_lamda =  lamda_*(lamda_+1.)/2.;
                      b_lamda =  lamda_*(lamda_-1.)/2.;
                      c_lamda =  -1.*lamda_*lamda_;
                  }
                  err_bkg3 += a_lamda*histTgammaBkgShapeUncPlus[i]->GetBinContent(bin_) + 
                      b_lamda*histTgammaBkg->GetBinContent(bin_)+
                      c_lamda*histTgammaBkgShapeUncMinus[i]->GetBinContent(bin_);
              }
              err_bkg3 = err_bkg3/histTgammaBkg->GetBinContent(bin_);
      }else{
          err_bkg3 = 1.; 
          double rnd_= 0.; 
          int bin_ = histTgammaBkg->GetXaxis()->FindBin((x0[0]+xf[0])/2);
          for(int i=0;i<NTgammaBkgShapeUnc;i++){
              //rnd_ = rnd3->Gaus(0, 10);
              //rnd_ = rnd3->Gaus(0, 1);
              //rnd_ = rnd3->Uniform(-10,10);
              rnd_ = nbkg3[i];
              if(rnd_>=0){
                  err_bkg3*= pow(histTgammaBkgShapeUncPlus[i]->GetBinContent(bin_)/histTgammaBkg->GetBinContent(bin_), rnd_);
              }else{
                  err_bkg3*= pow(histTgammaBkgShapeUncMinus[i]->GetBinContent(bin_)/histTgammaBkg->GetBinContent(bin_), -1.*rnd_);
              }   
          } 
      }
  //}
  if(DEBUGMODE)
      printf("[INTEGRAL] err_bkg3 (%f) and par[19] (%f)\n ",err_bkg3,par[19]);
  double bkg3=err_bkg3*BKGyields_af_selection3*(HISTCDFBKG3->GetBinContent(bin1)-HISTCDFBKG3->GetBinContent(bin2));

  bin1=HISTCDF3->GetXaxis()->FindBin(xf[0]);
  bin2=HISTCDF3->GetXaxis()->FindBin(x0[0]);
  if(bin1<1) bin1=1;
  if(bin1>HISTCDF3->GetNbinsX()) bin1=HISTCDF3->GetNbinsX();
  if(bin2<1) bin1=1;
  if(bin2>HISTCDF3->GetNbinsX()) bin2=HISTCDF3->GetNbinsX();
  double sig3=sys3*yields_af_selection3*xs*(HISTCDF3->GetBinContent(bin1)-HISTCDF3->GetBinContent(bin2));
  if(bkg3<0.) bkg3=0.0000001;

  vector <double> _estimate;
  if(BranchRatioForTgluon==0){
      _estimate.push_back(0.);
      _estimate.push_back(0.);
  }else{
      _estimate.push_back(bkg+sig);
      _estimate.push_back(bkg2+sig2);
  }
  _estimate.push_back(bkg3+sig3);

  if(bkg+sig < 0.) { printf("Negative estimated yields in muon %f\n", bkg+sig);_estimate.at(0) = 1e-10;}
  if(bkg2+sig2 < 0.) { printf("Negative estimated yields in electron %f\n", bkg2+sig2);_estimate.at(1) = 1e-10;}
  if(bkg3+sig3 < 0.) { printf("Negative estimated yields in Tgamma %f\n", bkg3+sig3);_estimate.at(2) = 1e-10;}

//printf("test %f  %f %f %f\n", bkg, bkg2, sig, sig2);
  return _estimate;
}

double GetCombinedSysInTGTGAnalysis(TFile* histfile, int ithBin);
////////////////////////////////////////////////////////////////////////////////
// main function
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{

  if(argc<=1) {
    cout << "Usage: stats signalmass [BR ResShapeType NSAMPLES AdditionRNDSeed useMCMC]" << endl;
    cout << "Example: stats 700 1 Combined 1000 0 1" << endl;
    return 0;
  }

  SIGMASS = atof(argv[1]);
  int masspoint = int(SIGMASS);
  int AdditionRNDSeed = 0;
  //if(argc>2) BR = atof(argv[2]);
  if(argc>2) BranchRatioForTgluon = atof(argv[2]);
  if(argc>3) ResShapeType = argv[3];
  if(argc>4) NSAMPLES = atoi(argv[4]);
  if(argc>5) AdditionRNDSeed = atoi(argv[5]);
  if(argc>6) useMCMC = atoi(argv[6]);
  double BranchRatioForTgamma = 1.0 - BranchRatioForTgluon;
  char buffer[128];

  if(masspoint<=850){
      PAR_MAX[0] = 10;
  }else if(masspoint>850&&masspoint<=1000){
      PAR_MAX[0] = 50;
  }else if(masspoint>1000){
      PAR_MAX[0] = 300;
  }

  for(int i=0;i<Nmass;i++){
      if(masspoint==SignalTotalUnc[i][0]){
          /*
                Tgamma : 0~15
                TGluon : (0~15) + 16 [ Ignore this due to very limited statistics, and contribution in BR(10%) is 0.01 of TgTg ]
                TGluonTgamma : 0~15 + 16*2
           */
          double Y_TgTg   = BranchRatioForTgamma*BranchRatioForTgamma*SignalTotalUnc[i][2];
          double Err_TgTg = (SignalTotalUnc[i][1]/100.)*Y_TgTg;
          double Y_TGTg   = 2*BranchRatioForTgluon*BranchRatioForTgamma*SignalTotalUnc[i+2*Nmass][2];
          double Err_TGTg = (SignalTotalUnc[i+2*Nmass][1]/100.)*Y_TGTg;
          PAR_ERR[NPARS-2] = sqrt(Err_TgTg*Err_TgTg+Err_TGTg*Err_TGTg)/(Y_TgTg+Y_TGTg);
          break;
      }
  }

  MorphingParameter = new RandomPrior(1/*logNormal*/, 0.001/*mean*/, 1./*width*/, 0./*min*/, 100/*max*/);
  //##################################################################################################################################
  // User Section 2
  //
  // (Change the code outside the User Sections only if you know what you are doing)
  //

  // input data file
  ostringstream inputfile;
  inputfile << "ExcitedQuarkAnalysis_Mu" <<".root";
  INPUTFILES.push_back(inputfile.str());

  ostringstream inputfile2;
  inputfile2 << "ExcitedQuarkAnalysis_Ele" <<".root";
  INPUTFILES2.push_back(inputfile2.str());

  // data histogram name
  string datahistname = "data_obs";

  // input signal files with resonance shapes
  string filename1 = inputfile.str();
  string filename2 = inputfile.str();
  string filename1_2 = inputfile2.str();
  string filename2_2 = inputfile2.str();

  sprintf(buffer,"data/Tgamma/output_theta_%1.1f.root",(float)BranchRatioForTgamma);
  if(BranchRatioForTgamma==0)
      sprintf(buffer,"data/Tgamma/output_theta_1.0.root"); // will be set as zero later
  string filename1_3(buffer);
  string filename2_3(buffer);
  INPUTFILES3.push_back(filename1_3);

  // signal histogram names
  // combined shape : hSigMC_{masspoint}_{BranchRatioForTgluon}
  sprintf(buffer,"_%1.1f",BranchRatioForTgluon);
  ostringstream histname1, histname2, histname3, histnameBKG3, datahistname3;
  histname1 << "hSigMC_"<< masspoint<<buffer;
  histname2 << "hSigMC_"<< masspoint<<buffer;
  histname3 << "Tgamma"<< masspoint<<"GeV__signal"<<masspoint;
  histnameBKG3 << "Tgamma"<< masspoint<<"GeV__pdf_bkg";
  datahistname3 << "Tgamma"<< masspoint<<"GeV__DATA";

  // loading bkg shape uncertainties for Tgamma channel
  rnd3 = new TRandom3(31415+AdditionRNDSeed);
  TFile *fileTgamma = new TFile(filename1_3.c_str());
  sprintf(buffer,"Tgamma%iGeV__pdf_bkg",masspoint);
  histTgammaBkg = (TH1D*) fileTgamma->Get(buffer); 
  for(int i=0;i<NTgammaBkgShapeUnc;i++){
      sprintf(buffer,"Tgamma%iGeV__pdf_bkg__%s__plus",masspoint,TgammaBkgShapeUnc[i].c_str());
      histTgammaBkgShapeUncPlus[i] = (TH1D*) fileTgamma->Get(buffer);
      sprintf(buffer,"Tgamma%iGeV__pdf_bkg__%s__minus",masspoint,TgammaBkgShapeUnc[i].c_str());
      histTgammaBkgShapeUncMinus[i] = (TH1D*) fileTgamma->Get(buffer);
  }

  // loading signal uncertainties
  TFile* histfile=new TFile(filename1.c_str());
  PAR_ERR[2] = GetCombinedSysInTGTGAnalysis(histfile,1);
  PAR_ERR[3] = GetCombinedSysInTGTGAnalysis(histfile,2);

  TFile* histfile2=new TFile(filename1_2.c_str());
  PAR_ERR[4] = GetCombinedSysInTGTGAnalysis(histfile2,1);
  PAR_ERR[5] = GetCombinedSysInTGTGAnalysis(histfile2,2);

  //New version already include stat in sys
  //PAR_ERR[2] = sqrt(pow(PAR_ERR[2], 2)+pow(PAR_ERR[3], 2));
  //PAR_ERR[4] = sqrt(pow(PAR_ERR[4], 2)+pow(PAR_ERR[5], 2));
  printf("signal err : %f(muon), %f(electron)\n", PAR_ERR[2], PAR_ERR[4]);

  // Information for theta package
  TH1D *hSigNormalMu  = (TH1D*)histfile->Get(histname1.str().c_str());
  TH1D *hSigNormalEle = (TH1D*)histfile2->Get(histname2.str().c_str());

  // End of User Section 2
  //##################################################################################################################################

  // PAR VALUES and COVARIANCE
  TH1D *hPARs[2];	// mu / ele
  sprintf(buffer,"hPARs_Mu_%i",masspoint);
  hPARs[0] = new TH1D(buffer,"",13,0,13);	// 3parameters + 3x3 covariance + [signal sys]
  sprintf(buffer,"hPARs_Ele_%i",masspoint);
  hPARs[1] = new TH1D(buffer,"",13,0,13);	// 3parameters + 3x3 covariance + [signal sys]

  //PAR Distribution
  TH1F* nui_dis[7];
  nui_dis[0] = new TH1F("mu_p1", "mu_p1", 500, 0, fabs(PAR_GUESSES[0+6])*2);
  nui_dis[1] = new TH1F("mu_p2", "mu_p2", 500, -fabs(PAR_GUESSES[1+6])*2.5, fabs(PAR_GUESSES[1+6])*2.5);
  nui_dis[2] = new TH1F("mu_p3", "mu_p3", 500, fabs(PAR_GUESSES[2+6])*0.5, fabs(PAR_GUESSES[2+6])*1.5);
  nui_dis[3] = new TH1F("ele_p1", "ele_p1", 500, 0, fabs(PAR_GUESSES[3+6])*10);
  nui_dis[4] = new TH1F("ele_p2", "ele_p2", 500, -fabs(PAR_GUESSES[4+6])*40, fabs(PAR_GUESSES[4+6])*40);
  nui_dis[5] = new TH1F("ele_p3", "ele_p3", 500, fabs(PAR_GUESSES[5+6])*0.5, fabs(PAR_GUESSES[5+6])*1.5);
  nui_dis[6] = new TH1F("sig", "sig", 500, fabs(PAR_GUESSES[0])*0, 5);
  FILE * paraFile;
  ostringstream para_file_name;
  //para_file_name << "parameter_h/com_para_" << masspoint <<"_"<< BranchRatioForTgluon << ".h";
  sprintf(buffer,"parameter_h/com_para_%i_%1.1f_%i.h",masspoint,BranchRatioForTgluon,AdditionRNDSeed);
  //paraFile = fopen (para_file_name.str().c_str(),"w");
  paraFile = fopen (buffer,"w");

  FILE * limitFile;
  ostringstream limit_file_name;
  //limit_file_name << "LOG/com_limit_" << masspoint <<"_"<< BranchRatioForTgluon << ".log";
  sprintf(buffer,"LOG/com_limit_%i_%1.1f_%i.log",masspoint,BranchRatioForTgluon,AdditionRNDSeed);
  //limitFile = fopen (limit_file_name.str().c_str(),"w");
  limitFile = fopen (buffer,"w");


  if(BonlyFitForSyst) shift = 0;

  if(!posS) PAR_MIN[POIINDEX] = -PAR_MAX[POIINDEX];

  // initialize the covariance matrix
  for(int i = 0; i<NPARS; ++i) { for(int j = 0; j<NPARS; ++j) COV_MATRIX[i][j]=0.; }

  HISTCDF=getSignalCDF(filename1.c_str(), histname1.str().c_str(), filename2.c_str(), histname2.str().c_str(), BR, 1., 1., yields_af_selection,0,"_TgM");
  HISTCDF2=getSignalCDF(filename1_2.c_str(), histname2.str().c_str(), filename2_2.c_str(), histname2.str().c_str(), BR, 1., 1., yields_af_selection2,0,"_TgE");

  HISTCDF3=getSignalCDF(filename1_3.c_str(), histname3.str().c_str(), filename2_3.c_str(), histname3.str().c_str(), BR, 1., 1., yields_af_selection3,1,"_Tr");

  HISTCDFBKG3=getSignalCDF(filename1_3.c_str(), histnameBKG3.str().c_str(), filename2_3.c_str(), histnameBKG3.str().c_str(), BR, 1., 1., BKGyields_af_selection3,1,"_TrBkg");
  if(true){
      //HISTCDF->Scale(BranchRatioForTgluon*BranchRatioForTgluon);
      //HISTCDF2->Scale(BranchRatioForTgluon*BranchRatioForTgluon);
      //yields_af_selection *= BranchRatioForTgluon*BranchRatioForTgluon;
      //yields_af_selection2 *= BranchRatioForTgluon*BranchRatioForTgluon;

      //HISTCDF3->Scale((1.-BranchRatioForTgluon)*(1.-BranchRatioForTgluon));
      //yields_af_selection3 *= (1.-BranchRatioForTgluon)*(1.-BranchRatioForTgluon);    // taken into account in template fit
      if(BranchRatioForTgluon==1)
          yields_af_selection3 = 0.;    // set TgTg in BR(0%) as 0
  }


  printf("Event Yields %f(muon), %f(electron)\nTgamma: Signal(%f) + Bkg(%f)\n", yields_af_selection, yields_af_selection2,
          yields_af_selection3,
          BKGyields_af_selection3);

  //assert(HISTCDF && SIGMASS>0);
  assert(((HISTCDF&&HISTCDF2)||(HISTCDFBKG3&&HISTCDF3)) && SIGMASS>0);

  // get the data
  std::string _hdataname1 = "normalized_data1";
  std::string _hdataname2 = "normalized_data2";
  std::string _hdataname3 = "normalized_data3";
  TH1D* data=getData(INPUTFILES, datahistname.c_str(), NBINS, BOUNDARIES, _hdataname1,0);
  TH1D* data2=getData(INPUTFILES2, datahistname.c_str(), NBINS, BOUNDARIES, _hdataname2,0);
  TH1D* data3=getData(INPUTFILES3, datahistname3.str().c_str(), NBINStr, BOUNDARIEStr, _hdataname3,1);

  /*Do a test*/
  /*
  data3->SetBinContent(1 ,2 /150.);
  data3->SetBinContent(2 ,22/150.);   // 22 : standard
  data3->SetBinContent(3 ,6 /150.);
  data3->SetBinContent(4 ,2 /150.);
  data3->SetBinContent(5 ,1 /150.);
  data3->SetBinContent(6 ,1 /150.);
  data3->SetBinContent(7 ,0 /150.);
  data3->SetBinContent(8 ,0 /150.);
  data3->SetBinContent(9 ,0 /150.);
  data3->SetBinContent(10,0 /150.);
  */

  if ( BranchRatioForTgluon == 1.0 ){
      // set Tgamma's data and Bkg as zero, of course signal should be zero as above
      for(int i=1;i<=data3->GetXaxis()->GetNbins();i++)
          data3->SetBinContent(i,0.);
      for(int i=1;i<=HISTCDFBKG3->GetXaxis()->GetNbins();i++)
          HISTCDFBKG3->SetBinContent(i,0.);
      BKGyields_af_selection3 = 0.;
      PAR_NUIS[18] = 0;
      for(int i=0;i<NTgammaBkgShapeUnc;i++)
          PAR_NUIS[19+i] = 0;
  }else if( BranchRatioForTgluon == 0.0 ){
      // set Tgluon's data as zero, of course signal should be zero as above
      for(int i=1;i<=data->GetXaxis()->GetNbins();i++)
          data->SetBinContent(i,0.);
      for(int i=1;i<=data2->GetXaxis()->GetNbins();i++)
          data2->SetBinContent(i,0.);
      // ignore fit parameters in a pure T+gamma channel
      for(int i = 0; i < NBKGPARS; i++) {
          PAR_TYPE[i+6] = 3;
          PAR_GUESSES[i+6] = 0;
          PAR_NUIS[i+12] = 0;
      }
      PAR_GUESSES[6] = 0;
      PAR_GUESSES[9] = 0;
      PAR_MAX[6] = 0;
      PAR_MAX[9] = 0;
      PAR_NUIS[2] = 0;
      PAR_NUIS[4] = 0;
  }
  if(PAR_ERR[2]<=0)    PAR_ERR[2] = 0.05;
  if(PAR_ERR[3]<=0)    PAR_ERR[3] = 0.05;
  if(PAR_ERR[4]<=0)    PAR_ERR[4] = 0.05;
  if(PAR_ERR[5]<=0)    PAR_ERR[5] = 0.05;
  printf("Data : Tgluon(%f, %f) + Tgamma(%f) \n",
          data->Integral(),
          data2->Integral(),
          data3->Integral());

  if(DEBUGMODE)
      for(int i=0; i<NPARS; i++) printf("[Parameter-%i] %s : %f, %f (%f -- %f) NUIS(%i)\n",
          i, PAR_NAMES[i], PAR_GUESSES[i], PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);
  // create the output file
  ostringstream outputfile;
  outputfile <<"nuisance/"<<  masspoint << "_" << BranchRatioForTgluon << "_" << ResShapeType << "_"<<AdditionRNDSeed<<".root";
  //outputfile <<"nuisance/"<< OUTPUTFILE.substr(0,OUTPUTFILE.find(".root")) << masspoint << "_" << BranchRatioForTgluon << "_" << ResShapeType << "_"<<AdditionRNDSeed<<".root";
  TFile* rootfile=new TFile(outputfile.str().c_str(), "RECREATE");  rootfile->cd();

  // POI value
  double POIval = 0;

  // setup an initial fitter to perform a signal+background fit
  printf("====================INITIAL Fit====================\n");
  Fitter *initfit = new Fitter(AdditionRNDSeed, data, data2, data3, INTEGRAL);
  initfit->setStrategy(2);//initial fit use SET STR 2 first
  for(int i=0; i<NPARS; i++) initfit->defineParameter(i, PAR_NAMES[i], PAR_GUESSES[i], PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);

  // do an initial signal+background fit first
  for(int i=0; i<NPARS; i++) if(PAR_TYPE[i]>=2 || PAR_MIN[i]==PAR_MAX[i]) initfit->fixParameter(i);
  initfit->doFit();
  POIval = initfit->getParameter(POIINDEX); // get the POI value for later use
  double POIval_err = initfit->getParError(POIINDEX);
  initfit->fixParameter(POIINDEX); // a parameter needs to be fixed before its value can be changed
  initfit->setParameter(POIINDEX, 0.0); // set the POI value to 0 to get the B component of the S+B fit (for calculating pulls and generating pseudo-data)
  initfit->setPrintLevel(0);
  int _mu = 1; int _ele = 2;
  initfit->calcPull("pull_bkg_init_mu", _mu)->Write();
  initfit->calcDiff("diff_bkg_init_mu", _mu)->Write();
  initfit->calcPull("pull_bkg_init_ele", _ele)->Write();
  initfit->calcDiff("diff_bkg_init_ele", _ele)->Write();
  initfit->write("fit_bkg_init");
  initfit->setParameter(POIINDEX, POIval);

  fprintf(paraFile, "double com_sig_%i = %e;\n", masspoint, POIval);
  fprintf(paraFile, "double com_mu_p1_%i = %e;\n", masspoint, initfit->getParameter(6));
  fprintf(paraFile, "double com_mu_p2_%i = %e;\n", masspoint, initfit->getParameter(7));
  fprintf(paraFile, "double com_mu_p3_%i = %e;\n", masspoint, initfit->getParameter(8));
  fprintf(paraFile, "double com_ele_p1_%i = %e;\n", masspoint, initfit->getParameter(9));
  fprintf(paraFile, "double com_ele_p2_%i = %e;\n", masspoint, initfit->getParameter(10));
  fprintf(paraFile, "double com_ele_p3_%i = %e;\n", masspoint, initfit->getParameter(11));
  fprintf(paraFile, "double com_sig_err_%i = %e;\n", masspoint, POIval_err);
  fprintf(paraFile, "double com_mu_p1_err_%i = %e;\n", masspoint, initfit->getParError(6));
  fprintf(paraFile, "double com_mu_p2_err_%i = %e;\n", masspoint, initfit->getParError(7));
  fprintf(paraFile, "double com_mu_p3_err_%i = %e;\n", masspoint, initfit->getParError(8));
  fprintf(paraFile, "double com_ele_p1_err_%i = %e;\n", masspoint, initfit->getParError(9));
  fprintf(paraFile, "double com_ele_p2_err_%i = %e;\n", masspoint, initfit->getParError(10));
  fprintf(paraFile, "double com_ele_p3_err_%i = %e;\n", masspoint, initfit->getParError(11));
  fclose(paraFile);

  // setup the limit values
  double observedLowerBound, observedUpperBound;
  vector<double> expectedLowerBounds;
  vector<double> expectedUpperBounds;

  cout << "*********** pe=0 (data) ***********" << endl;

  // setup the fitter with the input from the signal+background fit
  printf("====================SECOND Fit for observed====================\n");
  Fitter *fit_data = new Fitter(AdditionRNDSeed, data, data2, data3, INTEGRAL);
  fit_data->setStrategy(2);//observed fit use SET STR 2 also
  fit_data->setPOIIndex(POIINDEX);
  //fit_data->setPrintLevel(0);
  for(int i=0; i<NPARS; i++) fit_data->defineParameter(i, PAR_NAMES[i], initfit->getParameter(i), PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);

  // perform a signal+background fit possibly followed by a background-only fit with a fixed but non-zero signal
  for(int i=0; i<NPARS; i++) if(PAR_TYPE[i]>=2 || PAR_MIN[i]==PAR_MAX[i]) fit_data->fixParameter(i);
  if(BonlyFitForSyst) { fit_data->doFit();cout << "Data fit status: " << fit_data->getFitStatus() << endl; fit_data->fixParameter(POIINDEX); }
  fit_data->doFit(&COV_MATRIX[0][0], NPARS);
  cout << "Data fit status: " << fit_data->getFitStatus() << endl;
  POIval = fit_data->getParameter(POIINDEX); // get the POI value for later use
  fit_data->fixParameter(POIINDEX); // a parameter needs to be fixed before its value can be changed
  fit_data->setParameter(POIINDEX, 0.0); // set the POI value to 0 to get the B component of the S+B fit (for calculating pulls and generating pseudo-data)
  fit_data->setPrintLevel(0);
  fit_data->calcPull("pull_bkg_0_mu", _mu)->Write();
  fit_data->calcDiff("diff_bkg_0_mu", _mu)->Write();
  fit_data->calcPull("pull_bkg_0_ele", _ele)->Write();
  fit_data->calcDiff("diff_bkg_0_ele", _ele)->Write();
  fit_data->write("fit_bkg_0");

  //Fill int nuisance distribution
  for(int i = 0; i < 6; i++){
      nui_dis[i]->Fill(fit_data->getParameter(6+i));
  }
  nui_dis[6]->Fill(POIval);

  // calculate eigenvalues and eigenvectors
  for(int i = 0; i<NBKGPARS; ++i) { for(int j = 0; j<NBKGPARS; ++j) { covMatrix(i,j)=COV_MATRIX[i+shift][j+shift]; } }
  printf("[JACKY TEST] COV_MATRIX start\n");
  covMatrix.Print();

  if( BranchRatioForTgluon == 0.0 ){
      // ignore fit parameters in a pure T+gamma channel
      for(int i = 0; i<NBKGPARS; ++i) { for(int j = 0; j<NBKGPARS; ++j) { covMatrix(i,j)=0.0; } }   
  }

  hPARs[0]->Fill(0.,fit_data->getParameter(6));
  hPARs[0]->Fill(1.,fit_data->getParameter(7));
  hPARs[0]->Fill(2.,fit_data->getParameter(8));
  int counter=0;
  for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
          hPARs[0]->Fill(3.+counter, covMatrix(i,j) );
          counter++;
      }
  }
  //hPARs[0]->Fill(12.,(double)data->Integral());
  //hPARs[0]->Fill(13.,(double)hSigNormalMu->Integral());
  hPARs[0]->Fill(12.,(double)PAR_ERR[2]);

  //printf("(%1.4f, %1.4f, %1.4f) \n",data->Integral(), hSigNormalMu->Integral(), SigSysNumber->GetBinContent(1));
  for(int ihpars=1;ihpars<=hPARs[0]->GetNbinsX();ihpars++)
      printf("        %10.4f\n", hPARs[0]->GetBinContent(ihpars));

  hPARs[1]->Fill(0.,fit_data->getParameter(6+3));
  hPARs[1]->Fill(1.,fit_data->getParameter(7+3));
  hPARs[1]->Fill(2.,fit_data->getParameter(8+3));
  counter=0;
  for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
          hPARs[1]->Fill(3.+counter, covMatrix(i+3,j+3) );
          counter++;
      }
  }
  //hPARs[1]->Fill(12.,(double)data2->Integral());
  //hPARs[1]->Fill(13.,(double)hSigNormalEle->Integral());
  hPARs[1]->Fill(12.,(double)PAR_ERR[4]);

  //printf("(%1.4f, %1.4f, %1.4f) \n",data2->Integral(), hSigNormalEle->Integral(), SigSysNumber2->GetBinContent(1));
  for(int ihpars=1;ihpars<=hPARs[1]->GetNbinsX();ihpars++)
      printf("        %10.4f\n", hPARs[1]->GetBinContent(ihpars));

  hPARs[0]->SaveAs(TString::Format("parameter_h/hPARs_Mu_%i_%1.1f.root", masspoint,BranchRatioForTgluon));
  hPARs[1]->SaveAs(TString::Format("parameter_h/hPARs_Ele_%i_%1.1f.root", masspoint,BranchRatioForTgluon));

  printf("[JACKY TEST] COV_MATRIX end\n");
  const TMatrixDSymEigen eigen_data(covMatrix);
  eigenValues = eigen_data.GetEigenValues();
  eigenValues.Sqrt();
  eigenValues.Print();
  eigenVectors = eigen_data.GetEigenVectors();
  eigenVectors.Print();

  fit_data->setParLimits(POIINDEX, 0.0, PAR_MAX[POIINDEX]); // for posterior calculation, signal has to be positive
  TGraph* post_data = 0;
  if(useMCMC==0)
  {
      post_data=fit_data->calculatePosterior(NSAMPLES);
      post_data->Write("post_0");
      cout << "Call limit reached: " << (fit_data->callLimitReached() ? "True" : "False") << endl;
  }
  else
  {
      post_data=fit_data->calculatePosterior(0);
      pair<double, double> statonly_bounds=evaluateInterval(post_data, ALPHA, LEFTSIDETAIL);
      fit_data->setParLimits(0, 0.0, xsUpperBoundFactor*(statonly_bounds.second));
      post_data=fit_data->calculatePosterior(NSAMPLES, useMCMC);
      //fit_data->PrintAllMarginalized("plots.ps");
      //fit_data->PrintResults("results.txt");
      post_data->Write("post_0");
  }

  // evaluate the limit
  pair<double, double> bounds_data=evaluateInterval(post_data, ALPHA, LEFTSIDETAIL);
  observedLowerBound=bounds_data.first;
  observedUpperBound=bounds_data.second;

  if(DEBUGMODE)
      for(int i=0; i<NPARS; i++) printf("[End for Obs : Parameter-%i] %s : %f, %f (%f -- %f) NUIS(%i)\n",
          i, PAR_NAMES[i], PAR_GUESSES[i], PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);

  // reset the covariance matrix
  for(int i = 0; i<NPARS; ++i) { for(int j = 0; j<NPARS; ++j) COV_MATRIX[i][j]=0.; }

  //Varied paras of pseudata
  TVectorD datafit_eigenValues = TVectorD(NBKGPARS);
  TMatrixD datafit_eigenVectors = TMatrixD(NBKGPARS,NBKGPARS);
  datafit_eigenValues = eigenValues;
  datafit_eigenVectors = eigenVectors;
  double *pseudo_paras;
  pseudo_paras = fit_data->getParameters();
  TRandom3 *tr3;
  tr3 = new TRandom3(31415+AdditionRNDSeed);
  double nsigma[NBKGPARS] = {0.};
  double varied_para[NBKGPARS] = {0.};

  // perform the PEs
  for(int pe=1; pe<=NPES; ++pe) {

      //Varied paras of pseudata
      pseudo_paras = fit_data->getParameters();
      bool first_vari = true;

      //for(int i=0; i<NPARS; i++) printf("para: [%i] %s : %g\n",i, PAR_NAMES[i], fit_data->getParameter(i));
      //for(int i = 0; i < NBKGPARS; i++)std::cout << "before para: ["<< i+6 <<"] "<<pseudo_paras[i+6] << std::endl;

      while(pseudo_paras[0] < 0 || pseudo_paras[3] < 0 || first_vari){
          // Bkg for Tgluon channel
          for(int v1=0; v1<NBKGPARS; ++v1){
              nsigma[v1] = tr3->Gaus(0, 1);
              for(int k1=0; k1<NBKGPARS; ++k1) varied_para[k1]=nsigma[v1]*datafit_eigenValues(v1)*datafit_eigenVectors[k1][v1];
              //		    std::cout << "varied_para: ["<< k1 <<"] "<<varied_para[k1] << std::endl;
              for(int i = 0; i < NBKGPARS; i++) pseudo_paras[i+6] += varied_para[i];
          }
          /// Bkg for Tgamma channel
          for(int ipar = NPARS - NTgammaBkgShapeUnc - 1; ipar < NPARS; ipar++)
              pseudo_paras[ipar] = tr3->Gaus(0, 1);    
              //pseudo_paras[ipar] = tr3->Uniform(-1.*Nsigma_bkg_shape, Nsigma_bkg_shape);    
          first_vari = false;
      }
      if(BranchRatioForTgluon == 0.0){
          pseudo_paras[6] = 0;
          pseudo_paras[9] = 0;
      }else if(BranchRatioForTgluon == 1.0){
          for(int ipar = NPARS - NTgammaBkgShapeUnc - 1; ipar < NPARS; ipar++)
              pseudo_paras[ipar] = 0;    
      }
      for(int i = 0; i < NBKGPARS; i++) std::cout << "after para: ["<< i+6 <<"] "<<pseudo_paras[i+6] << std::endl;

      cout << "*********** pe=" << pe << " ***********" << endl;
      ostringstream pestr;
      pestr << "_" << pe;

      // setup the fitter with the input from the signal+background fit
      fit_data->setParameter(POIINDEX, 0.0); // set the POI value to 0 to get the B component of the S+B fit (for calculating pulls and generating pseudo-data)
      //TH1D* hist = fit_data->makePseudoData((string("data")+pestr.str()).c_str()).at(0);
      //TH1D* hist2 = fit_data->makePseudoData((string("data2")+pestr.str()).c_str()).at(1);
      //TH1D* hist3 = fit_data->makePseudoData((string("data3")+pestr.str()).c_str()).at(2);

      TH1D* hist = fit_data->makePseudoData((string("data")+pestr.str()).c_str(), pseudo_paras).at(0);
      hist->Write((string("data")+pestr.str()).c_str());
      TH1D* hist2 = fit_data->makePseudoData((string("data2")+pestr.str()).c_str(), pseudo_paras).at(1);
      hist2->Write((string("data2")+pestr.str()).c_str());
      pseudo_paras[0] = 0;
      TH1D* hist3 = fit_data->makePseudoData((string("data3")+pestr.str()).c_str(), pseudo_paras).at(2);
      hist3->Write((string("data3")+pestr.str()).c_str());
      fit_data->setParameter(POIINDEX, POIval);


      fprintf(limitFile,"[DEBUG] %ith-Pseudo-Exp ( %f, %f ) in bin(5, 7)\n",pe, hist3->GetBinContent(5), hist3->GetBinContent(7));

      if(BranchRatioForTgluon == 0.0){
          // set Tgluon's data as zero, of course signal should be zero as above
          for(int i=1;i<=hist->GetXaxis()->GetNbins();i++)
              hist->SetBinContent(i,0.);
          for(int i=1;i<=hist2->GetXaxis()->GetNbins();i++)
              hist2->SetBinContent(i,0.);
      }else if(BranchRatioForTgluon == 1.0){
          for(int i=1;i<=hist3->GetXaxis()->GetNbins();i++)
              hist3->SetBinContent(i,0.);
      }

      printf("====================PE %i Fit for expected====================\n", pe);
      Fitter *fit = new Fitter(AdditionRNDSeed, hist, hist2, hist3, INTEGRAL);
      fit->setPOIIndex(POIINDEX);
      fit->setPrintLevel(0);
      for(int i=0; i<NPARS; i++) fit->defineParameter(i, PAR_NAMES[i], fit_data->getParameter(i), PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);

      // perform a signal+background fit possibly followed by a background-only fit with a fixed but non-zero signal
      for(int i=0; i<NPARS; i++) if(PAR_TYPE[i]>=2 || PAR_MIN[i]==PAR_MAX[i]) fit->fixParameter(i);
      if(BonlyFitForSyst) { fit->doFit(); fit->fixParameter(POIINDEX); }
      fit->doFit(&COV_MATRIX[0][0], NPARS);
      string fitStatus = fit->getFitStatus();
      printf("[debug] fitStatus : %s\n",fitStatus.c_str());
      if((fitStatus.find("CONVERGED")==string::npos&&(fitStatus.find("RESET")==string::npos))) continue; // skip the PE if the fit did not converge
      //if(fitStatus.find("CONVERGED")==string::npos) continue; // skip the PE if the fit did not converge
      fit->fixParameter(POIINDEX); // a parameter needs to be fixed before its value can be changed
      double exp_POIval = fit->getParameter(POIINDEX);
      fit->setParameter(POIINDEX, 0.0); // set the POI value to 0 to get the B component of the S+B fit (for calculating pulls and generating pseudo-data)
      fit->calcPull((string("pull_bkg_mu_")+pestr.str()).c_str(), _mu)->Write();
      fit->calcDiff((string("diff_bkg_mu_")+pestr.str()).c_str(), _mu)->Write();
      fit->calcPull((string("pull_bkg_ele_")+pestr.str()).c_str(), _ele)->Write();
      fit->calcDiff((string("diff_bkg_ele_")+pestr.str()).c_str(), _ele)->Write();
      fit->write((string("fit_bkg")+pestr.str()).c_str());

      //Fill int nuisance distribution
      for(int i = 0; i < 6; i++){
          nui_dis[i]->Fill(fit->getParameter(6+i));
      }
      nui_dis[6]->Fill(exp_POIval);

      // calculate eigenvalues and eigenvectors
      for(int i = 0; i<NBKGPARS; ++i) { for(int j = 0; j<NBKGPARS; ++j) { covMatrix(i,j)=COV_MATRIX[i+shift][j+shift]; } }
      covMatrix.Print();

      if( BranchRatioForTgluon == 0.0 ){
          // ignore fit parameters in a pure T+gamma channel
          for(int i = 0; i<NBKGPARS; ++i) { for(int j = 0; j<NBKGPARS; ++j) { covMatrix(i,j)=0.0; } }   
      }

      const TMatrixDSymEigen eigen(covMatrix);
      eigenValues = eigen.GetEigenValues();
      eigenValues.Print();
      bool hasNegativeElement = false;
      for(int i = 0; i<NBKGPARS; ++i) {if(eigenValues(i)<0.) hasNegativeElement = true; }
      if(hasNegativeElement) continue; // this is principle should never happen. However, if it does, skip the PE
      eigenValues.Sqrt();
      eigenVectors = eigen.GetEigenVectors();
      eigenVectors.Print();

      fit->setParLimits(POIINDEX, 0.0, PAR_MAX[POIINDEX]); // for posterior calculation, signal has to be positive
      TGraph* post = 0;
      if(useMCMC==0)
      {
          post=fit->calculatePosterior(NSAMPLES);
          post->Write((string("post")+pestr.str()).c_str());
          cout << "Call limit reached in pe=" << pe << ": " << (fit->callLimitReached() ? "True" : "False") << endl;
      }
      else
      {
          post=fit->calculatePosterior(0);
          pair<double, double> statonly_bounds=evaluateInterval(post, ALPHA, LEFTSIDETAIL);
          fit->setParLimits(0, 0.0, xsUpperBoundFactor*(statonly_bounds.second));
          post=fit->calculatePosterior(NSAMPLES, useMCMC);
          post->Write((string("post")+pestr.str()).c_str());
      }
      //TGraph* post=fit->calculatePosterior(NSAMPLES);
      //cout << "Call limit reached in pe=" << pe << ": " << (fit->callLimitReached() ? "True" : "False") << endl;

      // evaluate the limit
      pair<double, double> bounds=evaluateInterval(post, ALPHA, LEFTSIDETAIL);
      if(bounds.first==0. && bounds.second>0.)
      {
          expectedLowerBounds.push_back(bounds.first);
          expectedUpperBounds.push_back(bounds.second);
      }

      sprintf(buffer,"Upper limit on %f",bounds.second);
      post->SetTitle(buffer);
      post->Write((string("post")+pestr.str()).c_str());
      // reset the covariance matrix
      for(int i = 0; i<NPARS; ++i) { for(int j = 0; j<NPARS; ++j) COV_MATRIX[i][j]=0.; }
  }

  //Saving the upper bound distribution (from pseudodata)
  TH1D* expected_limit_upper = new TH1D("exp_limit_upper", "exp_limit_upper", 100, 0, 0.5);
  for(unsigned int i = 0; i < expectedUpperBounds.size(); i++){expected_limit_upper->Fill(expectedUpperBounds[i]);}
  ostringstream exp_limit_upper;
  exp_limit_upper<< "exp_limit_dis/exp_limit_upper_"<< masspoint<<"_"<<BranchRatioForTgluon<<"_"<<AdditionRNDSeed<<".root";
  expected_limit_upper->SaveAs(exp_limit_upper.str().c_str());

  TH1D* expected_limit_lower = new TH1D("exp_limit_lower", "exp_limit_lower", 100, 0, 0.5);
  for(unsigned int i = 0; i < expectedLowerBounds.size(); i++){expected_limit_lower->Fill(expectedUpperBounds[i]);}
  ostringstream exp_limit_lower;
  exp_limit_lower<< "exp_limit_dis/exp_limit_lower_"<< masspoint<<"_"<<BranchRatioForTgluon<<"_"<<AdditionRNDSeed<<".root";
  expected_limit_lower->SaveAs(exp_limit_lower.str().c_str());

  nui_dis[0]->SaveAs(TString::Format("nuisance/mu_p1_%i_nui.root", masspoint));
  nui_dis[1]->SaveAs(TString::Format("nuisance/mu_p2_%i_nui.root", masspoint));
  nui_dis[2]->SaveAs(TString::Format("nuisance/mu_p3_%i_nui.root", masspoint));
  nui_dis[3]->SaveAs(TString::Format("nuisance/ele_p1_%i_nui.root", masspoint));
  nui_dis[4]->SaveAs(TString::Format("nuisance/ele_p2_%i_nui.root", masspoint));
  nui_dis[5]->SaveAs(TString::Format("nuisance/ele_p3_%i_nui.root", masspoint));
  nui_dis[6]->SaveAs(TString::Format("nuisance/sig_%i_nui.root", masspoint));

  ////////////////////////////////////////////////////////////////////////////////
  // print the results
  ////////////////////////////////////////////////////////////////////////////////

  cout << "**********************************************************************" << endl;
  for(unsigned int i=0; i<expectedLowerBounds.size(); i++) {
      cout << "expected bound(" << (i+1) << ") = [ " << expectedLowerBounds[i] << " , " << expectedUpperBounds[i] << " ]" << endl;
      fprintf(limitFile, "expected bound = [ %1.6f , %1.6f ] \n", expectedLowerBounds[i], expectedUpperBounds[i]);
  }

  cout << "\nobserved bound = [ " << observedLowerBound << " , " << observedUpperBound << " ]" << endl;

  fprintf(limitFile, "observed bound = [ %1.6f , %1.6f ] \n", observedLowerBound, observedUpperBound);

  if(LEFTSIDETAIL>0.0 && NPES>0) {
      cout << "\n***** expected lower bounds *****" << endl;
      double median;
      pair<double, double> onesigma;
      pair<double, double> twosigma;
      getQuantiles(expectedLowerBounds, median, onesigma, twosigma);
      cout << "median: " << median << endl;
      cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
      cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;

      fprintf(limitFile, "\n***** expected lower bounds *****\n");
      fprintf(limitFile, "median: %1.6f \n", median);
      fprintf(limitFile, "+/-1 sigma band = [ %1.6f , %1.6f ] \n", onesigma.first, onesigma.second);
      fprintf(limitFile, "+/-2 sigma band = [ %1.6f , %1.6f ] \n", twosigma.first, twosigma.second);
  }

  if(LEFTSIDETAIL<1.0 && NPES>0) {
      cout << "\n***** expected upper bounds *****" << endl;
      double median;
      pair<double, double> onesigma;
      pair<double, double> twosigma;
      getQuantiles(expectedUpperBounds, median, onesigma, twosigma);
      cout << "median: " << median << endl;
      cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
      cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;

      fprintf(limitFile, "\n***** expected upper bounds *****\n");
      fprintf(limitFile, "median: %1.6f \n", median);
      fprintf(limitFile, "+/-1 sigma band = [ %1.6f , %1.6f ] \n", onesigma.first, onesigma.second);
      fprintf(limitFile, "+/-2 sigma band = [ %1.6f , %1.6f ] \n", twosigma.first, twosigma.second);
  }
  fclose(limitFile);

  // close the output file
  rootfile->Close();    // casue "Segmentation fault"

  return 0;
}

double GetCombinedSysInTGTGAnalysis(TFile* histfile, int ithBin){
    int masspoint = (int) SIGMASS;
    ostringstream syshistTGTG, syshistTgTg, syshistTGTg;
    syshistTGTG << "hSigSysNumber_TGTG_" << masspoint ;
    syshistTgTg << "hSigSysNumber_TgTg_" << masspoint ;
    syshistTGTg << "hSigSysNumber_TGTg_" << masspoint ;
    ostringstream histTGTG, histTgTg, histTGTg;
    histTGTG << "hSigMC_TGTG_" << masspoint ;
    histTgTg << "hSigMC_TgTg_" << masspoint ;
    histTGTg << "hSigMC_TGTg_" << masspoint ;

    TH1D *SigSysNumberTGTG = (TH1D*)histfile->Get(syshistTGTG.str().c_str());
    TH1D *SigSysNumberTGTg = (TH1D*)histfile->Get(syshistTGTg.str().c_str());
    TH1D *SigSysNumberTgTg = (TH1D*)histfile->Get(syshistTgTg.str().c_str());
    TH1D *hSigMC_TGTG = (TH1D*)histfile->Get(histTGTG.str().c_str());
    TH1D *hSigMC_TGTg = (TH1D*)histfile->Get(histTGTg.str().c_str());
    TH1D *hSigMC_TgTg = (TH1D*)histfile->Get(histTgTg.str().c_str());

    double Y_TGTG = hSigMC_TGTG->Integral() * BranchRatioForTgluon*BranchRatioForTgluon;
    double Err_TGTG = Y_TGTG*SigSysNumberTGTG->GetBinContent(ithBin);
    double Y_TGTg = hSigMC_TGTg->Integral() * 2.* BranchRatioForTgluon*(1.-BranchRatioForTgluon);
    double Err_TGTg = Y_TGTg*SigSysNumberTGTg->GetBinContent(ithBin);
    double Y_TgTg = hSigMC_TgTg->Integral() * (1.-BranchRatioForTgluon) * (1.-BranchRatioForTgluon);
    double Err_TgTg = Y_TgTg*SigSysNumberTgTg->GetBinContent(ithBin);

    delete SigSysNumberTGTG;
    delete SigSysNumberTGTg;
    delete SigSysNumberTgTg;
    delete hSigMC_TGTG;
    delete hSigMC_TGTg;
    delete hSigMC_TgTg;

    return sqrt(Err_TGTG*Err_TGTG+Err_TGTg*Err_TGTg+Err_TgTg*Err_TgTg)/(Y_TGTG+Y_TGTg+Y_TgTg);
}
