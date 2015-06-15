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

#include "binneddata.hh"
#include "fit.hh"
#include "statistics.hh"

////////////////////////////////////////////////////////////////////////////////
// magic numbers
////////////////////////////////////////////////////////////////////////////////

// number of pseudoexperiments
const int NPES=0;

// number of samples of nuisance parameters for Bayesian MC integration
const int NSAMPLES=2E5;

// alpha (1-alpha=confidence interval)
const double ALPHA=0.05;

// left side tail
const double LEFTSIDETAIL=0.0;

// output file name
const char* OUTPUTFILE="stats.root";

// input file name
const char* INPUTFILE="dijet_mass_HT_fat_1p010fbm1.txt";

// histogram binning
const int NBINS=34;
double BOUNDARIES[NBINS] = { 838,
			     890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383,
			     1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132,
			     2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147,
                             3275, 3403, 3531, 3659, 3787, 3915 };

// parameters
double SIGMASS=0;
const int NPARS=8;
const int POIINDEX=0; // which parameter is "of interest"
const char* PAR_NAMES[8]    = { "xs", "lumi", "jes", "jer", "bkg norm",  "p1",  "p2", "p3" };
const double PAR_GUESSES[8] = {  0.1,  1010.,   1.0,   1.0,     5.6E-2,   7.4,   6.3,  0.2 };
const double PAR_MIN[8]     = {  0.0,    0.0,   0.0,   0.0,         0.,   0.0,   0.0,  0.0 };
const double PAR_MAX[8]     = { 1.E6,  5000.,   2.0,   2.0,       10.0, 100.0, 100.0, 10.0 };
const double PAR_ERR[8]     = { 0.01,   40.4,  0.04,  0.10,      1.E-3,   1.0,   1.0,  0.1 };
const int PAR_TYPE[8]       = {    1,      1,     1,     1,          0,     0,     0,    0 }; // 1 = signal, 0 = background
const int PAR_NUIS[8]       = {    0,      0,     0,     0,          0,     0,     0,    0 }; // 1 = nuisance parameter, 0 = not varied (the POI is not a nuisance parameter)

TH1D* HIST=0; // signal histogram
TH1D* HISTCDF=0; // signal CDF

////////////////////////////////////////////////////////////////////////////////
// function definition (not really used; only the integral is used)
////////////////////////////////////////////////////////////////////////////////
double FCN(double *x, double *par)
{
  double invmass=x[0];
  double xs=par[0];
  double lumi=par[1];
  double jes=par[2];
  double jer=par[3];
  double norm=par[4];
  double p1=par[5];
  double p2=par[6];
  double p3=par[7];

  double bkg = norm*pow(1.0-invmass/7000.0,p1)/pow(invmass/7000.0,p2+p3*log(invmass/7000.0));
  double mass=jes*(jer*(invmass-SIGMASS)+SIGMASS);
  int bin=HIST->GetXaxis()->FindBin(mass);
  double sig = xs*lumi*HIST->GetBinContent(bin);
  return bkg+sig;
}

////////////////////////////////////////////////////////////////////////////////
// function integral
////////////////////////////////////////////////////////////////////////////////
double INTEGRAL(double *x0, double *xf, double *par)
{
  double xs=par[0];
  double lumi=par[1];
  double jes=par[2];
  double jer=par[3];
  double norm=par[4];
  double p1=par[5];
  double p2=par[6];
  double p3=par[7];

  // uses Simpson's 3/8th rule to compute the background integral over a short interval
  // also use a power series expansion to determine the intermediate intervals since the pow() call is expensive

  double dx=(xf[0]-x0[0])/3./7000.;
  double x=x0[0]/7000.0;
  double logx=log(x);

  double a=pow(1-x,p1)/pow(x,p2+p3*logx);
  double b=dx*a/x/(x-1)*(p2+p1*x-p2*x-2*p3*(x-1)*logx);
  double c=0.5*dx*dx*a*( (p1-1)*p1/(x-1)/(x-1) - 2*p1*(p2+2*p3*logx)/(x-1)/x + (p2+p2*p2-2*p3+2*p3*logx*(1+2*p2+2*p3*logx))/x/x );
  double d=0.166666667*dx*dx*dx*a*( (p1-2)*(p1-1)*p1/(x-1)/(x-1)/(x-1) - 3*(p1-1)*p1*(p2+2*p3*logx)/(x-1)/(x-1)/x - (1+p2+2*p3*logx)*(p2*(2+p2) - 6*p3 + 4*p3*logx*(1+p2*p3*logx))/x/x/x + 3*p1*(p2+p2*p2-2*p3+2*p3*logx*(1+2*p2+2*p3*logx))/(x-1)/x/x );

  double bkg=(xf[0]-x0[0])*norm*(a+0.375*(b+c+d)+0.375*(2*b+4*c+8*d)+0.125*(3*b+9*c+27*d));

  if(xs==0.0) return bkg;
     
  double xprimef=jes*(jer*(xf[0]-SIGMASS)+SIGMASS);
  double xprime0=jes*(jer*(x0[0]-SIGMASS)+SIGMASS);
  int bin1=HISTCDF->GetXaxis()->FindBin(xprimef);
  int bin2=HISTCDF->GetXaxis()->FindBin(xprime0);
  if(bin1<1) bin1=1;
  if(bin1>HISTCDF->GetNbinsX()) bin1=HISTCDF->GetNbinsX();
  if(bin2<1) bin1=1;
  if(bin2>HISTCDF->GetNbinsX()) bin2=HISTCDF->GetNbinsX();
  double sig=xs*lumi*(HISTCDF->GetBinContent(bin1)-HISTCDF->GetBinContent(bin2));

  return bkg+sig;
}

////////////////////////////////////////////////////////////////////////////////
// main function
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  if(argc<=1) {
    std::cout << "Usage: stats signalmass" << std::endl;
    return 0;
  }

  // setup the signal histogram
  TFile* histfile=new TFile("Test_Resonance_Shapes.root");
  histfile->cd();
  SIGMASS = std::atof(argv[1]);
  int masspoint = static_cast<int>(SIGMASS);
  std::ostringstream histname, cdfname;
  histname << "h_qstar_" << masspoint;
  cdfname << "h_qstar_" << masspoint << "_cdf";
  HIST=dynamic_cast<TH1D*>(gROOT->FindObject(histname.str().c_str()));
  HISTCDF=dynamic_cast<TH1D*>(gROOT->FindObject(cdfname.str().c_str()));
  assert(HIST && HISTCDF && SIGMASS>0);
  HIST->Scale(5.); // proper normalization

  // create the output file
  TFile* rootfile=new TFile(OUTPUTFILE, "RECREATE");  rootfile->cd();

  // get the data
  TH1D* data=getData(INPUTFILE, "data_0", NBINS-1, BOUNDARIES);

  // setup an initial fitter to perform a background-only fit
  Fitter initfit(data, INTEGRAL);
  for(int i=0; i<NPARS; i++) initfit.defineParameter(i, PAR_NAMES[i], PAR_GUESSES[i], PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);

  // do an initial background-only fit, first
  for(int i=0; i<NPARS; i++) if(PAR_TYPE[i]==1) initfit.fixParameter(i);
  initfit.setParameter(POIINDEX, 0.0); // set the POI value to 0
  initfit.doFit();
  initfit.calcPull("pull_bkg_init")->Write();
  initfit.calcDiff("diff_bkg_init")->Write();
  initfit.write("fit_bkg_init");

  // setup the limit values
  double observedLowerBound, observedUpperBound;
  std::vector<double> expectedLowerBounds;
  std::vector<double> expectedUpperBounds;

  // perform the PEs (0 = data)
  for(int pe=0; pe<=NPES; ++pe) {

    std::cout << "*********** pe=" << pe << " ***********" << std::endl;
    std::ostringstream pestr;
    pestr << "_" << pe;

    // setup the fitter with the input from the background-only fit
    TH1D* hist = (pe==0) ? data : initfit.makePseudoData((std::string("data")+pestr.str()).c_str());
    Fitter fit(hist, INTEGRAL);
    fit.setPOIIndex(POIINDEX);
    fit.setPrintLevel(0);
    for(int i=0; i<NPARS; i++) fit.defineParameter(i, PAR_NAMES[i], initfit.getParameter(i), PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);

    // perform a background-only fit
    for(int i=0; i<NPARS; i++) if(PAR_TYPE[i]==1) fit.fixParameter(i);
    fit.setParameter(POIINDEX, 0.0); // set the POI value to 0
    fit.doFit();
    fit.calcPull((std::string("pull_bkg")+pestr.str()).c_str())->Write();
    fit.calcDiff((std::string("diff_bkg")+pestr.str()).c_str())->Write();
    fit.write((std::string("fit_bkg")+pestr.str()).c_str());

    // fix the ranges for the background parameters before calculating the posterior
    for(int i=0; i<NPARS; i++) {
      if(PAR_TYPE[i]==0 && PAR_NUIS[i]==1) {
	double val, err;
	fit.getParameter(i, val, err);
	fit.setParLimits(i, val-err, val+err);
      }
    }

    observedLowerBound=0.;
    observedUpperBound=fit.calculateUpperBoundWithCLs(NSAMPLES, ALPHA);
    

    // put the ranges back in place
    for(int i=0; i<NPARS; i++) {
      if(PAR_TYPE[i]==0 && PAR_NUIS[i]==1) {
	fit.setParLimits(i, PAR_MIN[i], PAR_MAX[i]);
      }
    }

    // evaluate the limit
    /*    std::pair<double, double> bounds=evaluateInterval(post, ALPHA, LEFTSIDETAIL);
    if(pe==0) {
      observedLowerBound=bounds.first;
      observedUpperBound=bounds.second;
    } else {
      expectedLowerBounds.push_back(bounds.first);
      expectedUpperBounds.push_back(bounds.second);
      }*/
  }

  ////////////////////////////////////////////////////////////////////////////////
  // print the results
  ////////////////////////////////////////////////////////////////////////////////

  std::cout << "**********************************************************************" << std::endl;
  for(unsigned int i=0; i<expectedLowerBounds.size(); i++)
    std::cout << "expected bound(" << (i+1) << ") = [ " << expectedLowerBounds[i] << " , " << expectedUpperBounds[i] << " ]" << std::endl;

  std::cout << "\nobserved bound = [ " << observedLowerBound << " , " << observedUpperBound << " ]" << std::endl;

  if(LEFTSIDETAIL>0.0 && NPES>0) {
    std::cout << "\n***** expected lower bounds *****" << std::endl;
    double median;
    std::pair<double, double> onesigma;
    std::pair<double, double> twosigma;
    getQuantiles(expectedLowerBounds, median, onesigma, twosigma);
    std::cout << "median: " << median << std::endl;
    std::cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << std::endl;
    std::cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << std::endl;
  }

  if(LEFTSIDETAIL<1.0 && NPES>0) {
    std::cout << "\n***** expected upper bounds *****" << std::endl;
    double median;
    std::pair<double, double> onesigma;
    std::pair<double, double> twosigma;
    getQuantiles(expectedUpperBounds, median, onesigma, twosigma);
    std::cout << "median: " << median << std::endl;
    std::cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << std::endl;
    std::cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << std::endl;
  }

  // close the output file
  rootfile->Close();

  return 0;
}
