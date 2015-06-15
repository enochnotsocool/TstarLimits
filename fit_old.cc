#include <cassert>
#include <sstream>
#include <cmath>
#include <limits>
#include <vector>

#include <TF1.h>
#include <TMath.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TGraph.h>

#include "fit.hh"
#include "binneddata.hh"
#include "statistics.hh"

Fitter* Fitter::theFitter_;
TRandom3* Fitter::rand_;

Fitter::Fitter(int AdditionRNDSeed)
{
  data_=0;
  data2_=0;
  data3_=0;
  functionIntegral_=0;
  printlevel_=3;
  //strategy_=2;
  strategy_=0;
  minuit_.SetErrorDef(0.5); // likelihood
  if(!rand_) rand_=new TRandom3(31415+AdditionRNDSeed);
  parameters_=0;
  poiIndex_=-1;
  callLimitReached_=false;
  poiBestFit_ = 0;
  poiUserError_ = 0;
}

Fitter::Fitter(int AdditionRNDSeed, TH1D* data, TH1D* data2, TH1D* data3, integral_ptr_t functionIntegral)
{
  data_=data;
  data2_=data2;
  data3_=data3;
  functionIntegral_=functionIntegral;
  printlevel_=3;
  //strategy_=2;
  strategy_=0;
  minuit_.SetErrorDef(0.5); // likelihood
  if(!rand_) rand_=new TRandom3(31415+AdditionRNDSeed);
  parameters_=0;
  poiIndex_=-1;
  callLimitReached_=false;
  poiBestFit_ = 0;
  poiUserError_ = 0;
}

Fitter::~Fitter()
{
  if(parameters_) delete[] parameters_;
}

void Fitter::doFit(void)
{
  // setup the fitter so the information can be retrieved by nll
  Fitter::theFitter_=this;

  // setup TMinuit
  minuit_.SetPrintLevel(printlevel_);

  // set the strategy
  std::ostringstream command;
  command << "SET STR " << strategy_;
  minuit_.Command(command.str().c_str());

  // do the fit
  minuit_.SetFCN(nll);
  Double_t arglist[1] = {20000.0};
  Int_t err = 0;
  minuit_.mnexcm("MIGRAD",arglist,1,err);

  return;
}

void Fitter::doFit(double* emat, int ndim)
{
  // setup the fitter so the information can be retrieved by nll
  Fitter::theFitter_=this;

  // setup TMinuit
  minuit_.SetPrintLevel(printlevel_);

  // set the strategy
  std::ostringstream command;
  command << "SET STR " << strategy_;
  minuit_.Command(command.str().c_str());

  // do the fit
  minuit_.SetFCN(nll);
  Double_t arglist[1] = {20000.0};
  Int_t err = 0;
  minuit_.mnexcm("MIGRAD",arglist,1,err);
  minuit_.mnemat(emat, ndim);

  return;
}

TH1D* Fitter::calcPull(const char* name, int which)
{
  // debug -- start --
  if(DEBUGMODE){
      for(int i=0; i<minuit_.GetNumPars(); i++) {
          double value=0;
          double err=0;
          getParameter(i, value, err);
          printf("[in calcPull %i] p(%i) : %f\n", minuit_.GetNumPars(), i, value);
      }
  }
  // debug -- end --

  const double alpha = 1 - 0.6827;
  TH1D* hPull;
  if(which == 1)  hPull=dynamic_cast<TH1D*>(data_->Clone(name));
  if(which == 2)  hPull=dynamic_cast<TH1D*>(data2_->Clone(name));
  if(which == 3)  hPull=dynamic_cast<TH1D*>(data3_->Clone(name));
  hPull->SetTitle("Pull Distribution");

  hPull->SetBinContent(0, 0.);
  hPull->SetBinContent(hPull->GetNbinsX()+1, 0.);
  hPull->SetBinError(0, 0.);
  hPull->SetBinError(hPull->GetNbinsX()+1, 0.);
  for(int bin=1; bin<=hPull->GetNbinsX(); bin++) {
    double binwidth=hPull->GetBinWidth(bin);
    double N = 0;
    if(which == 1) N=data_->GetBinContent(bin)*binwidth;
    if(which == 2) N=data2_->GetBinContent(bin)*binwidth;
    if(which == 3) N=data3_->GetBinContent(bin)*binwidth;
    double l = 0.5*TMath::ChisquareQuantile(alpha/2,2*N);
    double h = 0.5*TMath::ChisquareQuantile(1-alpha/2,2*(N+1));
    double el = N-l;
    double eh = h-N;

    double x0=hPull->GetBinLowEdge(bin);
    double xf=hPull->GetBinLowEdge(bin+1);
    double mu = 0;
    if(which == 1) mu=functionIntegral_(&x0, &xf, getParameters()).at(0);
    if(which == 2) mu=functionIntegral_(&x0, &xf, getParameters()).at(1);
    if(which == 3) mu=functionIntegral_(&x0, &xf, getParameters()).at(2);

    double err = (el + eh)/2.;
    if(N>=mu) err = el;
    if(N<mu) err = eh;
    
    double pull=(N-mu)/err;
    hPull->SetBinContent(bin, pull);
    hPull->SetBinError(bin, 1.0);
  }
  return hPull;
}

TH1D* Fitter::calcDiff(const char* name, int which)
{
  TH1D* hDiff;
  if(which == 1) hDiff=dynamic_cast<TH1D*>(data_->Clone(name));
  if(which == 2) hDiff=dynamic_cast<TH1D*>(data2_->Clone(name));
  if(which == 3) hDiff=dynamic_cast<TH1D*>(data3_->Clone(name));
  hDiff->SetTitle("Difference Distribution");

  hDiff->SetBinContent(0, 0.);
  hDiff->SetBinContent(hDiff->GetNbinsX()+1, 0.);
  hDiff->SetBinError(0, 0.);
  hDiff->SetBinError(hDiff->GetNbinsX()+1, 0.);
  for(int bin=1; bin<=hDiff->GetNbinsX(); bin++) {
    double binwidth=hDiff->GetBinWidth(bin);
    double N = 0;
    if(which == 1)N=data_->GetBinContent(bin)*binwidth;
    if(which == 2)N=data2_->GetBinContent(bin)*binwidth;
    if(which == 3)N=data3_->GetBinContent(bin)*binwidth;
    double err=Fitter::histError(N);

    double x0=hDiff->GetBinLowEdge(bin);
    double xf=hDiff->GetBinLowEdge(bin+1);
    double mu=0;
    if(which == 1)mu=functionIntegral_(&x0, &xf, getParameters()).at(0);
    if(which == 2)mu=functionIntegral_(&x0, &xf, getParameters()).at(1);
    if(which == 3)mu=functionIntegral_(&x0, &xf, getParameters()).at(2);

    double diff=(N-mu)/mu;
    hDiff->SetBinContent(bin, diff);
    hDiff->SetBinError(bin, err/mu);
  }
  return hDiff;
}

int Fitter::setParameter(int parno, double value)
{
  int err;
  double tmp[2];
  tmp[0]=parno+1;
  tmp[1]=value;

  minuit_.mnexcm("SET PAR", tmp, 2, err);
  return err;
}

int Fitter::setParLimits(int parno, double loLimit, double hiLimit)
{
  int err;
  double tmp[3];
  tmp[0]=parno+1;
  tmp[1]=loLimit;
  tmp[2]=hiLimit;

  minuit_.mnexcm("SET LIM", tmp, 3, err);
  return err;
}

double Fitter::getParameter(int parno) const
{
  double val, err;
  getParameter(parno, val, err);
  return val;
}

double Fitter::getParError(int parno) const
{
  double val, err;
  getParameter(parno, val, err);
  return err;
}

void Fitter::getParLimits(int parno, double &loLimit, double &hiLimit) const
{
  TString _name;
  double _val, _err;
  int _iuint;
  minuit_.mnpout(parno, _name, _val, _err, loLimit, hiLimit, _iuint);
  return;
}

double* Fitter::getParameters(void)
{
  // remove what is there
  if(parameters_) delete[] parameters_;

  int nPars=minuit_.GetNumPars();
  parameters_ = new double[nPars];
  for(int i=0; i<nPars; i++) {
    double value, err;
    getParameter(i, value, err);
    parameters_[i]=value;
  }

  return parameters_;
}

int Fitter::defineParameter(int parno, const char *name, double value, double error, double lo, double hi, int isNuisance)
{
  parameterIsNuisance_[parno]=isNuisance;
  if(poiIndex_>=0 && parno==poiIndex_)
    if(poiUserError_==0.) poiUserError_=error;
  return minuit_.DefineParameter(parno, name, value, error, lo, hi);
}

TGraph* Fitter::calculatePosterior(int nSamples)
{
  // we need a parameter of index defined
  assert(poiIndex_>=0);

  // set the number of samples for subsequent functions
  nSamples_=nSamples;
  nCalls_=0;
  nCallsMutation_=1;    // avoid the local minimum

  // fit for the best value of the POI
  int nPars=getNumPars();
  double loLimit, hiLimit;
  getParLimits(poiIndex_, loLimit, hiLimit);
  for(int i=0; i<nPars; i++) if(i==poiIndex_) floatParameter(i); else fixParameter(i);
  setParameter(poiIndex_, 0.1*(loLimit+hiLimit));
  doFit();
  fixParameter(poiIndex_);
  poiBestFit_=getParameter(poiIndex_);

  // setup Fitter
  Fitter::theFitter_=this;

  // now float all parameters
  for(int i=0; i<nPars; i++) if(i!=poiIndex_) floatParameter(i);

  // evalulate NLL at the POI fit value for future normalizations
  double nllNormalization=evalNLL();

  // recursively evaluate the posterior
  std::map<double, double> fcnEvaluations;
  evaluateForPosterior(loLimit, poiBestFit_, hiLimit, nllNormalization, fcnEvaluations);

  // dump the info into a graph
  int cntr=0;
  double maximumVal=-9999.;
  TGraph* graph=new TGraph(fcnEvaluations.size());
  for(std::map<double, double>::const_iterator it=fcnEvaluations.begin(); it!=fcnEvaluations.end(); ++it) {
    graph->SetPoint(cntr++, it->first, it->second);
    if(it->second>maximumVal) maximumVal=it->second;
  }

  // identify trivially small points on the left
  std::vector<int> pointsToRemove;
//   for(int i=0; i<graph->GetN()-1; i++) {
//     double x, y, nextx, nexty;
//     graph->GetPoint(i, x, y);
//     graph->GetPoint(i+1, nextx, nexty);
// 
//     if(y/maximumVal<1.E-3 && nexty/maximumVal<1.E-3) pointsToRemove.push_back(i);
//     else break;
//   }

  // identify trivially small points on the right
  for(int i=graph->GetN()-1; i>=1; i--) {
    double x, y, nextx, nexty;
    graph->GetPoint(i, x, y);
    graph->GetPoint(i-1, nextx, nexty);

    if(y/maximumVal<1.E-3 && nexty/maximumVal<1.E-3) pointsToRemove.push_back(i);
    else break;
  }

  // sort the points to remove from first to last
  std::sort(pointsToRemove.begin(), pointsToRemove.end());

  // remove the points
  for(int i=pointsToRemove.size()-1; i>=0; i--)
    graph->RemovePoint(pointsToRemove[i]);

  return graph;
}

double Fitter::evalNLL(void)
{
  Fitter::theFitter_=this;
  int a;
  double f;
  nll(a, 0, f, getParameters(), 0);
  return f;
}

double Fitter::histError(double val)
{
  const double alpha = 1 - 0.6827;

  if(val<25. && val>0.) return (0.5*TMath::ChisquareQuantile(1-alpha/2,2*(val+1))-0.5*TMath::ChisquareQuantile(alpha/2,2*val))/2.0; // this is not exactly correct since it symmetrizes what are otherwise asymmetric error bars
  else if(val==0.) return 0.5*TMath::ChisquareQuantile(1-alpha/2,2*(val+1)); // special case of 0 events with one-sided error bar
  else return sqrt(val); // for val>25 error bars with correct coverage are essentially symmetric
}

void Fitter::nll(int &, double *, double &f, double *par, int)
{
  assert(Fitter::theFitter_);
  TH1D* data=Fitter::theFitter_->data_;
  TH1D* data2=Fitter::theFitter_->data2_;
  TH1D* data3=Fitter::theFitter_->data3_;
  integral_ptr_t functionIntegral=Fitter::theFitter_->functionIntegral_;
  if(data->GetNbinsX() != data2->GetNbinsX() ) 
  {printf("data1 and data2 has different #bin, please check"); return;}

  f=0.0;
  // Tgluon Channel
  for(int bin=1; bin<=data->GetNbinsX(); bin++) {
    double binwidth=data->GetBinWidth(bin);
    double N=data->GetBinContent(bin)*binwidth;
    double N2=data2->GetBinContent(bin)*binwidth;

    double x0=data->GetBinLowEdge(bin);
    double xf=data->GetBinLowEdge(bin+1);
    double mu=functionIntegral(&x0, &xf, par).at(0);
    double mu2=functionIntegral(&x0, &xf, par).at(1);

    if(N==0.0) f += mu;
    else f -= (N*TMath::Log(mu) - TMath::LnGamma(N+1) - mu);
    if(N2==0.0) f += mu2;
    else f -= (N2*TMath::Log(mu2) - TMath::LnGamma(N2+1) - mu2);
	//printf("f if : %f\n", f);
	if(fabs(f) > 1000) return;
  }
  // Tgamma Channel
  for(int bin=1; bin<=data3->GetNbinsX(); bin++) {
    double binwidth=data3->GetBinWidth(bin);
    double N3=data3->GetBinContent(bin)*binwidth;

    double x0=data3->GetBinLowEdge(bin);
    double xf=data3->GetBinLowEdge(bin+1);
    double mu3=functionIntegral(&x0, &xf, par).at(2);

    if(N3==0.0) f += mu3;
    else f -= (N3*TMath::Log(mu3) - TMath::LnGamma(N3+1) - mu3);
    if(DEBUGMODE)
        printf("f if : %f , N3 : %f , mu3 : %f (%i), p19-20(%f, %f)\n", f, N3, mu3,bin, par[18], par[19]);

    //if(bin == 7)
    //    printf("f if : %f , N3 : %f , mu3 : %f (%i), p0,19-20(%f, %f, %f)\n", f, N3, mu3,bin, par[0], par[18], par[19]);
	if(fabs(f) > 1000) return;
  }
  return;
}

//TH1D* Fitter::makePseudoData(const char* name, double* parameters)
//std::pair <TH1D*,TH1D*> Fitter::makePseudoData(const char* name, double* parameters)
std::vector <TH1D*> Fitter::makePseudoData(const char* name, double* parameters)
{
  if(!parameters) parameters=getParameters();

  // start with a copy of the original dataset
  TH1D* hData=dynamic_cast<TH1D*>(data_->Clone(name));
  TH1D* hData2=dynamic_cast<TH1D*>(data2_->Clone(name));
  TH1D* hData3=dynamic_cast<TH1D*>(data3_->Clone(name));

  for(int bin=1; bin<=hData->GetNbinsX(); ++bin) {
    double lobin=hData->GetBinLowEdge(bin);
    double hibin=hData->GetBinLowEdge(bin+1);
    double integral=functionIntegral_(&lobin, &hibin, parameters).at(0);
    double val=rand_->Poisson(integral);
    double binWidth=hData->GetBinWidth(bin);
    hData->SetBinContent(bin, val/binWidth);
    hData->SetBinError(bin, histError(val)/binWidth);
  }
  for(int bin=1; bin<=hData2->GetNbinsX(); ++bin) {
    double lobin=hData2->GetBinLowEdge(bin);
    double hibin=hData2->GetBinLowEdge(bin+1);
    double integral=functionIntegral_(&lobin, &hibin, parameters).at(1);
    double val=rand_->Poisson(integral);
    double binWidth=hData2->GetBinWidth(bin);
    hData2->SetBinContent(bin, val/binWidth);
    hData2->SetBinError(bin, histError(val)/binWidth);
  }
  for(int bin=1; bin<=hData3->GetNbinsX(); ++bin) {
    double lobin=hData3->GetBinLowEdge(bin);
    double hibin=hData3->GetBinLowEdge(bin+1);
    double integral=functionIntegral_(&lobin, &hibin, parameters).at(2);
    double val=rand_->Poisson(integral);
    double binWidth=hData3->GetBinWidth(bin);
    hData3->SetBinContent(bin, val/binWidth);
    hData3->SetBinError(bin, histError(val)/binWidth);
  }
  std::vector <TH1D*> _pseudata;
  _pseudata.push_back(hData);
  _pseudata.push_back(hData2);
  _pseudata.push_back(hData3);
  //return hData;
  return _pseudata;
}


void Fitter::evaluateForPosterior(double lo, double mid, double hi, double nllNormalization, std::map<double, double>& fcnEval_)
{
  if((nCalls_++)>1000) {
    callLimitReached_=true;
    return;
  }

  // get the low value
  std::map<double, double>::iterator findit;
  findit = fcnEval_.find(lo);
  double loVal;
  if(findit==fcnEval_.end()) {
    loVal=computeLikelihoodWithSystematics(lo, nllNormalization);
    fcnEval_[lo]=loVal;
  } else {
    loVal=fcnEval_[lo];
  }

  // get the middle value
  findit = fcnEval_.find(mid);
  double midVal;
  if(findit==fcnEval_.end()) {
    midVal=computeLikelihoodWithSystematics(mid, nllNormalization);
    fcnEval_[mid]=midVal;
  } else {
    midVal=fcnEval_[mid];
  }

  // get the high value
  findit = fcnEval_.find(hi);
  double hiVal;
  if(findit==fcnEval_.end()) {
    hiVal=computeLikelihoodWithSystematics(hi, nllNormalization);
    fcnEval_[hi]=hiVal;
  } else {
    hiVal=fcnEval_[hi];
  }

  double maximumValX = 0.;
  double maximumVal = -999.;
  for(std::map<double, double>::const_iterator it=fcnEval_.begin(); it!=fcnEval_.end(); ++it)
    if(maximumVal<it->second)
    {
      maximumValX=it->first;
      maximumVal=it->second;
    }

  // for debugging
  //std::cout << "nCalls: " << nCalls_<<" ( nCallsMutation : "<<nCallsMutation_<<" )" << std::endl
  //          << "lo, mid, high: " << lo << ", " << mid << ", " << hi << std::endl
  //          << "loval, midval, hival: " << loVal << ", " << midVal << ", " << hiVal << std::endl
  //          << "maximumValX, maximumVal: " << maximumValX << ", " << maximumVal << std::endl << std::endl;

  if(fabs(loVal-midVal)>0.05*maximumVal || fabs(hiVal-midVal)>0.05*maximumVal) {
      if(fabs(hi-mid)/hi>0.01 && fabs(hi-mid)/poiBestFit_>0.01 && fabs(hi-mid)>poiUserError_) evaluateForPosterior(mid, 0.5*(mid+hi), hi, nllNormalization, fcnEval_); // imortant to go to the mid-high range first to get a nice falling posteriror tail in case the number of calls limit is reached
      if(fabs(lo-mid)/mid>0.01 && fabs(lo-mid)/poiBestFit_>0.01 && fabs(lo-mid)>poiUserError_) evaluateForPosterior(lo, 0.5*(lo+mid), mid, nllNormalization, fcnEval_);

      if(fabs(hi-mid)/hi>0.01 && fabs(hi-mid)/poiBestFit_>0.01 && fabs(hi-mid)>poiUserError_) evaluateForPosterior(mid, (2*mid+hi)/3., hi, nllNormalization, fcnEval_); // imortant to go to the mid-high range first to get a nice falling posteriror tail in case the number of calls limit is reached
      if(fabs(hi-mid)/hi>0.01 && fabs(hi-mid)/poiBestFit_>0.01 && fabs(hi-mid)>poiUserError_) evaluateForPosterior(mid, (mid+2*hi)/3., hi, nllNormalization, fcnEval_); // imortant to go to the mid-high range first to get a nice falling posteriror tail in case the number of calls limit is reached


      if(fabs(lo-mid)/mid>0.01 && fabs(lo-mid)/poiBestFit_>0.01 && fabs(lo-mid)>poiUserError_) evaluateForPosterior(lo, (2*lo+mid)/3., mid, nllNormalization, fcnEval_);
      if(fabs(lo-mid)/mid>0.01 && fabs(lo-mid)/poiBestFit_>0.01 && fabs(lo-mid)>poiUserError_) evaluateForPosterior(lo, (lo+2*mid)/3., mid, nllNormalization, fcnEval_);

  }

  if(RUNMUTATION&&rand_->Gaus(0,1) > 0.0 && (double)nCallsMutation_/(double)nCalls_ <0.1 ){   // 10% mutation rate
      nCallsMutation_++;
      double loLimitRND = 0;  
      double hiLimitRND = 0;
      getParLimits(poiIndex_, loLimitRND, hiLimitRND);
      evaluateForPosterior(loLimitRND, rand_->Uniform(loLimitRND,hiLimitRND), hiLimitRND, nllNormalization, fcnEval_);
  }

  return;
}

double Fitter::computeLikelihoodWithSystematics(double poiVal, double nllNormalization)
{
  // setup parameters
  int a;
  double f;
  double *pars=getParameters();
  pars[poiIndex_]=poiVal;

  // no systematics
  if(nSamples_==0) {
    nll(a, 0, f, pars, 0);
	//printf("liklihood : %f, %f, %f\n", poiVal, -f+nllNormalization, TMath::Exp(-f+nllNormalization));
    return TMath::Exp(-f+nllNormalization);
  }

  // setup nuisance parameter priors
  std::map<int, RandomPrior*> priors;
  for(std::map<int, int>::const_iterator it=parameterIsNuisance_.begin(); it!=parameterIsNuisance_.end(); ++it) {
    if(it->second>0) {
      assert(it->first!=poiIndex_);
      double parval, parerr;
      double lolim, uplim;
      getParameter(it->first, parval, parerr);
      getParLimits(it->first, lolim, uplim);
      if(DEBUGMODE)
          printf("Nuisance prior(%i) lolim,uplim,parval,parerr : %f, %f, %f, %f\n", it->first, lolim, uplim, parval, parerr);
      if(parameterIsNuisance_[it->first] != 4){
        if(lolim<parval-5*parerr) lolim=parval-5*parerr;
        if(uplim>parval+5*parerr) uplim=parval+5*parerr;
        //printf("parapara af%f, %f, %f, %f\n", lolim, parval, uplim, parerr);
	  }
      priors[it->first]=new RandomPrior(it->second, parval, parerr, lolim, uplim);
    }
  }

  // calculate average likelihood value over nuisance parameters
  double total=0.0;
  for(int sample=0; sample<=nSamples_; sample++) {
    for(std::map<int, RandomPrior*>::const_iterator it=priors.begin(); it!=priors.end(); ++it)
      pars[it->first]=it->second->getRandom();
    nll(a,0,f,pars,0);
    double like=TMath::Exp(-f+nllNormalization);
    if(isnan(f)) like = 0;
    //printf("Per. LH : %f, %f, %f\n", f,nllNormalization, like);

    //if(like>10) { // for debugging
    //  int nPars=minuit_.GetNumPars();
    //  std::cout << "sample=" << sample << std::endl;
    //  for(int i=0; i<nPars; i++) {
    //    std::cout << "pars[" << i << "]=" << pars[i] << std::endl;
    //  }
    //  std::cout << "like=" << like << "; f=" << f << "; norm=" << nllNormalization << std::endl;
    //  assert(0);
    //}

    total+=like;
  }

  // remove nuisance parameter priors
  for(std::map<int, RandomPrior*>::const_iterator it=priors.begin(); it!=priors.end(); ++it) delete it->second;
  //printf("Avg. LH : %f, %f, %i alpha_s(%f), normalNLL(%f)\n", total/nSamples_, total, nSamples_, pars[0], nllNormalization);
  return total/nSamples_;
}

double Fitter::calculateUpperBoundWithCLs(int nSamples, double alpha)
{
  // we need a parameter of index defined
  assert(poiIndex_>=0);

  // set the number of samples for subsequent functions
  nSamples_=nSamples;

  // fit for the best value of the POI
  int nPars=getNumPars();
  double loLimit, hiLimit;
  getParLimits(poiIndex_, loLimit, hiLimit);
  for(int i=0; i<nPars; i++) if(i==poiIndex_) floatParameter(i); else fixParameter(i);
  setParameter(poiIndex_, 0.1*(loLimit+hiLimit));
  doFit();
  double poiBestFit=getParameter(poiIndex_);
  double poiBestFitErr=getParError(poiIndex_);
  fixParameter(poiIndex_);

  // now float all parameters
  for(int i=0; i<nPars; i++) if(i!=poiIndex_) floatParameter(i);

  // scan the CLs in steps of error, first
  double poiVal;
  double prevPoiVal=0.0;
  double prevCLsVal=nSamples;
  double CLsVal=0.0;
  for(poiVal=poiBestFit; poiVal<hiLimit; poiVal+=poiBestFitErr) {

    // calculate CLs
    std::vector<double> CLb, CLsb;
    std::pair<int, int> numden=calculateCLs_(poiVal, CLb, CLsb);
    double A=static_cast<double>(numden.first)/nSamples;
    double B=static_cast<double>(numden.second)/nSamples;
    double Aerr=sqrt(A*(1-A)/nSamples);
    double Berr=sqrt(B*(1-B)/nSamples);
    CLsVal = B==0 ? nSamples*2 : A/B;
    double CLsErr = Aerr==0 || Berr==0 ? 0. : CLsVal*sqrt(Aerr*Aerr/A/A+Berr*Berr/B/B);
    double diff=fabs(CLsVal-alpha);

    std::cout << "Scan: CLsVal=" << CLsVal << "; poiVal=" << poiVal << std::endl;

    if(diff<CLsErr) return poiVal;

    if((prevCLsVal>=alpha && CLsVal<=alpha) || (prevCLsVal<=alpha && CLsVal>=alpha)) break;

    prevPoiVal=poiVal;
    prevCLsVal=CLsVal;
  }
  if(poiVal>=hiLimit) return -999.;

  // now, try to converge on best CLs point
  double poiValLo=prevPoiVal;
  double poiValHi=poiVal;
  double CLsValLo=prevCLsVal;
  double CLsValHi=CLsVal;
  int cntr=0;
  do {
    poiVal=(alpha-CLsValLo)*(poiValHi-poiValLo)/(CLsValHi-CLsValLo)+poiValLo;
    std::vector<double> CLb, CLsb;
    std::pair<int, int> numden=calculateCLs_(poiVal, CLb, CLsb);
    double A=static_cast<double>(numden.first)/nSamples;
    double B=static_cast<double>(numden.second)/nSamples;
    double Aerr=sqrt(A*(1-A)/nSamples);
    double Berr=sqrt(B*(1-B)/nSamples);
    double CLsVal = B==0 ? nSamples*2 : A/B;
    double CLsErr = Aerr==0 || Berr==0 ? 0. : CLsVal*sqrt(Aerr*Aerr/A/A+Berr*Berr/B/B);
    double diff=fabs(CLsVal-alpha);

    std::cout << "poiVal=" << poiVal << "; poiValLo=" << poiValLo << "; poiValHi=" << poiValHi
	      << "; CLsVal=" << CLsVal << "; CLsErr=" << CLsErr << "; CLSValLo=" << CLsValLo << "; CLsValHi=" << CLsValHi << std::endl;
    std::cout << "A=" << A << "; Aerr=" << Aerr
	      <<"; B=" << B << "; Berr=" << Berr << std::endl;

    if(diff>CLsErr) {
      if(alpha>CLsVal) {
	poiValHi=poiVal;
	CLsValHi=CLsVal;
      } else {
	poiValLo=poiVal;
	CLsValLo=CLsVal;
      }
    } else {
      break;
    }

  }while(cntr<100);

  return poiVal;
}

std::pair<int, int> Fitter::calculateCLs_(double poiVal, std::vector<double>& CLb, std::vector<double>& CLsb)
{
  // setup parameters
  int a;
  double f;
  double *pars=getParameters();

  // setup nuisance parameter priors
  std::map<int, RandomPrior*> priors;
  for(std::map<int, int>::const_iterator it=parameterIsNuisance_.begin(); it!=parameterIsNuisance_.end(); ++it) {
    if(it->second>0) {
      assert(it->first!=poiIndex_);
      double parval, parerr;
      double lolim, uplim;
      getParameter(it->first, parval, parerr);
      getParLimits(it->first, lolim, uplim);
      if(lolim<parval-5*parerr) lolim=parval-5*parerr;
      if(uplim>parval+5*parerr) uplim=parval+5*parerr;
      priors[it->first]=new RandomPrior(it->second, parval, parerr, lolim, uplim);
    }
  }

  TH1D* theData=data_;
  TH1D* theData2=data2_;
  TH1D* theData3=data3_;

  for(int i=0; i<nSamples_; i++) {

    // create the pseudodata
    TH1D* CLSB_pdata;
    TH1D* CLB_pdata;
    TH1D* CLSB_pdata2;
    TH1D* CLB_pdata2;
    TH1D* CLSB_pdata3;
    TH1D* CLB_pdata3;
    if(i>0) {
      for(std::map<int, RandomPrior*>::const_iterator it=priors.begin(); it!=priors.end(); ++it)
	pars[it->first]=it->second->getRandom();

      pars[poiIndex_]=poiVal;
      data_=theData;
      data2_=theData2;
      data3_=theData3;
      CLSB_pdata=makePseudoData("CLSB_pdata", pars).at(0);
      CLSB_pdata2=makePseudoData("CLSB_pdata", pars).at(1);
      CLSB_pdata3=makePseudoData("CLSB_pdata", pars).at(2);
      pars[poiIndex_]=0.0;
      CLB_pdata=makePseudoData("CLB_pdata", pars).at(0);
      CLB_pdata2=makePseudoData("CLB_pdata", pars).at(1);
      CLB_pdata3=makePseudoData("CLB_pdata", pars).at(2);
    } else {
      CLSB_pdata=theData;
      CLB_pdata=theData;
      CLSB_pdata2=theData2;
      CLB_pdata2=theData2;
      CLSB_pdata3=theData3;
      CLB_pdata3=theData3;
    }

    // CL_{s+b}
    data_=CLSB_pdata;
    data2_=CLSB_pdata2;
    data3_=CLSB_pdata3;
    pars[poiIndex_]=poiVal;
    nll(a, 0, f, pars, 0);
    double llnum=f;
    pars[poiIndex_]=0.0;
    nll(a, 0, f, pars, 0);
    double llden=f;
    CLsb.push_back(2*llnum-2*llden);

    // CL_{b}
    data_=CLB_pdata;
    data2_=CLB_pdata2;
    data3_=CLB_pdata3;
    pars[poiIndex_]=poiVal;
    nll(a, 0, f, pars, 0);
    llnum=f;
    pars[poiIndex_]=0.0;
    nll(a, 0, f, pars, 0);
    llden=f;
    CLb.push_back(2*llnum-2*llden);

    // remove the pseudodata
    if(i>0) {
      delete CLSB_pdata;
      delete CLB_pdata;
    }
  }

  // remove nuisance parameter priors
  for(std::map<int, RandomPrior*>::const_iterator it=priors.begin(); it!=priors.end(); ++it) delete it->second;

  // set the data back
  data_=theData;
  data2_=theData2;
  data3_=theData3;

  // calculate the CLs
  int nNum=0, nDen=0;
  for(unsigned int i=1; i<CLb.size(); i++) {

    bool numPass=(CLsb[i]>=CLsb[0]);
    bool denPass=(CLb[i]>CLb[0]);

    nNum+=numPass;
    nDen+=denPass;
  }
  return std::pair<int, int>(nNum, nDen);
}
