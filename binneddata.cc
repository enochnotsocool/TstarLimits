#include "binneddata.hh"
#include "fit.hh"

#include <cassert>
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1D.h>

TH1D* getData(const std::vector<std::string>& filenames, const char* histname, int nbins, double* bins, std::string _dataname, int ModeForCMP/*0:exclude before CMP, 1: include before CMP*/)
{
  //TH1D* hist=new TH1D("normalized_data", "normalized data", nbins, bins);
  TH1D* hist=new TH1D(_dataname.c_str(), _dataname.c_str(), nbins, bins);

  unsigned nFiles = filenames.size();

  TFile *files[nFiles];
  TH1D *hists[nFiles];

  // open input files and get input histograms
  for(unsigned i=0; i<nFiles; ++i)
  {
    files[i] = new TFile(filenames[i].c_str());
    hists[i] = (TH1D*)files[i]->Get(histname);
  }

  // normalize the data in the histogram
  for(int i=1; i<=hist->GetNbinsX(); ++i)
  {

    int bin_min = hists[0]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinLowEdge(i)+0.5);
    int bin_max = hists[0]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinUpEdge(i)-0.5);
    double val=hists[0]->Integral(bin_min,bin_max);
    if(hist->GetXaxis()->GetBinLowEdge(i)<CriticalMassPoint && ModeForCMP==0)
        val = 0.;
    double err=Fitter::histError(val);
    double width=hist->GetBinWidth(i);
    hist->SetBinContent(i, val/width);
    hist->SetBinError(i, err/width);

    // for simultaneous fit in 3 distribution (used in di-b-jet analysis, currently commented out)
    //if(hist->GetXaxis()->GetBinUpEdge(i)<=5890.)
    //{
    //  double offset = 0.;
    //  int bin_min = hists[0]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinLowEdge(i)-offset+0.5);
    //  int bin_max = hists[0]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinUpEdge(i)-offset-0.5);
    //  double val=hists[0]->Integral(bin_min,bin_max);
    //  double err=Fitter::histError(val);
    //  double width=hist->GetBinWidth(i);
    //  hist->SetBinContent(i, val/width);
    //  hist->SetBinError(i, err/width);
    //}
    //else if(hist->GetXaxis()->GetBinUpEdge(i)>5890. && hist->GetXaxis()->GetBinUpEdge(i)<=10890. && nFiles>1)
    //{
    //  double offset = 5000.;
    //  int bin_min = hists[1]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinLowEdge(i)-offset+0.5);
    //  int bin_max = hists[1]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinUpEdge(i)-offset-0.5);
    //  double val=hists[1]->Integral(bin_min,bin_max);
    //  double err=Fitter::histError(val);
    //  double width=hist->GetBinWidth(i);
    //  hist->SetBinContent(i, val/width);
    //  hist->SetBinError(i, err/width);
    //}
    //else if(hist->GetXaxis()->GetBinUpEdge(i)>10890. && nFiles>2)
    //{
    //  double offset = 10000.;
    //  int bin_min = hists[2]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinLowEdge(i)-offset+0.5);
    //  int bin_max = hists[2]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinUpEdge(i)-offset-0.5);
    //  double val=hists[2]->Integral(bin_min,bin_max);
    //  double err=Fitter::histError(val);
    //  double width=hist->GetBinWidth(i);
    //  hist->SetBinContent(i, val/width);
    //  hist->SetBinError(i, err/width);
    //}
  }

  // do some checks
  if(hist->GetBinContent(0)>0.0 || hist->GetBinContent(hist->GetNbinsX()+1)>0.0)
  {
    std::cerr << "There is an underflow/overflow.  This should be fixed.  Please rebin!" << std::endl;
    assert(0);
  }
  if(hist->GetBinContent(1)<=0.0)
  {
    std::cerr << "The first data bin is empty. This will give incorrect results.  Please rebin!" << std::endl;
    assert(0);
  }

  for(unsigned i=0; i<nFiles; ++i)
  {
    files[i]->Close();
  }
  
  return hist;
}


TH1D* getSignalCDF(const char* filename1, const char* histname1, const char* filename2, const char* histname2, const double BR, const double eff_h, const double eff_l, double& _events, 
        int ModeForCMP/*0:exclude before CMP, 1: include before CMP*/,
        const std::string& postfix 
        )
{
  TH1D* h_pdf=new TH1D("h_pdf", "Resonance Shape", 1550, 0, 1550);
  TH1D* h_cdf=new TH1D(("h_cdf"+postfix).c_str(), "Resonance Shape CDF", 1550, 0, 1550);


  TFile* histfile1=new TFile(filename1);
  TFile* histfile2=new TFile(filename2);

  TH1D* hist1=(TH1D*)histfile1->Get(histname1);
  TH1D* hist2=(TH1D*)histfile2->Get(histname2);


  hist1->Scale(eff_h*BR);
  hist2->Scale(eff_l*(1-BR));
  hist1->Add(hist2); // this adds together two resonance shapes with appropriate branching fractions and efficiencies
//  _events = hist1->Integral();
  double StartingMP = 50.;
  if(ModeForCMP==0)
      StartingMP = CriticalMassPoint;
  _events = hist1->Integral(hist1->GetXaxis()->FindBin(StartingMP),hist1->GetXaxis()->FindBin( 1550));

  for(int i=1; i<=hist1->GetNbinsX(); i++)
  {

    int bin_min = h_pdf->GetXaxis()->FindBin(hist1->GetXaxis()->GetBinLowEdge(i)+0.5);
    int bin_max = h_pdf->GetXaxis()->FindBin(hist1->GetXaxis()->GetBinUpEdge(i)-0.5);
    double bin_content = hist1->GetBinContent(i)/double(bin_max-bin_min+1);
    if(bin_content<0.) bin_content=0.; // some extrapolated resonance shapes can have a small negative bin content in some of the bins. this protects against such cases
    for(int b=bin_min; b<=bin_max; b++){
       h_pdf->SetBinContent(b, bin_content);
    }
  }


  h_pdf->Scale(1./h_pdf->Integral()); // normalize final PDF

  // produce CDF from the normalized PDF
  for(int i=1; i<=h_cdf->GetNbinsX(); i++)
  {
    int bin_min = h_pdf->GetXaxis()->FindBin(h_cdf->GetXaxis()->GetBinLowEdge(i)+0.5);
    int bin_max = h_pdf->GetXaxis()->FindBin(h_cdf->GetXaxis()->GetBinUpEdge(i)-0.5);

    double curr = 0;
    for(int b=bin_min; b<=bin_max; b++){
       curr+=h_pdf->GetBinContent(b);
    }

    double prev=h_cdf->GetBinContent(i-1);

    h_cdf->SetBinContent(i, prev+curr);
  }

  //printf("Ready read input from %s/%s (%f / %f / %f / %f)\n",filename1, histname1, hist1->Integral(), h_pdf->Integral(),
  //        h_cdf->Integral(), _events);

  histfile1->Close();
  histfile2->Close();
  delete h_pdf;

  return h_cdf;
}
