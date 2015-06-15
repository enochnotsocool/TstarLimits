#ifndef __BINNED_DATA_HH__
#define __BINNED_DATA_HH__

#include "fit.hh"

const double CriticalMassPoint = 350.;  // [CMP] : Tgluon excluding mass points before CMP, but Tgamma should include

class TH1D;

TH1D* getData(const std::vector<std::string>& filenames, const char* histname, int nbins, double* bins, std::string _dataname,int ModeForCMP/*0:exclude before CMP, 1: include before CMP*/);

TH1D* getSignalCDF(const char* filename1, const char* histname1, const char* filename2, const char* histname2, const double BR, const double eff_h, const double eff_l, double& _events,int ModeForCMP/*0:exclude before CMP, 1: include before CMP*/, const std::string& postfix = "");

#endif
