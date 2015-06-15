#ifndef __STATISTICS_HH__
#define __STATISTICS_HH__

#include <vector>

#include "fit.hh"

class TF1;
class TH1D;
class TGraph;

// lognormal distribution where par[0] is the median and par[1] is the variance
double lognormal(double *x, double *par);
// gaussian distribution where par[0] is the median and par[1] is the variance
double gaussian(double *x, double *par);
// gamma distribution where par[0] is the median and par[1] is the variance
double gamma(double *x, double *par);

// evaluate the interval given a posterior distribution, alpha, and left side tail
std::pair<double, double> evaluateInterval(TGraph* posterior, double alpha=0.05, double leftsidetail=0.0);

// get the quantiles for a set of numbers/limits
void getQuantiles(std::vector<double>& limits, double &median, std::pair<double, double>& onesigma, std::pair<double, double>& twosigma);

// get a random lognormal value given a value for the median and variance
class RandomPrior
{
public:
  RandomPrior(int priorType, double median, double variance, double min, double max);
  virtual ~RandomPrior();

  double getRandom(void) const;
  double getXmin(void) const;
  double getXmax(void) const;
  double eval(double x) const;
  int getPriorType(void) const;

private:

  TF1* priorfcn_;
  static int counter_;
  int priorType_;
};


#endif
