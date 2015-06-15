/*
   - Author : Yeng-Ming (Jacky) Tzeng
   - e-mail : ymtzeng@cern.ch
   - create time : 11/09/2014
   - version : 0.0

    Usage :
    1). source root environment
    2). g++ -o trim_posterior trim_posterior.cc `root-config --cflags` `root-config --glibs`

    Solutions to problems
    1). Postierior : bad with negative signal strength
         com_limit_600_0.1_12114.log:expected bound = [ 0.000000 , 0.029427 ] 
         nuisance/600_0.1_Combined_12114.root

        solution : require no negative        

    2). Bad fit on muon channel : due to bad PE (the most of bins almost close to zero)
        com_limit_1000_1.0_12123.log:expected bound = [ 0.000000 , 0.649900 ]
        nuisance/1000_1_Combined_12123.root

         - solution : the first bin should be larger than 10


   */


#include <iostream>
#include <utility>
#include <cmath>
#include <stdlib.h>
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"

using namespace std;

// alpha (1-alpha=confidence interval)
const double ALPHA=0.05;
// left side tail
const double LEFTSIDETAIL=0.0;
char buffer[256];
std::pair<double, double> evaluateInterval(TGraph* posterior, double alpha, double leftsidetail);

int main(int argc, char *argv[]){

    if(argc<5) {
        cout << "Usage: trim_posterior [file] [mass] [BR] [Rnd] [Npe]" << endl;
        cout << "Example: trim_posterior nuisance/1000_0.5_Combined_123.root 1000 0.5 123 5" << endl;
        return 0;
    }

    string filename = argv[1];
    int mass = atoi(argv[2]);
    float BR = atof(argv[3]);   // BR for Tgluon    ; 0 : 100% Tgamma
    int Rnd  = atoi(argv[4]);
    int Npe  = atoi(argv[5]);

    TFile *file = new TFile(filename.c_str());

    if(true){

        // setup the limit values
        double observedLowerBound, observedUpperBound;
        vector<double> expectedLowerBounds;
        vector<double> expectedUpperBounds;

        // evaluate the limit
        TGraph *post_data = (TGraph*) file->Get("post_0");
        pair<double, double> bounds_data=evaluateInterval(post_data, ALPHA, LEFTSIDETAIL);
        observedLowerBound=bounds_data.first;
        observedUpperBound=bounds_data.second;

        FILE * limitFile;
        sprintf(buffer,"LOG_Trim/com_limit_%i_%1.1f_%i.log",mass,BR,Rnd);
        limitFile = fopen (buffer,"w");

        // evaluate the limit
        for(int ipe=1;ipe<=Npe;ipe++){
            sprintf(buffer,"post_%i",ipe);
            TGraph *post = (TGraph*) file->Get(buffer);
            if(!post)
                continue;
            pair<double, double> bounds=evaluateInterval(post, ALPHA, LEFTSIDETAIL);

            // solution1 to probelm 1
            bool IsNegativeSignalStrength=false;
            for(int ipost = 1; ipost <= post->GetN(); ipost++){
                double xpost = 999;
                double ypost = 999;
                post->GetPoint(ipost, xpost, ypost); 
                if(xpost<0)
                    IsNegativeSignalStrength = true;
            }

            // solution 2 to problem 2
            sprintf(buffer,"data_%i",ipe);  // muon channel in Tgloun
            TH1D *PE_data = (TH1D*) file->Get(buffer);
            bool IsBadPE = false;
            if(PE_data->GetBinContent(1) < 10) 
                IsBadPE = true;

            sprintf(buffer,"data2_%i",ipe); // electron channel in Tgloun
            TH1D *PE2_data = (TH1D*) file->Get(buffer);
            if(PE2_data->GetBinContent(1) < 10) 
                IsBadPE = true;

            if (BR==0)  // exception for 100% Tgamma channel
                IsBadPE = false;

            if(bounds.first==0. && bounds.second>0. && !IsNegativeSignalStrength && !IsBadPE)
            {
                expectedLowerBounds.push_back(bounds.first);
                expectedUpperBounds.push_back(bounds.second);
            }
            delete post;
            delete PE_data;
            delete PE2_data;
        }

        for(unsigned int i=0; i<expectedLowerBounds.size(); i++) {
            cout << "expected bound(" << (i+1) << ") = [ " << expectedLowerBounds[i] << " , " << expectedUpperBounds[i] << " ]" << endl;
            fprintf(limitFile, "expected bound = [ %1.6f , %1.6f ] \n", expectedLowerBounds[i], expectedUpperBounds[i]);
        }

        cout << "\nobserved bound = [ " << observedLowerBound << " , " << observedUpperBound << " ]" << endl;

        fprintf(limitFile, "observed bound = [ %1.6f , %1.6f ] \n", observedLowerBound, observedUpperBound);
        fclose(limitFile);
    }

    delete file;

    return 0;
}

std::pair<double, double> evaluateInterval(TGraph* posterior, double alpha, double leftsidetail)
{
    double lowerCutOff = leftsidetail * alpha;
    double upperCutOff = 1. - (1.- leftsidetail) * alpha;

    double upper = 0, lower = 0;

    // normalize the interval, first
    double normalization=0.0;
    for(int i=0; i<posterior->GetN()-1; i++) {
        double firstx, firsty;
        double nextx, nexty;
        posterior->GetPoint(i, firstx, firsty);
        posterior->GetPoint(i+1, nextx, nexty);

        double intervalIntegral=(nextx-firstx)*0.5*(firsty+nexty);
        normalization+=intervalIntegral;
    }
    // now compute the intervals
    double integral=0.0;
    for(int i=0; i<posterior->GetN()-1; i++) {
        double firstx, firsty;
        double nextx, nexty;
        posterior->GetPoint(i, firstx, firsty);
        posterior->GetPoint(i+1, nextx, nexty);

        double intervalIntegral=(nextx-firstx)*0.5*(firsty+nexty)/normalization;
        // interpolate lower
        if(integral<=lowerCutOff && (integral+intervalIntegral)>=lowerCutOff) {
            lower=firstx;
        }
        if(integral<=upperCutOff && (integral+intervalIntegral)>=upperCutOff) {
            double m=(nexty-firsty)/(nextx-firstx);
            upper = firstx+(-firsty+sqrt(firsty*firsty+2*m*(upperCutOff-integral)*normalization))/m;
        }
        integral+=intervalIntegral;
    }

    std::pair<double, double> p(lower, upper);
    return p;
}

