float BranchingRatio = 1.0;
#include "TFile.h"
#include "TH2F.h"
#include "math.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TBox.h"
#include "TPolyLine.h"
#include "TStyle.h"
#include "TSpline.h"
#include <iostream>


const int NxsecPlot = 49;
double pythia_xsec_xPlot[NxsecPlot] = {
    300.,
    325.,
    350.,
    375.,
    400.,
    425.,
    450.,
    475.,
    500.,
    525.,
    550.,
    575.,
    600.,
    625.,
    650.,
    675.,
    700.,
    725.,
    750.,
    775.,
    800.,
    825.,
    850.,
    875.,
    900.,
    925.,
    950.,
    975.,
    1000.,
    1025.,
    1050.,
    1075.,
    1100.,
    1125.,
    1150.,
    1175.,
    1200.,
    1225.,
    1250.,
    1275.,
    1300.,
    1325.,
    1350.,
    1375.,
    1400.,
    1425.,
    1450.,
    1475.,
    1500.
};
double nlo_xsec_yPlot[NxsecPlot] = {
    298.834, //300. GeV
    159.382, //325. GeV
    88.8963, //350. GeV
    51.5036, //375. GeV
    30.83, //400. GeV
    18.9847, //425. GeV
    11.9843, //450. GeV
    7.7328, //475. GeV
    5.0875, //500. GeV
    3.40571, //525. GeV
    2.31567, //550. GeV
    1.59673, //575. GeV
    1.1151, //600. GeV
    0.787801, //625. GeV
    0.562451, //650. GeV
    0.405446, //675. GeV
    0.294853, //700. GeV
    0.216161, //725. GeV
    0.159653, //750. GeV
    0.118728, //775. GeV
    0.088853, //800. GeV
    0.0668855, //825. GeV
    0.0506221, //850. GeV
    0.0385055, //875. GeV
    0.0294248, //900. GeV
    0.0225825, //925. GeV
    0.0174004, //950. GeV
    0.013457, //975. GeV
    0.0104431, //1000. GeV
    0.00812992, //1025. GeV
    0.00634794, //1050. GeV
    0.0049701, //1075. GeV
    0.00390116, //1100. GeV
    0.00306934, //1125. GeV
    0.00242006, //1150. GeV
    0.00191195, //1175. GeV
    0.00151328, //1200. GeV
    0.00119975, //1225. GeV
    0.000952637, //1250. GeV
    0.000757473, //1275. GeV
    0.00060306, //1300. GeV
    0.000480664, //1325. GeV
    0.000383501, //1350. GeV
    0.000306258, //1375. GeV
    0.000244766, //1400. GeV
    0.000195757, //1425. GeV
    0.000156656, //1450. GeV
    0.000125425, //1475. GeV
    0.00010046 //1500. GeV
}; 

const int Nxsec = 12+1;
double pythia_xsec_x[Nxsec] = {
    600.,
    650.,
    700.,
    750.,
    800.,
    850.,
    900.,
    950.,
    1000.,
    1100.,
    1200.,
    1300.,
    1400.,
    //1500.
};

double calimits(double *S95,double *xsec ){
    double exclude = 0.;
    for(int idx=0;idx<Nxsec-1;idx++) {
        if ((S95[idx]-xsec[idx])*(S95[idx+1]-xsec[idx+1])>0.) continue;
        exclude = (0.-log(S95[idx]/xsec[idx]))/(log(S95[idx+1]/xsec[idx+1])-log(S95[idx]/xsec[idx]))*(pythia_xsec_x[idx+1]-pythia_xsec_x[idx]) + pythia_xsec_x[idx];
    } 
    // Using an extrapolation
    /*
    if (exclude==0){
        // L1 (theoritical) -> logy = a1 X + b1
        double a1 = (log(xsec[Nxsec-2]) - log(xsec[Nxsec-1]))/(pythia_xsec_x[Nxsec-2] - pythia_xsec_x[Nxsec-1]);
        double b1 = (-1*pythia_xsec_x[Nxsec-1]* log(xsec[Nxsec-2]) + pythia_xsec_x[Nxsec-2]*log(xsec[Nxsec-1]))/(pythia_xsec_x[Nxsec-2] - pythia_xsec_x[Nxsec-1]);

        // L2 (expected) -> logy = a2 X + b2
        double a2 = (log(S95[Nxsec-2]) - log(S95[Nxsec-1]))/(pythia_xsec_x[Nxsec-2] - pythia_xsec_x[Nxsec-1]);
        double b2 = (-1*pythia_xsec_x[Nxsec-1]* log(S95[Nxsec-2]) + pythia_xsec_x[Nxsec-2]*log(S95[Nxsec-1]))/(pythia_xsec_x[Nxsec-2] - pythia_xsec_x[Nxsec-1]);
        exclude = (b2 - b1)/(a1 - a2);
    }   
    */
    if (exclude==0&&(S95[Nxsec-1]<=xsec[Nxsec-1]&&S95[Nxsec-2]<=xsec[Nxsec-2])){
        // L1 (theoritical) -> logy = a1 X + b1
        double a1 = (log(xsec[Nxsec-2]) - log(xsec[Nxsec-1]))/(pythia_xsec_x[Nxsec-2] - pythia_xsec_x[Nxsec-1]);
        double b1 = (-1*pythia_xsec_x[Nxsec-1]* log(xsec[Nxsec-2]) + pythia_xsec_x[Nxsec-2]*log(xsec[Nxsec-1]))/(pythia_xsec_x[Nxsec-2] - pythia_xsec_x[Nxsec-1]);

        // L2 (expected) -> logy = a2 X + b2
        double a2 = (log(S95[Nxsec-2]) - log(S95[Nxsec-1]))/(pythia_xsec_x[Nxsec-2] - pythia_xsec_x[Nxsec-1]);
        double b2 = (-1*pythia_xsec_x[Nxsec-1]* log(S95[Nxsec-2]) + pythia_xsec_x[Nxsec-2]*log(S95[Nxsec-1]))/(pythia_xsec_x[Nxsec-2] - pythia_xsec_x[Nxsec-1]);
        exclude = (b2 - b1)/(a1 - a2);
    }else if (exclude==0&&(S95[0]>=xsec[0]&&S95[1]>=xsec[1])){
        // xsec <-> 1
        // Xsec_x <-> pythia_xsec_x
        // S95 <-> S95
        // L1 (theoritical) -> logy = a1 X + b1
        double a1 = (log(xsec[0]) - log(xsec[1]))/(pythia_xsec_x[0] - pythia_xsec_x[1]);
        double b1 = (-1*pythia_xsec_x[1]* log(xsec[0]) + pythia_xsec_x[0]*log(xsec[1]))/(pythia_xsec_x[0] - pythia_xsec_x[1]);

        // L2 (expected) -> logy = a2 X + b2
        double a2 = (log(S95[0]) - log(S95[1]))/(pythia_xsec_x[0] - pythia_xsec_x[1]);
        double b2 = (-1*pythia_xsec_x[1]* log(S95[0]) + pythia_xsec_x[0]*log(S95[1]))/(pythia_xsec_x[0] - pythia_xsec_x[1]);
        exclude = (b2 - b1)/(a1 - a2); 
    }    


    return exclude;
}

double nlo_xsec_y[Nxsec] = {
    1.1151, //600. GeV
    0.562451, //650. GeV
    0.294853, //700. GeV
    0.159653, //750. GeV
    0.088853, //800. GeV
    0.0506221, //850. GeV
    0.0294248, //900. GeV
    0.0174004, //950. GeV
    0.0104431, //1000. GeV
    0.00390116, //1100. GeV
    0.00151328, //1200. GeV
    0.00060306, //1300. GeV
    0.000244766, //1400. GeV
    //0.00010046 //1500. GeV
}; 

// Scale errors include renormalization and factorization (https://twiki.cern.ch/twiki/bin/view/Sandbox/CrossSectionsCalculationTool)
/*
double nlo_xsec_y_perr_scal[17] = {
8.09,
5.06,
3.25,
2.13,
1.44,
0.68,
0.34,
0.178,
0.097
};
*/

const int NMassPoints = 12+1;
double limitmid_xsec_x[NMassPoints] = {
    600,
    650,
    700,
    750,
    800,
    850,
    900,
    950,
    1000,
    1100,
    1200,
    1300,
    1400,
    //1500
};     

double limit_data_xsec_y[NMassPoints] = {
// Observed
0.008675, //Bayesian_limit_Fit_datacard_Tgamma600GeV_1.0.txt.log
0.016745, //Bayesian_limit_Fit_datacard_Tgamma650GeV_1.0.txt.log
0.030505, //Bayesian_limit_Fit_datacard_Tgamma700GeV_1.0.txt.log
0.052181, //Bayesian_limit_Fit_datacard_Tgamma750GeV_1.0.txt.log
0.0909105, //Bayesian_limit_Fit_datacard_Tgamma800GeV_1.0.txt.log
0.14799, //Bayesian_limit_Fit_datacard_Tgamma850GeV_1.0.txt.log
0.250258, //Bayesian_limit_Fit_datacard_Tgamma900GeV_1.0.txt.log
0.408365, //Bayesian_limit_Fit_datacard_Tgamma950GeV_1.0.txt.log
0.710698, //Bayesian_limit_Fit_datacard_Tgamma1000GeV_1.0.txt.log
1.80344, //Bayesian_limit_Fit_datacard_Tgamma1100GeV_1.0.txt.log
4.89664, //Bayesian_limit_Fit_datacard_Tgamma1200GeV_1.0.txt.log
12.8046, //Bayesian_limit_Fit_datacard_Tgamma1300GeV_1.0.txt.log
33.192, //Bayesian_limit_Fit_datacard_Tgamma1400GeV_1.0.txt.log
}; 

double limit_m2sig_xsec_y[NMassPoints] = {

// 2.5%
0.00254055, //Bayesian_limit_Fit_datacard_Tgamma600GeV_1.0.txt.log
0.00496384, //Bayesian_limit_Fit_datacard_Tgamma650GeV_1.0.txt.log
0.00993615, //Bayesian_limit_Fit_datacard_Tgamma700GeV_1.0.txt.log
0.0178206, //Bayesian_limit_Fit_datacard_Tgamma750GeV_1.0.txt.log
0.0324785, //Bayesian_limit_Fit_datacard_Tgamma800GeV_1.0.txt.log
0.0573803, //Bayesian_limit_Fit_datacard_Tgamma850GeV_1.0.txt.log
0.101018, //Bayesian_limit_Fit_datacard_Tgamma900GeV_1.0.txt.log
0.170919, //Bayesian_limit_Fit_datacard_Tgamma950GeV_1.0.txt.log
0.309391, //Bayesian_limit_Fit_datacard_Tgamma1000GeV_1.0.txt.log
0.841573, //Bayesian_limit_Fit_datacard_Tgamma1100GeV_1.0.txt.log
2.33447, //Bayesian_limit_Fit_datacard_Tgamma1200GeV_1.0.txt.log
6.55901, //Bayesian_limit_Fit_datacard_Tgamma1300GeV_1.0.txt.log
17.366, //Bayesian_limit_Fit_datacard_Tgamma1400GeV_1.0.txt.log
}; 
double limit_m1sig_xsec_y[NMassPoints] = {
// 16.0%
0.00320169, //Bayesian_limit_Fit_datacard_Tgamma600GeV_1.0.txt.log
0.00584216, //Bayesian_limit_Fit_datacard_Tgamma650GeV_1.0.txt.log
0.0112021, //Bayesian_limit_Fit_datacard_Tgamma700GeV_1.0.txt.log
0.0196855, //Bayesian_limit_Fit_datacard_Tgamma750GeV_1.0.txt.log
0.0352851, //Bayesian_limit_Fit_datacard_Tgamma800GeV_1.0.txt.log
0.0623161, //Bayesian_limit_Fit_datacard_Tgamma850GeV_1.0.txt.log
0.108376, //Bayesian_limit_Fit_datacard_Tgamma900GeV_1.0.txt.log
0.179253, //Bayesian_limit_Fit_datacard_Tgamma950GeV_1.0.txt.log
0.326408, //Bayesian_limit_Fit_datacard_Tgamma1000GeV_1.0.txt.log
0.86984, //Bayesian_limit_Fit_datacard_Tgamma1100GeV_1.0.txt.log
2.41829, //Bayesian_limit_Fit_datacard_Tgamma1200GeV_1.0.txt.log
6.78211, //Bayesian_limit_Fit_datacard_Tgamma1300GeV_1.0.txt.log
17.8959, //Bayesian_limit_Fit_datacard_Tgamma1400GeV_1.0.txt.log

};     

double limitmid_xsec_y[NMassPoints] = {
// 50.0%
0.005026, //Bayesian_limit_Fit_datacard_Tgamma600GeV_1.0.txt.log
0.008371, //Bayesian_limit_Fit_datacard_Tgamma650GeV_1.0.txt.log
0.014863, //Bayesian_limit_Fit_datacard_Tgamma700GeV_1.0.txt.log
0.025241, //Bayesian_limit_Fit_datacard_Tgamma750GeV_1.0.txt.log
0.046919, //Bayesian_limit_Fit_datacard_Tgamma800GeV_1.0.txt.log
0.078059, //Bayesian_limit_Fit_datacard_Tgamma850GeV_1.0.txt.log
0.132636, //Bayesian_limit_Fit_datacard_Tgamma900GeV_1.0.txt.log
0.215042, //Bayesian_limit_Fit_datacard_Tgamma950GeV_1.0.txt.log
0.387317, //Bayesian_limit_Fit_datacard_Tgamma1000GeV_1.0.txt.log
1.01393, //Bayesian_limit_Fit_datacard_Tgamma1100GeV_1.0.txt.log
2.79528, //Bayesian_limit_Fit_datacard_Tgamma1200GeV_1.0.txt.log
7.5844, //Bayesian_limit_Fit_datacard_Tgamma1300GeV_1.0.txt.log
19.8815, //Bayesian_limit_Fit_datacard_Tgamma1400GeV_1.0.txt.log
};

double limit_p1sig_xsec_y[NMassPoints] = {

// 84.0%
0.00747526, //Bayesian_limit_Fit_datacard_Tgamma600GeV_1.0.txt.log
0.0133331, //Bayesian_limit_Fit_datacard_Tgamma650GeV_1.0.txt.log
0.0223599, //Bayesian_limit_Fit_datacard_Tgamma700GeV_1.0.txt.log
0.0373546, //Bayesian_limit_Fit_datacard_Tgamma750GeV_1.0.txt.log
0.0707257, //Bayesian_limit_Fit_datacard_Tgamma800GeV_1.0.txt.log
0.119508, //Bayesian_limit_Fit_datacard_Tgamma850GeV_1.0.txt.log
0.191912, //Bayesian_limit_Fit_datacard_Tgamma900GeV_1.0.txt.log
0.309126, //Bayesian_limit_Fit_datacard_Tgamma950GeV_1.0.txt.log
0.546215, //Bayesian_limit_Fit_datacard_Tgamma1000GeV_1.0.txt.log
1.39321, //Bayesian_limit_Fit_datacard_Tgamma1100GeV_1.0.txt.log
4.04633, //Bayesian_limit_Fit_datacard_Tgamma1200GeV_1.0.txt.log
10.9721, //Bayesian_limit_Fit_datacard_Tgamma1300GeV_1.0.txt.log
29.3467, //Bayesian_limit_Fit_datacard_Tgamma1400GeV_1.0.txt.log
};     

double limit_p2sig_xsec_y[NMassPoints] = {

// 97.5%
0.0123398, //Bayesian_limit_Fit_datacard_Tgamma600GeV_1.0.txt.log
0.0221078, //Bayesian_limit_Fit_datacard_Tgamma650GeV_1.0.txt.log
0.0325536, //Bayesian_limit_Fit_datacard_Tgamma700GeV_1.0.txt.log
0.0572421, //Bayesian_limit_Fit_datacard_Tgamma750GeV_1.0.txt.log
0.105403, //Bayesian_limit_Fit_datacard_Tgamma800GeV_1.0.txt.log
0.167538, //Bayesian_limit_Fit_datacard_Tgamma850GeV_1.0.txt.log
0.27909, //Bayesian_limit_Fit_datacard_Tgamma900GeV_1.0.txt.log
0.450859, //Bayesian_limit_Fit_datacard_Tgamma950GeV_1.0.txt.log
0.764974, //Bayesian_limit_Fit_datacard_Tgamma1000GeV_1.0.txt.log
1.97288, //Bayesian_limit_Fit_datacard_Tgamma1100GeV_1.0.txt.log
5.37964, //Bayesian_limit_Fit_datacard_Tgamma1200GeV_1.0.txt.log
17.1737, //Bayesian_limit_Fit_datacard_Tgamma1300GeV_1.0.txt.log
43.7423, //Bayesian_limit_Fit_datacard_Tgamma1400GeV_1.0.txt.log

}; 




void show_limit_combinedFit()
{

    char buffer[128];
    for(int i=0;i<NMassPoints;i++){

        limit_data_xsec_y[i]*=nlo_xsec_y[i];
        limitmid_xsec_y[i]*=nlo_xsec_y[i];
        limit_p1sig_xsec_y[i]*=nlo_xsec_y[i];
        limit_p2sig_xsec_y[i]*=nlo_xsec_y[i];
        limit_m1sig_xsec_y[i]*=nlo_xsec_y[i];
        limit_m2sig_xsec_y[i]*=nlo_xsec_y[i];

    }
    double box2sig_xsec_x[2*NMassPoints] ;
    double box2sig_xsec_y[2*NMassPoints] ;
    double box1sig_xsec_x[2*NMassPoints] ;
    double box1sig_xsec_y[2*NMassPoints] ;
    for(int i=0;i<2*NMassPoints;i++)
    {
        if(i<NMassPoints){
            box2sig_xsec_x[i] = limitmid_xsec_x[i];
            box2sig_xsec_y[i] = limit_p2sig_xsec_y[i];
            box1sig_xsec_x[i] = limitmid_xsec_x[i];
            box1sig_xsec_y[i] = limit_p1sig_xsec_y[i];
        }else{
            box2sig_xsec_x[i] = limitmid_xsec_x[2*NMassPoints-1-i];
            box2sig_xsec_y[i] = limit_m2sig_xsec_y[2*NMassPoints-1-i];
            box1sig_xsec_x[i] = limitmid_xsec_x[2*NMassPoints-1-i];
            box1sig_xsec_y[i] = limit_m1sig_xsec_y[2*NMassPoints-1-i];
        }
    };     

    gROOT->ProcessLine(".L setTDRStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

   //std::cout<<"expected : "<<calimits(limitmid_xsec_y,nlo_xsec_y)<<std::endl;
   //std::cout<<"expected upper : "<<calimits(limit_m2sig_xsec_y,nlo_xsec_y)<<std::endl;
   //std::cout<<"expected lower : "<<calimits(limit_p2sig_xsec_y,nlo_xsec_y)<<std::endl;
   //std::cout<<"observed : "<<calimits(limit_data_xsec_y,nlo_xsec_y)<<std::endl;

   //TH2F *frame = new TH2F("frame","frame",14,230,570,10,0.05,50.); 
   //TH2F *frame = new TH2F("frame","frame",14,470,830,10,0.01,20.); 

   //TH2F *frame = new TH2F("frame","frame",18,430,980,10,0.005,4.); 
   //TH2F *frame = new TH2F("frame","frame",18,430,980,10,nlo_xsec_y[Nxsec-1],nlo_xsec_y[0]); 
   //TH2F *frame = new TH2F("frame","frame",18,570,1330,20,0.0001,10.); 
   //TH2F *frame = new TH2F("frame","frame",18,620,1430,20,0.0001,10.); 
   TH2F *frame = new TH2F("frame","frame",18,620-50,1430+100-100,20,0.0001,0.5*10); 

   TCanvas *c1 = new TCanvas("c1","c1",0,0,700,500);
   c1->SetFrameFillColor(10);
   c1->SetFrameBorderMode(0);
   c1->SetFillColor(10);
   //c1->SetFillStyle(4000);
   c1->SetBorderMode(0);
   c1->SetFrameLineWidth(3);
   c1->SetBottomMargin(0.17);
   c1->SetRightMargin(0.05);
   c1->SetLeftMargin(0.14);
   c1->SetTopMargin(0.08);

   frame->GetYaxis()->SetLabelSize(0.06);
   frame->GetYaxis()->SetLabelOffset(0.01);
   frame->GetYaxis()->SetTitleSize(0.07);
   frame->GetYaxis()->SetTitleOffset(0.90);
   frame->GetYaxis()->SetTitle("#sigma(pp #rightarrow t*#bar{t*}X) [pb]");

   frame->GetXaxis()->SetLabelSize(0.06);
   frame->GetXaxis()->SetLabelOffset(0.01);
   frame->GetXaxis()->SetTitleSize(0.07);
   frame->GetXaxis()->SetTitleOffset(0.99);
   frame->GetXaxis()->SetTitle("M_{t*} [GeV/c^{2}]");
   frame->SetStats(kFALSE);
   frame->SetTitle("");

   frame->SetLineWidth(3);
   frame->SetMarkerStyle(21);		
   frame->Draw("");	


   c1->SetLogy();	
   //c1->SetGridx();		
   //c1->SetGridy();

   gStyle->SetHatchesSpacing(2);

   //box1
   TPolyLine *bx2 = new TPolyLine(2*NMassPoints,box2sig_xsec_x,box2sig_xsec_y);
   bx2->SetFillColor(kGreen+2);
   bx2->SetFillStyle(3144);
   bx2->Draw("f");

   //box2
   TPolyLine *bx3 = new TPolyLine(2*NMassPoints,box1sig_xsec_x,box1sig_xsec_y);
   bx3->SetFillColor(kYellow-3);
   bx3->SetFillStyle(1001);
   bx3->Draw("f");

   TSpline3 *sp5 = new TSpline3("sp5",pythia_xsec_xPlot,nlo_xsec_yPlot,NxsecPlot);
	sp5->SetLineWidth(4);
	sp5->SetLineStyle(5);
   sp5->SetLineColor(593);
   sp5->Draw("same");

   TPolyLine *sp2 = new TPolyLine(NMassPoints,limitmid_xsec_x,limitmid_xsec_y);
	sp2->SetLineWidth(4);
	sp2->SetLineStyle(7);
   sp2->SetLineColor(kRed-2);
   sp2->Draw("same");

   TPolyLine *sp6 = new TPolyLine(NMassPoints,limitmid_xsec_x,limit_data_xsec_y);
   sp6->SetLineWidth(4);
   sp6->SetLineColor(kBlack);
   sp6->Draw("same");

   //TLatex *tex1 = new TLatex(420,0.55,"NLO #sqrt{s}=7 TeV (Berger and Cao)");
   //TLatex *tex1 = new TLatex(420,0.55,"NLO (Berger and Cao)");
   //TLatex *tex1 = new TLatex(420,0.55,"NNLO (HATHOR)");
   //TLatex *tex1 = new TLatex(420,0.55,"#sigma(T#bar{T}) @ NNLO [25]");
   //TLatex *tex1 = new TLatex(265.9202,8.839836,"Prediction");
   TLatex *tex1 = new TLatex(700,0.151,"Spin 3/2 t*");
   tex1->SetTextSize(0.055);
   tex1->SetTextFont(22);
   tex1->SetTextColor(593);
   //tex1->SetTextAngle(-16);
   tex1->SetTextAngle(322.369);
   //tex1->Draw();	

	tex1 = new TLatex(475,limitmid_xsec_y[1]*1.2,"expected limit");
	tex1->SetTextSize(0.055);
	tex1->SetTextFont(22);
	tex1->SetTextColor(kRed-2);
	tex1->SetTextAngle(-1);
	//tex1->Draw();	
	
	tex1 = new TLatex(475,limit_data_xsec_y[1]*0.6,"observed limit");
	tex1->SetTextSize(0.055);
	tex1->SetTextFont(22);
	tex1->SetTextColor(kBlack);
	tex1->SetTextAngle(-1);
	//tex1->Draw();
	
	//TLatex *tex1 = new TLatex(fabs(limitmid_xsec_x[0]-limitmid_xsec_x[NMassPoints-1])*0.12+300,0.11,"Limit at 95% CL: M_{t*} > 435 GeV/c^{2}");
    sprintf(buffer,"Expected limit at 95%% CL: M_{t*} > %1.f^{+%1.f}_{%1.f} GeV/c^{2}",
            calimits(limitmid_xsec_y,nlo_xsec_y),
            calimits(limit_m2sig_xsec_y,nlo_xsec_y) - calimits(limitmid_xsec_y,nlo_xsec_y),
            calimits(limit_p2sig_xsec_y,nlo_xsec_y) - calimits(limitmid_xsec_y,nlo_xsec_y)
            );
	tex1 = new TLatex(fabs(limitmid_xsec_x[0]-limitmid_xsec_x[NMassPoints-1])*0.12+550+50-50,0.00018,buffer);
	tex1->SetTextSize(0.055);
	tex1->SetTextFont(22);
	tex1->SetTextColor(kGreen-2);
	tex1->Draw();	

    sprintf(buffer,"Observed limit at 95%% CL: M_{t*} > %1.f GeV/c^{2}",calimits(limit_data_xsec_y,nlo_xsec_y));
	tex1 = new TLatex(fabs(limitmid_xsec_x[0]-limitmid_xsec_x[NMassPoints-1])*0.12+550+50-50,0.0004,buffer);
	tex1->SetTextSize(0.055);
	tex1->SetTextFont(22);
	tex1->SetTextColor(kGreen-2);
	tex1->Draw();	
	
	tex1 = new TLatex(5.4*fabs(limitmid_xsec_x[0]-limitmid_xsec_x[NMassPoints-1])/6.+limitmid_xsec_x[0],45*0.4*0.4*1.05,"2 #sigma");
	tex1->SetTextSize(0.050);
	tex1->SetTextFont(22);
	tex1->SetTextColor(kBlack);
	//tex1->Draw();	

	tex1 = new TLatex(5.4*fabs(limitmid_xsec_x[0]-limitmid_xsec_x[NMassPoints-1])/6.+limitmid_xsec_x[0],0.4*0.2*1.05,"1 #sigma");
	tex1->SetTextSize(0.050);
	tex1->SetTextFont(22);
	tex1->SetTextColor(kBlack);
	//tex1->Draw();
	
	TBox *bx1_1 = new TBox(4.9*fabs(limitmid_xsec_x[0]-limitmid_xsec_x[NMassPoints-1])/6.+limitmid_xsec_x[0],
            0.4*0.4*1.5,
            5.2*fabs(limitmid_xsec_x[0]-limitmid_xsec_x[NMassPoints-1])/6.+limitmid_xsec_x[0],
            0.4*0.4);
	bx1_1->SetFillColor(kGreen+2);
	bx1_1->SetFillStyle(3144);
	//bx1->Draw();
	
	TBox *bx1_2 = new TBox(4.9*fabs(limitmid_xsec_x[0]-limitmid_xsec_x[NMassPoints-1])/6.+limitmid_xsec_x[0],
            0.4*0.2*1.5,
            5.2*fabs(limitmid_xsec_x[0]-limitmid_xsec_x[NMassPoints-1])/6.+limitmid_xsec_x[0],
            0.4*0.2);
	bx1_2->SetFillColor(kYellow-3);
	bx1_2->SetFillStyle(1001);
	//bx2->Draw();

	
    if(BranchingRatio == 0.0){
        sprintf(buffer,"CMS, L = 19.5 fb^{-1}, #sqrt{s} = 8 TeV");
        tex1 = new TLatex((limitmid_xsec_x[0]+limitmid_xsec_x[NMassPoints-1]-100)/2.,0.5*1.15*10,buffer);
    }else if (BranchingRatio != 0.0 && BranchingRatio != 1){
        sprintf(buffer,"CMS, L = 19.7 (19.5) fb^{-1}, #sqrt{s} = 8 TeV");
        tex1 = new TLatex((limitmid_xsec_x[0]+limitmid_xsec_x[NMassPoints-1]-100)*4./9.,0.5*1.15*10,buffer);
    }else{
        sprintf(buffer,"CMS, L = 19.7 fb^{-1}, #sqrt{s} = 8 TeV");
        tex1 = new TLatex((limitmid_xsec_x[0]+limitmid_xsec_x[NMassPoints-1]-100)/2.,0.5*1.15*10,buffer);
    }

	tex1->SetTextSize(0.055);
	tex1->Draw();	

    TSpline3 *spE = new TSpline3("",pythia_xsec_xPlot,nlo_xsec_yPlot,NxsecPlot);

    TLegend *legend_nm = new TLegend(0.48,0.70,0.90+0.05,0.85);
    //legend_nm->AddEntry(sp2,"expected limit","l");
    //legend_nm->AddEntry(bx1_2,"1 #sigma","f");
    //legend_nm->AddEntry(sp6,"observed limit","l");
    //legend_nm->AddEntry(bx1_1,"2 #sigma","f");
    legend_nm->AddEntry(sp5,"Spin 3/2 t*","l");
    sprintf(buffer,"BR_{t*#rightarrow t#gamma}(%1.1f)",BranchingRatio);
    legend_nm->AddEntry(spE,buffer,"");
    legend_nm->AddEntry(sp6,"Observed limit","l");
    legend_nm->AddEntry(bx1_2,"68\% CL","f");
    //legend_nm->AddEntry(bx1_2,"1 #sigma","f");
    legend_nm->AddEntry(sp2,"Expected limit","l");
    //legend_nm->AddEntry(bx1_1,"2 #sigma","f");
    legend_nm->AddEntry(bx1_1,"95\% CL","f");


    legend_nm->SetBorderSize(0);
    legend_nm->SetFillColor(0);
    legend_nm->SetNColumns(2);
    legend_nm->SetTextSize(0.045);
    legend_nm->SetTextSizePixels(14);
    legend_nm->Draw();
    c1->SaveAs("LimitCombine.pdf");
    c1->SaveAs("LimitCombine.png");
}
