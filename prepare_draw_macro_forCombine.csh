#!/bin/csh

#Observed Limit: r < 1.3737
#Expected  2.5%: r < 0.1674
#Expected 16.0%: r < 0.2273
#Expected 50.0%: r < 0.3271
#Expected 84.0%: r < 0.4849
#Expected 97.5%: r < 0.7020

#set BRlist="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
#set BRlist="0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"

#set BRlist="0.0 0.3 0.5 0.7 0.8 0.9 1.0"
#set BRlist="0.0 0.1 0.2 0.3 0.5 0.7 1.0"
set BRlist="0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
#set BRlist="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"

#set BRlist="0.0 0.3 0.5 0.7 1.0"
#set limitValue=`printf "Observed 2.5%% 16.0%% 50.0%% 84.0%% 97.5%%"`
set limitValue=`printf "Observed Minus2sigma Minus1sigma Median Plus1sigma Plus2sigma"`

if ( -e templateX.cc ) then
rm templateX.cc
endif
cat .template1.cc >& templateX.cc

foreach BR(${BRlist})

if ( -e obs_limit_TPrime.txt ) then 
rm obs_limit_TPrime.txt
endif

if ( -e exp_limit_TPrime.txt ) then
rm exp_limit_TPrime.txt
endif

cat LIMIT_LOG/obs_limit_TPrime_*_${BR}.txt >& obs_limit_TPrime.txt
cat LIMIT_LOG/exp_limit_TPrime_*_${BR}.txt >& exp_limit_TPrime.txt

#set  mass_list=`cat obs_limit_TPrime.txt| grep -v x | awk '{print $1}' | sort -n| grep -v 1500`
#set  Nmass=`cat obs_limit_TPrime.txt| grep -v x | awk '{print $1}' | sort -n|wc|awk '{print $1-1}'`
set  mass_list=`cat obs_limit_TPrime.txt| grep -v x | awk '{print $1}' | sort -n`
set  Nmass=`cat obs_limit_TPrime.txt| grep -v x | awk '{print $1}' | sort -n|wc|awk '{print $1}'`

set br_=`echo ${BR} | sed 's/\./p/g'`
echo "const int N_${br_} = ${Nmass};"  >> templateX.cc

echo "// Nmass with BR(tgamma) = ${BR}"  >> templateX.cc
echo "double Xsec_x_${br_}[N_${br_}] = {"  >> templateX.cc
@ counter=1
foreach ml(${mass_list})
    if ( ${counter} != ${Nmass} ) then
        echo ${ml}"," >> templateX.cc
    else
        echo ${ml} >> templateX.cc
    endif
@ counter++
end
echo "};" >> templateX.cc
echo "" >> templateX.cc

@ counter=1
foreach lv(${limitValue})

    echo "// ${lv} with BR(tgamma) = ${BR}"  >> templateX.cc
    echo "double ${lv}_${br_}[N_${br_}] = {"  >> templateX.cc
    @ subcounter=1
    foreach ml(${mass_list})

    if ( ${subcounter} != ${Nmass} ) then
    paste obs_limit_TPrime.txt exp_limit_TPrime.txt  | grep -v x| awk '{print $2, $6, $8, $5, $9, $7, "Mass"$1}' | grep "Mass${ml}" | awk '{print $'${counter}' ", //'${ml}'"}'  >> templateX.cc
    else
    paste obs_limit_TPrime.txt exp_limit_TPrime.txt  | grep -v x| awk '{print $2, $6, $8, $5, $9, $7, "Mass"$1}' | grep "Mass${ml}" | awk '{print $'${counter}' " //'${ml}'"}'  >> templateX.cc
    endif

    @ subcounter++
    end

    echo "};" >> templateX.cc
    echo "" >> templateX.cc

@ counter++
end

end

rm *_limit_TPrime.txt

####### Part II ##########
#cat .template2.cc >> templateX.cc
echo "void templateX(){" >> templateX.cc

echo "// LimitsInfo[ BR, ${limitValue}]" >> templateX.cc
set NLimitsInfo=`echo "BR ${limitValue}" | awk '{print NF}'`
set NBRlist=`echo ${BRlist} | awk '{print NF}'`
echo "const int NLimitsInfo = ${NLimitsInfo};" >> templateX.cc
echo "const int NBRlist = ${NBRlist};" >> templateX.cc
echo "double LimitsInfo[NBRlist][NLimitsInfo] = {" >> templateX.cc

@ counterBR=1
foreach BR(${BRlist})

    set br_=`echo ${BR} | sed 's/\./p/g'`

    echo "{ ${BR}, // BR = ${BR}" >> templateX.cc
    @ counter=2
    foreach lv(${limitValue})
        #echo "std::cout<<"'"'" ${lv} : "'"'"<<calimits(${lv}_${br_}, Xsec_x_${br_},N_${br_})<<std::endl;" >> templateX.cc
        if ( ${counter} != ${NLimitsInfo} ) then
            echo "calimits(${lv}_${br_}, Xsec_x_${br_},N_${br_})," >> templateX.cc
        else 
            echo "calimits(${lv}_${br_}, Xsec_x_${br_},N_${br_})" >> templateX.cc
        endif
        @ counter++

    end

    if ( ${counterBR} != ${NBRlist} ) then
        echo "}," >> templateX.cc
    else
        echo "}" >> templateX.cc
    endif

    @ counterBR++
end
echo "};" >> templateX.cc

echo "TCanvas *c1 = new TCanvas("'"'"c1"'"'","'"'"c1"'"'",0,0,700,500);\
c1->SetFrameFillColor(10);\
c1->SetFrameBorderMode(0);\
c1->SetFillColor(10);\
c1->SetBorderMode(0);\
c1->SetFrameLineWidth(3);\
c1->SetBottomMargin(0.17);\
c1->SetRightMargin(0.05);\
c1->SetLeftMargin(0.14);\
c1->SetTopMargin(0.08);" >> templateX.cc
echo "    TGraphAsymmErrors *G_obs = new TGraphAsymmErrors(NBRlist);\
          for(int ibr=0;ibr<NBRlist;ibr++)\
          G_obs->SetPoint(ibr,LimitsInfo[ibr][1], LimitsInfo[ibr][0] );\
          //G_obs->Draw("'"'"APL"'"'"); " >>  templateX.cc 
echo "    TGraphAsymmErrors *G_exp = new TGraphAsymmErrors(NBRlist);\
         TGraphAsymmErrors *G_expBand = new TGraphAsymmErrors(NBRlist);\
          for(int ibr=0;ibr<NBRlist;ibr++){\
          G_exp->SetPoint(ibr,LimitsInfo[ibr][4], LimitsInfo[ibr][0] );\
          G_expBand->SetPoint(ibr,LimitsInfo[ibr][4], LimitsInfo[ibr][0] );\
          G_expBand->SetPointError(ibr, fabs(LimitsInfo[ibr][4]- LimitsInfo[ibr][6]), fabs(LimitsInfo[ibr][4]- LimitsInfo[ibr][2]),\
                               0.05,0.05);\
          }    \
          G_exp->SetLineWidth(4);\
          G_exp->SetLineStyle(7);\
          G_exp->SetLineColor(kRed-2);\
          G_expBand->GetYaxis()->SetRangeUser(0,1);\
          G_expBand->SetTitle("");\
          G_expBand->GetYaxis()->SetTitle("'"'"BR(t* #rightarrow t#gamma)"'"'");\
          G_expBand->GetXaxis()->SetLabelSize(0.06);\
          G_expBand->GetXaxis()->SetLabelOffset(0.01);\
          G_expBand->GetXaxis()->SetTitleSize(0.07);\
          G_expBand->GetXaxis()->SetTitleOffset(0.99);\
          G_expBand->GetYaxis()->SetLabelSize(0.06);\
          G_expBand->GetYaxis()->SetLabelOffset(0.01);\
          G_expBand->GetYaxis()->SetTitleSize(0.07);\
          G_expBand->GetYaxis()->SetTitleOffset(0.90);\
          G_expBand->GetXaxis()->SetTitle("'"'"Exclusion limit [GeV]"'"'");\
          G_expBand->SetFillColor(kWhite); \
          G_expBand->Draw("'"'"a2"'"'"); \
          // 2 sigma band\
          double box2sig_xsec_x[2*NBRlist];\
          double box2sig_xsec_y[2*NBRlist];\
          for(int i=0;i<NBRlist;i++){\
              box2sig_xsec_y[i] = LimitsInfo[i][0];\
              box2sig_xsec_x[i] = LimitsInfo[i][6];\
          }        \
for(int i=0;i<NBRlist;i++){\
    box2sig_xsec_y[i+NBRlist] = LimitsInfo[NBRlist-i-1][0];\
    box2sig_xsec_x[i+NBRlist] = LimitsInfo[NBRlist-i-1][2];\
}\
TPolyLine *bx2 = new TPolyLine(2*NBRlist,box2sig_xsec_x,box2sig_xsec_y);\
bx2->SetFillColor(kGreen+2);\
bx2->SetFillStyle(3144);\
bx2->Draw("'"'"f"'"'");\
// 1 sigma band\
double box1sig_xsec_x[2*NBRlist];\
double box1sig_xsec_y[2*NBRlist];\
for(int i=0;i<NBRlist;i++){\
    box1sig_xsec_y[i] = LimitsInfo[i][0];\
    box1sig_xsec_x[i] = LimitsInfo[i][5];\
}\
for(int i=0;i<NBRlist;i++){\
    box1sig_xsec_y[i+NBRlist] = LimitsInfo[NBRlist-i-1][0];\
    box1sig_xsec_x[i+NBRlist] = LimitsInfo[NBRlist-i-1][3];\
}\
TPolyLine *bx1 = new TPolyLine(2*NBRlist,box1sig_xsec_x,box1sig_xsec_y);\
bx1->SetFillColor(kYellow-3);\
bx1->SetFillStyle(1001);\
bx1->Draw("'"'"f"'"'");\
    // exclusion\
    double box3sig_xsec_x[2+NBRlist];\
    double box3sig_xsec_y[2+NBRlist];\
    for(int i=0;i<NBRlist;i++){\
        box3sig_xsec_y[i] = LimitsInfo[i][0];\
        box3sig_xsec_x[i] = LimitsInfo[i][1];\
    }\
box3sig_xsec_y[NBRlist] = LimitsInfo[NBRlist-1][0];\
box3sig_xsec_x[NBRlist] = G_expBand->GetXaxis()->GetXmin();\
box3sig_xsec_y[NBRlist+1] = LimitsInfo[0][0];\
box3sig_xsec_x[NBRlist+1] = G_expBand->GetXaxis()->GetXmin();\
TPolyLine *bx3 = new TPolyLine(2+NBRlist,box3sig_xsec_x,box3sig_xsec_y);\
bx3->SetFillColor(kGray+2);\
bx3->SetFillStyle(3004);\
bx3->Draw("'"'"f"'"'");\
          G_exp->Draw("'"'"L"'"'"); \
          \
            G_obs->SetLineWidth(4);\
          G_obs->SetLineColor(kBlack);\
          G_obs->Draw("'"'"PL, same"'"'"); "                >> templateX.cc

echo " TLatex *tex1 = new TLatex(850-50,1.03,"'"'"CMS 19.7 (19.5) fb^{-1}  #sqrt{s} = 8 TeV"'"'");\
          tex1->SetTextSize(0.055);\
                    tex1->Draw();">> templateX.cc

echo "            TLegend *legend_nm = new TLegend(0.157,0.605,0.5,0.92);\
                    legend_nm->AddEntry(G_obs,"'"'"Observed limit"'"'","'"'"l"'"'");\
                    legend_nm->AddEntry(G_exp,"'"'"Expected limit"'"'","'"'"l"'"'");\
                    //legend_nm->AddEntry(G_expBand,"'"'"2 sigma band"'"'","'"'"f"'"'");\
                    legend_nm->AddEntry(bx1,"'"'"1 sigma band"'"'","'"'"f"'"'");\
                    legend_nm->AddEntry(bx2,"'"'"2 sigma band"'"'","'"'"f"'"'");\
                    legend_nm->AddEntry(bx3,"'"'"Exclusion area (obs.)"'"'","'"'"f"'"'");\
                    legend_nm->SetBorderSize(0);\
                    legend_nm->SetFillColor(0);\
                    legend_nm->SetFillStyle(0);\
                    legend_nm->SetNColumns(1);\
                    legend_nm->SetTextSize(0.04);\
                    legend_nm->SetTextSizePixels(25);\
                    legend_nm->Draw();" >> templateX.cc

    echo "c1->SaveAs("'"'"CombinedResult.pdf"'"'");" >> templateX.cc

cat >> templateX.cc << EOF
printf("%%[Combination]\n");
for(int i=0;i<NBRlist;i++){
    printf("%1.1f & %1.0f & %1.0f & [%1.0f, %1.0f] & [%1.0f, %1.0f] \\\\\\\\ \n",
            LimitsInfo[i][0],
            LimitsInfo[i][0+1],
            LimitsInfo[i][3+1],
            LimitsInfo[i][2+1],
            LimitsInfo[i][4+1],
            LimitsInfo[i][1+1],
            LimitsInfo[i][5+1]
          );   
}

printf("\n\n%%[100%% Tgamma]\n");
for(int i=0;i<N_1p0;i++){
    printf("%1.0f & %1.4f & %1.4f & [%1.4f, %1.4f] & [%1.4f, %1.4f] \\\\\\\\ \n",
            Xsec_x_1p0[i],
            Observed_1p0[i],
            Median_1p0[i],
            Minus1sigma_1p0[i],
            Plus1sigma_1p0[i],
            Minus2sigma_1p0[i],
            Plus2sigma_1p0[i]
          );
}
EOF


echo "}" >> templateX.cc

