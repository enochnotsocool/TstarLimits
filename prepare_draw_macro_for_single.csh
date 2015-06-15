#!/bin/csh

#Observed Limit: r < 1.3737
#Expected  2.5%: r < 0.1674
#Expected 16.0%: r < 0.2273
#Expected 50.0%: r < 0.3271
#Expected 84.0%: r < 0.4849
#Expected 97.5%: r < 0.7020

set BR=0.0
#set BR=0.5

if ( $1 != "" ) then
set BR="$1"
endif


set limitValue=`printf "Observed 2.5%% 16.0%% 50.0%% 84.0%% 97.5%%"`
#set  mass_list=`ls LIMIT_LOG/ | grep _limit_Fit | grep "_${BR}\.txt\.log" | awk -F"datacard_Tgamma" '{print $2" "$1"datacard_Tgamma"$2}' | sort -n|awk '{print $2}' |  grep -v 550 | grep -v 600 | grep -v 1500`
set  mass_list=`ls LIMIT_LOG/ | grep _limit_Fit | grep "_${BR}\.txt\.log" | awk -F"datacard_Tgamma" '{print $2" "$1"datacard_Tgamma"$2}' | sort -n|awk '{print $2}' |  grep -v 550 |  grep -v 1500`

## .show_limit_combinedFit_part1.cc
if ( -e show_limit_combinedFit.cc ) then
rm show_limit_combinedFit.cc
endif
touch show_limit_combinedFit.cc
echo "float BranchingRatio = ${BR};" >> show_limit_combinedFit.cc

@ counter=1
foreach lv(${limitValue})

    cat .show_limit_combinedFit_part${counter}.cc >> show_limit_combinedFit.cc
    echo "// ${lv}" >> show_limit_combinedFit.cc
    foreach ml(${mass_list})
    grep ${lv} LIMIT_LOG/${ml} | awk '{print $5", //'${ml}'"}'  >> show_limit_combinedFit.cc

    end

    echo ""

@ counter++
end


cat .show_limit_combinedFit_part7.cc >> show_limit_combinedFit.cc

root -l -b -q show_limit_combinedFit.cc
set br_=`echo ${BR} | sed 's/\./p/g'`

if ( -e PLOTS/show_limit_combinedFit_${br_}.cc ) then
rm PLOTS/show_limit_combinedFit_${br_}.cc
endif
sed "s/show_limit_combinedFit/show_limit_combinedFit_${br_}/g" show_limit_combinedFit.cc >& PLOTS/show_limit_combinedFit_${br_}.cc
mv LimitCombine.pdf PLOTS/LimitCombine_${br_}.pdf

echo "Please check plots and macro in PLOTS folder"
