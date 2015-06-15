#!/bin/csh

## 1). Postierior : bad with negative signal strength
##      - require no negative        
## 2). Bad fit on muon channel : due to bad PE (the most of bins almost close to zero)
##      - solution : the first bin should be larger than 10

#set filelist=`ls LOG| grep com_limit`
set filelist=`ls LOG| grep com_limit|grep "1000_0\.1"`
# LOG/com_limit_950_1.0_12137.log
# nuisance/1300_0.4_Combined_12113.root

if ( ! -e LOG_Trim ) then
mkdir LOG_Trim
endif

foreach fl(${filelist})

    set mass=`echo ${fl} | awk -F"_" '{print $3}'`
    set br=`echo ${fl} | awk -F"_" '{print $4}'`
    set rnd=`echo ${fl} | awk -F"_" '{print $5}' | awk -F"." '{print $1}'`
    set Npe=`cat LOG/${fl}  | grep "expected bound"|wc| awk '{print $1}'`

    if ( -e LOG_Trim/com_limit_${mass}_${br}_${rnd}.log ) then
        echo "Skip LOG_Trim/com_limit_${mass}_${br}_${rnd}.log"
        continue
    endif

    if ( "${br}" == "0.0" ) then
        set br=0
    else if ( "${br}" == "1.0" ) then
        set br=1
    endif


    echo "Loading nuisance/${mass}_${br}_Combined_${rnd}.root"

    if ( -e nuisance/${mass}_${br}_Combined_${rnd}.root ) then
        ## start to trim
        ## ./trim_posterior nuisance/1000_1_Combined_12123.root 1000 1 12123 5
        ./trim_posterior nuisance/${mass}_${br}_Combined_${rnd}.root ${mass} ${br} ${rnd} ${Npe}
        sleep 0.01
    else
        echo "[WARNING] No such file : nuisance/${mass}_${br}_Combined_${rnd}.root"
    endif
end
