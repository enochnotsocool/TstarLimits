#!/bin/csh

#set Need2Sleep=`seq 0 1 5400`
#foreach nsleep(${Need2Sleep})
#    echo ${nsleep} | awk '{print "Need "5400-$1" seconds to submit!"}'
#    sleep 1
#end

## FOR SLC6

set WORKSPACE=`pwd`
#set BRs=`seq 0.0 0.1 1.0 | awk '{printf("%1.1f\n", $1)}'`
#set BRs=`echo "0.0 0.3 0.5 0.7 0.8 0.9 1.0"`
#set BRs=`echo "0.1 0.2 0.4 0.6"`
#set BRs=`seq 0.0 1 1.0 | awk '{printf("%1.1f\n", $1)}'`
#set BRs=`seq 0.0 1 1.0 | awk '{printf("%1.1f\n", $1)}'`


#set BRs="0.0 1.0"
#set BRs="0.2 0.8"
#set BRs="0.4 0.5 0.6"
#set BRs="0.1 0.9"
#set BRs="0.3 0.7"

set BRs="0.0 1.0 0.2 0.8 0.4 0.5 0.6 0.1 0.9 0.3 0.7"


set masses=`cat .list.log`

#set RandomSeeds=`seq 58281 1 58300`
#set RandomSeeds=`seq 58781 1 58800`
#set RandomSeeds=`seq 59581 1 59600`
#set RandomSeeds=`seq 69581 1 69590`
set RandomSeeds=`seq 69681 1 69690`

if ( ! -e thisroot.csh ) then
    set source_Path="/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/root/5.27.06b-cms23/bin/thisroot.csh"
    set s_Path=`echo $source_Path | sed 's/\//\\\//g'`
    cat $source_Path  | sed 's/set ARGS=/set ARGS=\`echo "source '"$s_Path"'"\`/g' | sed 's/($_)//g' >& thisroot.csh
endif

if ( ! -e ${WORKSPACE}/scripts ) then
mkdir ${WORKSPACE}/scripts
endif

if ( -e .runlist_slc5 ) then
rm .runlist_slc5
endif
touch .runlist_slc5

echo "Running on `hostname`" >> .runlist_slc5
echo "" >> .runlist_slc5
echo "" >> .runlist_slc5


set counter=0

foreach rnd(${RandomSeeds})

foreach br(${BRs})

    foreach m(${masses})

        if ( -e ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh ) then
            rm ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh
        endif 
        touch ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

        echo "#\!/bin/csh \
            cd ${WORKSPACE}/../CMSSW_4_2_8/src/ ; cmsenv\
            cd ${WORKSPACE}\
            ./stats ${m} ${br} Combined 10000 ${rnd} ">> ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh
        chmod +x ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

#echo "bsub -q 8nh  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh"
#bsub -q 8nh  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

        set lxplusNumber=999
        set IsArmed=0
        while ( ${IsArmed} == "0" ) 
            set lxplusNumber=`date +%s | awk '{print 401+$1%48}'`
#set lxplusNumber=`date +%s | awk '{printf("%04i", 001+$1%245)}'`
            echo "       [In While Loop] ${m} ${br} lxplus${lxplusNumber} is armed (??)" >> .runlist_slc5
            set IsArmed=`ssh lxplus${lxplusNumber} "echo test" | wc|awk '{print $1}'`
            echo "[In While Loop] ${m} ${br} lxplus${lxplusNumber} is armed (${IsArmed})"
        end 

        echo "[Ready for Submit] ssh lxplus${lxplusNumber}  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh "
        ssh lxplus${lxplusNumber}  "${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh" &
        echo "ssh lxplus${lxplusNumber}  "\""${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh"\""" >> .runlist_slc5

        set counter=`echo ${counter} | awk '{print $1+1}'`


    end ## masses loop

end ## BRs loop

        set IsSleep=`echo ${counter} | awk '{if($1%200==199){printf("1");}else{printf("0");}}'`    ## SLC5
        if ( "${IsSleep}" == "1" ) then
        sleep 5400
        endif

end ## RandomSeeds loop

