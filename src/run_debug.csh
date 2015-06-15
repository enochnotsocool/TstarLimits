#!/bin/csh


set WORKSPACE=`pwd`
#set BRs=`seq 0.0 0.1 1.0 | awk '{printf("%1.1f\n", $1)}'`
set BRs=0.0
#set masses=`seq 450 50 950`
#set masses=`cat .list.log`
set masses=1000

set RandomSeeds=`seq 1 1 200`

if ( ! -e ${WORKSPACE}/scripts ) then
mkdir ${WORKSPACE}/scripts
endif

if ( -e .runlist ) then
rm .runlist
endif
touch .runlist


foreach rnd(${RandomSeeds})

foreach br(${BRs})

    foreach m(${masses})

        if ( -e ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh ) then
            rm ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh
        endif 
        touch ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

        echo "#\!/bin/csh \
            cd ${WORKSPACE}/../CMSSW_4_2_8/src; cmsenv;cd - \
            cd ${WORKSPACE}\
            ./stats ${m} ${br} Combined 1000 ${rnd} ">> ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh
        chmod +x ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

#echo "bsub -q 8nh  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh"
#bsub -q 8nh  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

        set lxplusNumber=999
        set IsArmed=0
        while ( ${IsArmed} == "0" ) 
        set lxplusNumber=`date +%s | awk '{print 401+$1%40}'`
        set IsArmed=`ssh lxplus${lxplusNumber} "echo test" | wc|awk '{print $1}'`
        echo "[In While Loop] ${m} ${br} lxplus${lxplusNumber} is armed (${IsArmed})"
        end 

        echo "[Ready for Submit] ssh lxplus${lxplusNumber}  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh "
        ssh lxplus${lxplusNumber}  "${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh" &
        echo "ssh lxplus${lxplusNumber}  "\""${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh"\""" >> .runlist


    end ## masses loop

end ## BRs loop

end ## RandomSeeds loop

