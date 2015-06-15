#!/bin/csh


set WORKSPACE=`pwd`
#set BRs=`seq 0.0 0.1 1.0 | awk '{printf("%1.1f\n", $1)}'`
set BRs=0.0
#set masses=`seq 450 50 950`
set masses=`cat .list.log`

set RandomSeeds=`seq 1 1 10`

if ( ! -e ${WORKSPACE}/scripts ) then
mkdir ${WORKSPACE}/scripts
endif

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
            ./stats ${m} ${br} Combined 10000 ${rnd} ">> ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh
        chmod +x ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

        echo "bsub -q 8nh  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh"
        bsub -q 8nh  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

    end ## masses loop

end ## BRs loop

end ## RandomSeeds loop

