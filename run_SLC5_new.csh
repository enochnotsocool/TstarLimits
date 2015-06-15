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


#set BRs="0.0 0.5 1.0"
#set BRs="0.2 0.8"
#set BRs="0.4 0.6"
#set BRs="0.1 0.9"
#set BRs="0.3 0.7"

#set BRs="0.0 1.0 0.2 0.8 0.4 0.5 0.6 0.1 0.9 0.3 0.7"

#set BRs="0.2 0.8 0.4 0.6 0.1 0.9 0.3 0.7"
#set BRs="0.5 0.8 0.9"
#set BRs="0.1 0.2"
#set BRs="0.3 0.4"
#set BRs="0.5"

#set masses=`cat .list.log`

#set BRs="0.0"
#set masses="1000 1100"

#set RandomSeeds=`seq 12101 1 12180`
#set RandomSeeds=`seq 12201 1 12280`
#set RandomSeeds=`seq 12301 1 12380`
#set RandomSeeds=`seq 22301 1 22380`
#set RandomSeeds=`seq 22401 1 22480`
#set RandomSeeds=`seq 23401 1 23480`
#set RandomSeeds=`seq 22801 1 22880`
set RandomSeeds="22325"

#set RandomSeeds=`seq 12101 1 12140`
#set RandomSeeds=`seq 12201 1 12240`
set BRs="0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
set masses="600 650 700 750 800 850 900 950 1000 1100 1200 1300 1400"

#set masses="1000 1100 1200"
#set masses="600 650 700 750 800 850 900 950  1300 1400"

#set BRs=" 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 "
#set BRs="0.0"
#set BRs="0.5"
#set BRs="0.1 0.9"
#set BRs="0.2"
#set BRs="0.8"
#set BRs="0.3"
#set BRs="0.7"
#set BRs="0.4"
#set BRs="0.6"
#set masses="600 650 700 750 800 850 900 950 1000 1100 1200 1300 1400"

set BRs="0.1"
set masses="1000 "
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
            cd ${WORKSPACE}/../ ; source env.csh\
            cd ${WORKSPACE}\
            ./stats ${m} ${br} Combined 10000 ${rnd} 1\
            rm -r ${WORKSPACE}/LSF*">> ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh
        chmod +x ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

#        echo "bsub -q 1nd  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh"
#        bsub -q 1nd  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh
        echo "qsub -q cms   ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh"
        qsub -q cms   ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

#        set lxplusNumber=999
#        set IsArmed=0
#        while ( ${IsArmed} == "0" ) 
#            set lxplusNumber=`date +%s | awk '{print 401+$1%48}'`
#            echo "       [In While Loop] ${m} ${br} lxplus${lxplusNumber} is armed (??)" >> .runlist_slc5
#            set IsArmed=`ssh lxplus${lxplusNumber} "echo test" | wc|awk '{print $1}'`
#            echo "[In While Loop] ${m} ${br} lxplus${lxplusNumber} is armed (${IsArmed})"
#        end 
#
#        echo "[Ready for Submit] ssh lxplus${lxplusNumber}  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh "
#        ssh lxplus${lxplusNumber}  "${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh" &
#        echo "ssh lxplus${lxplusNumber}  "\""${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh"\""" >> .runlist_slc5
#
#        set counter=`echo ${counter} | awk '{print $1+1}'`


    end ## masses loop
    sleep 1

end ## BRs loop

#        set IsSleep=`echo ${counter} | awk '{if($1%200==199){printf("1");}else{printf("0");}}'`    ## SLC5
#        if ( "${IsSleep}" == "1" ) then
#        sleep 5400
#        endif

end ## RandomSeeds loop

