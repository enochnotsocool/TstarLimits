#!/bin/csh

## FOR SLC6

set WORKSPACE=`pwd`
set BRs=`seq 0.0 0.1 1.0 | awk '{printf("%1.1f\n", $1)}'`
#set BRs=`echo "0.0 0.3 0.5 0.7 0.8 0.9 1.0"`
#set BRs=`echo "0.1 0.2 0.4 0.6"`
#set BRs=`seq 0.0 1 1.0 | awk '{printf("%1.1f\n", $1)}'`
#set BRs=`seq 0.0 1 1.0 | awk '{printf("%1.1f\n", $1)}'`


#set BRs="0.0 1.0"
#set BRs="0.2 0.8"
#set BRs="0.4 0.5 0.6"
#set BRs="0.1 0.9"
#set BRs="0.3 0.7"

set masses=`cat .list.log`

set RandomSeeds=`seq 18281 1 18300`

if ( ! -e thisroot.csh ) then
    set source_Path="/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/root/5.27.06b-cms23/bin/thisroot.csh"
    set s_Path=`echo $source_Path | sed 's/\//\\\//g'`
    cat $source_Path  | sed 's/set ARGS=/set ARGS=\`echo "source '"$s_Path"'"\`/g' | sed 's/($_)//g' >& thisroot.csh
endif

if ( ! -e ${WORKSPACE}/scripts ) then
mkdir ${WORKSPACE}/scripts
endif

if ( -e .runlist ) then
rm .runlist
endif
touch .runlist

echo "Running on `hostname`" >> .runlist
echo "" >> .runlist
echo "" >> .runlist

hostname >> .black_list

set counter=0

foreach rnd(${RandomSeeds})

foreach br(${BRs})

    foreach m(${masses})

        if ( -e ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh ) then
            rm ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh
        endif 
        touch ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

        echo "#\!/bin/csh \
            source ${WORKSPACE}/thisroot.csh \
            cd ${WORKSPACE}\
            ./stats ${m} ${br} Combined 10000 ${rnd} ">> ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh
        chmod +x ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

#echo "bsub -q 8nh  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh"
#bsub -q 8nh  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh

        set lxplusNumber=999
        set IsArmed=0
        while ( ${IsArmed} == "0" ) 
            #set lxplusNumber=`date +%s | awk '{print 401+$1%40}'`
            set lxplusNumber=`date +%s | awk '{printf("%04i", 001+$1%245)}'`
            set IsInBlackList=`grep ${lxplusNumber} .black_list| wc|awk '{print $1}'`

            if ( ${IsInBlackList} == "0" ) then
                echo "       [In While Loop] ${m} ${br} lxplus${lxplusNumber} is armed (??)" >> .runlist
                set IsArmed=`ssh lxplus${lxplusNumber} "echo test" | wc|awk '{print $1}'`
                echo "[In While Loop] ${m} ${br} lxplus${lxplusNumber} is armed (${IsArmed})"
            endif
        end 

        echo "[Ready for Submit] ssh lxplus${lxplusNumber}  ${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh "
        ssh lxplus${lxplusNumber}  "${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh" &
        echo "ssh lxplus${lxplusNumber}  "\""${WORKSPACE}/scripts/job_M${m}_BR${br}_rnd${rnd}.csh"\""" >> .runlist

        set counter=`echo ${counter} | awk '{print $1+1}'`

    end ## masses loop

end ## BRs loop

#set IsSleep=`echo ${rnd} | awk '{if($1%3==2){printf("1");}else{printf("0");}}'`
#set IsSleep=`echo ${rnd} | awk '{if($1%6==5){printf("1");}else{printf("0");}}'`
#set IsSleep=`echo ${counter} | awk '{if($1%200==199){printf("1");}else{printf("0");}}'`    ## SLC5
set IsSleep=`echo ${counter} | awk '{if($1%750==749){printf("1");}else{printf("0");}}'`     ## SLC6

if ( "${IsSleep}" == "1" ) then
#sleep 1800
sleep 3600
#sleep 5400
endif

end ## RandomSeeds loop




set HN=`hostname`
grep -v ${HN} .black_list >& .black_list_tmp
mv .black_list_tmp .black_list
