#!/bin/csh

if ( -e .log_BL1 ) then 
rm .log_BL1
endif

if ( -e .log_BL2 ) then 
rm .log_BL2
endif

if ( -e .black_list ) then
rm .black_list
endif

cat .runlist | grep "In While Loop" | awk '{print $6}' | sort -u>& .log_BL1
cat .runlist | grep "ssh" | awk '{print $2}' | sort -u>& .log_BL2
grep -v -f .log_BL2 .log_BL1 >& .black_list
