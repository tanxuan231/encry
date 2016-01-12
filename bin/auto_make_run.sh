#!/bin/bash

cur_pwd=`pwd`

cd "../ldecod"
make $2 || exit
cd $cur_pwd

if [[ $1 -eq 1 ]]
then
if [ $2 == "DBG=1" ]
then
	./ldecod.dbg.exe
else
	./ldecod.exe
fi
fi
