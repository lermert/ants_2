#!/bin/bash

for file in $1/* ;

do
	
	filesize=$(du $file | awk '{print $1}');
	
	if [ $filesize == 0 ] ; then
	
		
		rm $file;
		
	fi;
done
