#!/bin/bash


#Set service base path, change to your own service host
SERVICEBASE="http://service.scedc.caltech.edu"

# Set all service specific locations using the service base
TIMESERIESWS="${SERVICEBASE}/fdsnws/dataselect/1" 
METADATAWS="${SERVICEBASE}/fdsnws/station/1" 
EVENTWS="${SERVICEBASE}/fdsnws/event/1" 
SACPZWS="${SERVICEBASE}/scedcws/sacpz/1" 
RESPWS="${SERVICEBASE}/scedcws/resp/1"

export SERVICEBASE TIMESERIESWS METADATAWS EVENTWS SACPZWS RESPWS

/home/lermert/code/ants_2/ants_2/tools_ext/FetchData "$@"

