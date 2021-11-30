#!/bin/bash
# set -vx
# Simple script to switch to using the reduced precision input files
# Finds all 4.2_RC strings in input*.cfg files and appends _LITE.
# Repeat with the -r flag to reverse this process.
#
# This needs to be run in the SETTE directory. 
#
#########################################################################################
######################### Start of function definitions #################################
##
# 
reverse=0
  if [ $# -gt 0 ]; then
    while getopts r option; do 
       case $option in
          r) reverse=1;;
          h | *) echo ''
                 echo 'sette_use_LITE.sh : ' 
                 echo '     Switch to using the reduced precision, LITE input files'
                 echo ' [-r] :'
                 echo '     Switch back to using the normal 4.2_RC files'
                 exit 42;;
       esac
    done
    shift $((OPTIND - 1))
  fi
#
  if [  $reverse == 0 ] ; then
    for file in $( grep -l -e '4.2_RC\.' -e '4.2_RC$' in*cfg )
    do
      perl -0777 -pi -e 's@4.2_RC@4.2_RC_LITE@g'  $file
    done
  else
    for file in $( grep -l -e '4.2_RC_LITE\.' -e '4.2_RC_LITE$' in*cfg )
    do
      perl -0777 -pi -e 's@4.2_RC_LITE@4.2_RC@g'  $file
    done
  fi
#
exit
