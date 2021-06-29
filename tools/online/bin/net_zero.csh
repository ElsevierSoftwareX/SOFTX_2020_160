#!/bin/tcsh -f

#setenv CWB_UPARAMETERS_FILE $1
ln -s $1 config/user_parameters.C
${HOME_CWB}/scripts/cwb_inet.csh 1 >& log &
#rm config/user_parameters.C
