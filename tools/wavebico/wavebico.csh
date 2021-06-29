if ( $1 == '' ) then
  echo "wavebico 'config.C'"
  exit
endif

if ( $1 != '' ) then
  setenv WB_CONFIG_FILE $1
else
  unsetenv WB_CONFIG_FILE 
endif

root -n -l -b ${CWB_ROOTLOGON_FILE} ${HOME_WAVEBICO}/MakeBicoherence.C\(\"${WB_CONFIG_FILE}\"\)

unsetenv WB_CONFIG_FILE
