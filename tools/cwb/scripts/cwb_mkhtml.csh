#!/bin/tcsh -f

onintr irq_ctrlc

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_mkhtml
  exit
endif

if ( "`echo $1 | cut -d'.' -f2`" == "texi" )  then
  # convert texi file

  if (($2 != '') && ($2 != 'pdf') && ($2 != 'header') && ($2 != 'nheader') && ($2 != 'wheader')) then
    echo ""
    echo "cwb_mkhtml" : \'$2\' " bad option"
    echo ""
    exit
  endif

  set fPath = $1 
  set fDir = `echo $fPath| sed 's/.texi//g' `

  echo    '\input texinfo.tex    @c -*-texinfo-*-'         >! $fDir.tmp
  echo    '@setfilename '$fDir.info                        >> $fDir.tmp
  echo    '@afourpaper '                                   >> $fDir.tmp
  echo    '@syncodeindex fn cp '                           >> $fDir.tmp
  echo    '@syncodeindex ky cp '                           >> $fDir.tmp
  echo    ' '                                              >> $fDir.tmp
  echo    '@ifhtml '                                       >> $fDir.tmp
  echo    '@html '                                         >> $fDir.tmp
  echo    '<script type="text/javascript" src="cwb.js"></script>'  >> $fDir.tmp
  echo    '@end html '                                     >> $fDir.tmp
  echo    '@end ifhtml '                                   >> $fDir.tmp
  echo    ' '                                              >> $fDir.tmp
  cat     $fPath                                           >> $fDir.tmp
  echo    ' '                                              >> $fDir.tmp
  echo    '@bye '                                          >> $fDir.tmp

  # create pdf file
  if ($2 == 'pdf') then
    texi2pdf -I $HOME_WAT/tools/cwb/macros/texi_includes -o $fDir.pdf $fDir.tmp
    rm $fDir.tmp
    rm $fDir.aux $fDir.cp  $fDir.fn $fDir.ky 
    rm $fDir.pg  $fDir.tp  $fDir.vr $fDir.log
    exit
  endif

  makeinfo --html --paragraph-indent=0 -I $HOME_WAT/tools/cwb/macros/texi_includes -o $fDir $fDir.tmp
  rm $fDir.tmp
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml.C error : process terminated"
    echo ""
    exit
  endif

  # strip MathJax Tag from files
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${HOME_CWB}/macros/StripMathJax.C\(\"$fDir\"\)
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml.C error : process terminated"
    echo ""
    exit
  endif

  # add cwb header banner
  if (($2 == 'header') || ($2 == 'nheader') || ($2 == 'wheader')) then
    if (($2 == 'header') || ($2 == 'wheader')) then
      # wide header 
      root -n -l -b ${CWB_ROOTLOGON_FILE} ${HOME_CWB}/macros/AddHtmlHeaderFooter.C\(\"$fDir\",false\)
    else 
      # narrow header 
      root -n -l -b ${CWB_ROOTLOGON_FILE} ${HOME_CWB}/macros/AddHtmlHeaderFooter.C\(\"$fDir\",true\)
    endif
  endif

  cp ${HOME_WAT}//html/etc/html/ROOT.css $fDir
  cp ${HOME_WAT}//html/etc/html/ROOT.js  $fDir
  cp ${HOME_WAT}//html/etc/html/tabber.css $fDir
  cp ${HOME_WAT}//html/etc/html/tabber.js  $fDir

else
  # convert file to html
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${HOME_CWB}/macros/cwb_mkhtml_file.C\(\"$1\",\""$2"\"\)
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_file.C error : process terminated"
    echo ""
    exit 1
  endif
  exit 0
endif

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

