#!/bin/bash

mv log/${argv[JOB]}*  ../log/ ;
cd output ; mv *${argv[JOB]}*  ../../output/ ;

#exit 0
