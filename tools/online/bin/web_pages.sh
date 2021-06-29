#!/bin/bash 

dir=$1
opt=$2

#source $HOME/.bash_profile

cd $dir

while true; do
    echo "================"
    date
    ./web_pages.py $opt
done

