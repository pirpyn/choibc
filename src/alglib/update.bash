#!/bin/bash

dir=$(dirname $(readlink -m $0))
cd $dir
# Go to https://www.alglib.net/download.php and get url
url="https://www.alglib.net/translator/re/alglib-3.16.0.cpp.gpl.tgz"
wget $url
tar xvf alglib-*.tgz
files=(optimization ap alglibinternal linalg alglibmisc solvers stdafx)
for file in ${files[@]}; do
    cp cpp/src/$file.* .
done
rm -f $dir/*.tgz $dir/*.tgz.*
rm -rf $dir/cpp