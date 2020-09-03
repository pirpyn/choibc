#!/bin/bash

dir=$(dirname $(readlink -m $0))

# Go to https://www.alglib.net/download.php and get url
last_url="https://www.alglib.net/translator/re/alglib-3.16.0.cpp.gpl.tgz"
url="$(curl -s https://www.alglib.net/download.php | awk '/<a.*\.cpp\.gpl\.tgz/ { print gensub(".*href=\"(.*)\".*","https://www.alglib.net\\1","g",$0);}')"
if [[ $last_url == $url ]]; then
    echo "No update available"
    exit
else
    sed -i "s|$last_url|$url|" $(readlink -m $0)
fi

cd $dir

wget $url
tar xvf alglib-*.tgz
files=(optimization ap alglibinternal linalg alglibmisc solvers stdafx)
for file in ${files[@]}; do
    cp cpp/src/$file.* .
    sed -i 's/^M$//' $file.*
done
rm -f $dir/*.tgz $dir/*.tgz.*
rm -rf $dir/cpp