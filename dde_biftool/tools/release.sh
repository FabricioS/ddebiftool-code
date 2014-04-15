#!/bin/bash
#
# $Id: release.sh 39 2013-06-13 12:31:32Z jansie $
#
set -e 
base=$1 
version=$2
if [[ $# -lt 3 ]]; then
    dotest=0
else
    dotest=1
    cmd=$3
fi
name="dde_biftool_v"$version
destdir="$name"
destdir=$base/$destdir
zip=$name".zip"
files="*/*.*"
rm -rf $destdir
svn export ^/trunk/dde_biftool  $destdir --native-eol CRLF
#
# compile manual
cd $destdir/manual
tex="TWcover manual ExtensionChanges-v3"
echo '\newcommand{\version}{'$version'}' >version.tex
for x in $tex; do
    pdflatex $x && pdflatex $x && pdflatex $x && \
    bibtex $x && pdflatex $x && pdflatex $x
done
#
# join cover and manual body
pdfjoin=`which pdftk`
if [[ -x $pdfjoin ]]; then
    $pdfjoin A=TWcover.pdf B=manual.pdf cat A B output all.pdf
else
    gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER \
    -sOutputFile=all.pdf TWcover.pdf manual.pdf
fi
mv all.pdf $destdir/manual.pdf
mv ExtensionChanges-v3.pdf Addendum_Manual_DDE-BIFTOOL_2_03.pdf $destdir
# clean up manual creation
cd $destdir
#
# set (c) line in all m and html files that have (c) or Id
python $destdir/tools/c_insert.py $destdir $version
#
# set version readme
sed -i -- "s/|version|/$version/g" Readme.html
html2text Readme.html >Readme.txt
#
# if testing required perform tests of demos
if [[ $dotest -eq 1 ]]; then
    mkdir $destdir/test
    python $destdir/tools/test_demos.py $destdir $cmd
fi
#
# remove folders not intended for distibution
rm -rf  tools manual FilesChangedAndAdded_V203 system test


# zip
cd $base
rm -rf $zip
zip -r $zip $name


