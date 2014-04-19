#!/bin/bash
#
# $Id: release.sh 39 2013-06-13 12:31:32Z jansie $
#
set -e 
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 basefolder version [testprog]"
    echo " * calls svn export putting export into maindir=basefolder/dde_biftool_v{version}"
    echo " * compiles manual, cover and other docs and copies them into maindir,"
    echo " * replaces Id lines and updates (c) lines with {version}(commit)"
    echo " * insert {version} into Readme.html and creates Readme.txt from Readme.html"
    echo " * if testprog is given applies testprog to all demos (in temporary folder test)"
    echo "    possible choices for testprog: matlab octave"
    echo " * removes folders not intended for distribution"
    echo " * zips maindir into dde_biftool_v{version}.zip"
    echo " "
    echo " used programs:"
    echo " bash, svn, pdflatex, bibtex, python (tested with 2.6), unix2dos, html2text"
    echo " for testing {matlab, octave}"
    exit
fi
base=`cd $1;pwd`
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
license=$destdir/tools/license.txt
zip=$name".zip"
files="*/*.*"
rm -rf $destdir
svn export ^/trunk/dde_biftool  $destdir --native-eol CRLF
#
# set year in license.txt
cd $destdir
year=`date +%G`
sed -i -- "s/|year|/$year/g" $license
#
# compile manual
cd $destdir/manual
cp -p $license ./license.txt
tex="manual ExtensionChanges-v3"
echo '\newcommand{\version}{'$version'}' >version.tex

for x in $tex; do
    pdflatex $x && pdflatex $x && pdflatex $x && \
    bibtex $x && pdflatex $x && pdflatex $x
done
mv manual.pdf ExtensionChanges-v3.pdf Addendum_Manual_DDE-BIFTOOL_2_03.pdf $destdir
cd $destdir
#
# set (c) line in all m and html files that have (c) or Id
python $destdir/tools/c_insert.py $destdir $version
#
# insert license and set version in readme
nr=`awk -- '/\|license\|/{print NR}' Readme.html`
awk -- "NR<$nr"'{print $0}' Readme.html >tmp.txt
awk -- 'BEGIN{RS="\r\n"}{if(NF==0)print "<br><br>";else print $0}' $license >>tmp.txt
awk -- "NR>$nr"'{print $0}' Readme.html >>tmp.txt
sed  -- "s/|version|/$version/g" tmp.txt > Readme.html
rm -f tmp.txt
html2text Readme.html >Readme.txt
unix2dos Readme.txt
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


