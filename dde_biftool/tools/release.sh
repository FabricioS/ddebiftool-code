#!/bin/bash
#
# $Id$
#
set -e 
#
# Check arguments
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 tagfolder releasefolder version [testprog]"
    echo " * calls svn export putting export of trunk into exportdir=tagfolder/dde_biftool_r{revision}"
    echo " * copies export to maindir=releasefolder/dde_biftool_v{version}"
    echo " * compiles manual, cover and other docs and copies them into maindir,"
    echo " * replaces Id lines and updates (c) lines with {version}(commit)"
    echo " * insert {version} and {revision} into Readme.html and creates Readme.txt from Readme.html"
    echo " * if testprog is given applies testprog to all demos (in temporary folder test)"
    echo "    possible choices for testprog: {matlab octave}"
    echo " * removes folders not intended for distribution"
    echo " * zips maindir into dde_biftool_v{version}.zip"
    echo " "
    echo " used programs:"
    echo " bash, svn, pdflatex, bibtex, python (tested with 2.6), unix2dos, html2text"
    echo " for testing {matlab, octave}"
    exit
fi
curdir=`pwd`
tag=`cd $1;pwd`
cd $curdir
base=`cd $2;pwd`
version=$3
# Redirect stdout ( > ) into a named pipe ( >() ) running "tee"
exec > >(tee "$base/logfile_v$version.txt")
if [[ $# -lt 4 ]]; then
    dotest=0
else
    dotest=1
    cmd=$4
fi
# obtain current revision number
cd $tag
pwd
revision=`svn info | awk -- '/Revision/{print $2}'`
# generate names for folders and files
name="dde_biftool_v"$version
exportdir=$tag"/"$name"_r"$revision
destdir="$name"
destdir=$base/$destdir
license=$destdir/tools/license.txt
zip=$name".zip"
files="*/*.*"
rm -rf $exportdir
svn export ^/trunk/dde_biftool  $exportdir --native-eol CRLF
rm -rf $destdir
cp -urp $exportdir $destdir
#
# set year in license.txt
cd $destdir
year=`date +%G`
sed -i -- "s/|year|/$year/g" $license
#
# compile manuals
cd $destdir/manual
cp -p $license ./license.txt
tex="manual Changes-v3 Extra_psol_extension"
echo '\newcommand{\version}{'$version'}' >version.tex

for x in $tex; do
    pdflatex $x && pdflatex $x && pdflatex $x && \
    bibtex $x && pdflatex $x && pdflatex $x && mv $x".pdf" $destdir
done
mv Addendum_Manual_DDE-BIFTOOL_2_03.pdf $destdir
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
# remove folders not intended for distibution (except tools)
rm -rf  manual FilesChangedAndAdded_V203 system test
#
# create index.html in folders which don't have them
i=0
for i in `find "$destdir" -mindepth 1 -type d`; do
    index=$i"/index.html"
    if [[ ! ( -f $index ) ]]; then
	echo creating $index
	sh tools/create_index.html tools/template_index.html $version $index
    fi
done
rm -rf tools

# zip
cd $base
rm -rf $zip
zip -r $zip $name


