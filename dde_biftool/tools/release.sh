#!/bin/bash
#
# $Id: release.sh 39 2013-06-13 12:31:32Z jansie $
#
set -e 
version=$1
name="ddebiftool_v"$version
destdir="tags/$name"
zip=$name".zip"
files="*/*.*"
rm -rf $destdir
svn export ^/trunk/dde_biftool  $destdir --native-eol CRLF
base=`pwd`
destdir=$base/$destdir
#
# compile manual
cd $destdir/manual
tex="TWcover manual ExtensionChanges-v3"
echo '\newcommand{\version}{'$version'}' >version.tex
for x in $tex; do
    pdflatex $x && pdflatex $x && pdflatex $x && \
    bibtex $x && pdflatex $x && pdflatex $x
done
pdftk A=TWcover.pdf B=manual.pdf cat A B output all.pdf
mv all.pdf $destdir/manual.pdf
mv ExtensionChanges-v3.pdf Addendum_Manual_DDE-BIFTOOL_2_03.pdf $destdir
# clean up manual creation
cd $destdir
#
#
python $destdir/tools/c_insert.py $destdir $version

# remove folders not intended for distibution
rm -rf  tools manual FilesChangedAndAdded_V203 system

#correct author name
#find . -type f -print0 | xargs -0 sed -i '/$Id/s/jansie/Jan Sieber/g'
#find . -type f -print0 | xargs -0 sed -i '/$Id/s/js543/Jan Sieber/g'

#insert release number
#find . -type f -print0 | xargs -0 sed -i '/$Id/s/$Id/$Id('"$version"')/g'

# zip
cd $base/tags
rm -rf $zip
zip -r $zip $name
# copy zip file to Dropbox
#cp -p $zip $HOME/Dropbox/Public/ddebiftool_ext.zip

