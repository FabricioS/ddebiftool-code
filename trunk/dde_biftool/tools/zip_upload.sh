#!/bin/bash
#
# $Id$
#
set -e 
#
# Check arguments
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 username zipfile"
    echo " put (using cp) zipfile to sourceforge filespace for download"
    exit
fi
username=$1
zipfile=$2
scp  $zipfile $username",ddebiftool@web.sourceforge.net:/home/frs/project/ddebiftool/"
