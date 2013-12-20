#!/bin/bash

atom=$1
ion=$2

if [ -z ${atom} ] ; then export atom="1" ; fi
if [ -z ${ion} ] ; then export ion="1" ; fi

let ion=ion-1
echo "atom/ion:" ${atom} ${ion}

if [ ${atom} -lt 10 ]  ; then atom="0"${atom} ; fi
if [ ${ion} -lt 10 ]  ; then ion="0"${ion} ; fi
echo "atom/ion after cvt:" ${atom} ${ion}


wget --recursive --page-requisites --html-extension --convert-links --domains kurucz.harvard.edu --no-parent kurucz.harvard.edu/atoms/${atom}${ion}/ 

