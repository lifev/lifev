#! /bin/sh

DIRNAME=$(dirname $0)

DESTDIR=${DIRNAME}/../..

if test -e ${DESTDIR}/.project ; then
  echo ".project already exists. The file will not be copied"
else
  cp ${DIRNAME}/project ${DESTDIR}/.project
fi

if test -e ${DESTDIR}/.cproject ; then
  echo ".cproject already exists. The file will not be copied"
else
  cp ${DIRNAME}/cproject ${DESTDIR}/.cproject
fi

if test -e ${DESTDIR}/.settings ; then
  echo ".settings already exists. The directory will not be copied"
else
  cp -r ${DIRNAME}/settings ${DESTDIR}/.settings
fi


