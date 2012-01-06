#! /bin/sh

DIRNAME=$(dirname $0)

DESTDIR=${DIRNAME}/../..

if test -e ${DESTDIR}/life.kdevelop ; then
  echo "life.kdevelop already exists. The file will not be copied"
else
  cp ${DIRNAME}/life.kdevelop ${DESTDIR}/life.kdevelop
fi
