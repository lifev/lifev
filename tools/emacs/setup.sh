#! /bin/sh

DIRNAME=$(dirname $0)

DESTDIR=${DIRNAME}/../..

cd ${DESTDIR}
if test -e .emacs-dirvars ; then
  echo ".emacs-dirvars already exists. The file will not be linked"
else    
  ln -s tools/emacs/emacs-dirvars .emacs-dirvars
fi

if test -e  Templates; then
  echo "directory Templates already exists. The directiry will not be linked."
else 
 ln -s tools/emacs/Templates .
fi
