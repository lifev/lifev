#!/bin/sh

logfile="ltxcompile.log"
# redirect errors to a file in user's home directory if we can
if cp /dev/null "$logfile" 2> /dev/null ; then
    chmod a+r "$logfile"
#    exec > "$logfile" 2>&1
fi

showlog(){
  # print 5 lines of leading context before each error
  # pint the line of ltxcompile.log
  cat ltxcompile.log | grep -B 5 -n Error
}
trap "showlog" 1 2 3 4 5 6 7 8 9 10 11

export TEXINPUTS=@CMAKE_CURRENT_SOURCE_DIR@:@CMAKE_SOURCE_DIR@/doc/tex//:$TEXINPUTS
export BIBINPUTS=@CMAKE_CURRENT_SOURCE_DIR@:@CMAKE_SOURCE_DIR@/doc/tex//:$BIBINPUTS

echo export TEXINPUTS=$TEXINPUTS | tee $logfile
echo export BIBINPUTS=$BIBINPUTS | tee -a $logfile

prefix=`basename $1 .tex`

#ltxcmdlines="--interaction nonstopmode --file-line-error-style"
ltxcmdlines="--interaction nonstopmode"

echo prefix=$prefix | tee -a $logfile

pdflatex $ltxcmdlines $prefix | tee -a $logfile
bibtex $prefix                | tee -a $logfile
makeindex $prefix             | tee -a $logfile
pdflatex $ltxcmdlines $prefix | tee -a $logfile
pdflatex $ltxcmdlines $prefix | tee -a $logfile

rm -f $prefix.aux $prefix.bbl $prefix.blg $prefix.log $prefix.out
