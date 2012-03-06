#! /bin/bash

cd life
for DIR in life*; do
    echo $DIR
    cd $DIR
    mv -f CMakeLists.txt CMakeLists.txt.old
    CPPFILES=$(find . -name \*.cpp | sed 's,./,,')
    HPPFILES=$(ls *.hpp)
    cat > CMakeLists.txt << eof
add_library ( $DIR $CPPFILES )
install ( TARGETS $DIR DESTINATION lib )
install ( FILES $HPPFILES DESTINATION include/life/$DIR)
eof
if [ $DIR == lifefunctions ]; then
cat >> CMakeLists.txt << eof
install ( FILES bessel/bessel.hpp DESTINATION include/life/$DIR/bessel)
eof
fi
cd ..
done
cd ..
