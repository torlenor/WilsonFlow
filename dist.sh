#!/bin/bash
MAJOR_VERSION=0
MINOR_VERSION=0
REVISION_VERSION=2

NAME="wilsonflow"
FILES="flow.cpp init.cpp options.cpp matrixhandling.cpp su3.cc CArray.hpp init.h measuretime.h su3.h ConfigData.hpp matrixhandling.h options.h CMakeLists.txt"

export DISTDIR="${NAME}_${MAJOR_VERSION}.${MINOR_VERSION}.${REVISION_VERSION}"

if [ -d $DISTDIR ]
then
	rm -rf $DISTDIR &>/dev/null
fi
mkdir $DISTDIR
echo -e "#define MAJOR_VERSION $MAJOR_VERSION\n#define MINOR_VERSION $MINOR_VERSION\n#define REVISION_VERSION $REVISION_VERSION" > ${DISTDIR}/version.h
for f in $FILES ;do
  cp $f $DISTDIR
done

tar cfz ${DISTDIR}.tgz $DISTDIR

rm -rf $DISTDIR

echo "$DISTDIR.tgz ready to ship!"
