#!/bin/bash
#To extract a release, do the following:
export SUBREL=5.0.4
rm -rf ../ramd-${SUBREL}*
git archive --format=tar --prefix=ramd-${SUBREL}/ master | gzip -9 > ../ramd-${SUBREL}.tgz
cd .. &&  rm -rf tmp && mkdir tmp && cd tmp
tar xf ../ramd-${SUBREL}.tgz
rm -r ramd-${SUBREL}/TRJ-Analysis-R
zip -rm ../ramd-only-${SUBREL}.zip ramd-${SUBREL}
cd .. && rmdir tmp
