#!/bin/bash
g++ -DHAVE_CONFIG_H -I.  -I./cudd -I./mtr -I./epd -I./st   -Wall -Wextra -std=c++0x -g -O3 -MT cplusplus/cplusplus_testobj-testobj.o -MD -MP -MF cplusplus/.deps/cplusplus_testobj-testobj.Tpo -c -o cplusplus/cplusplus_testobj-testobj.o `test -f 'cplusplus/testobj.cc' || echo './'`cplusplus/testobj.cc
mv -f cplusplus/.deps/cplusplus_testobj-testobj.Tpo cplusplus/.deps/cplusplus_testobj-testobj.Po
/bin/sh ./libtool  --tag=CXX   --mode=compile g++ -DHAVE_CONFIG_H -I.  -I./cudd -I./mtr -I./epd -I./st   -Wall -Wextra -std=c++0x -g -O3 -MT cplusplus/cplusplus_libobj_la-cuddObj.lo -MD -MP -MF cplusplus/.deps/cplusplus_libobj_la-cuddObj.Tpo -c -o cplusplus/cplusplus_libobj_la-cuddObj.lo `test -f 'cplusplus/cuddObj.cc' || echo './'`cplusplus/cuddObj.cc
mv -f cplusplus/.deps/cplusplus_libobj_la-cuddObj.Tpo cplusplus/.deps/cplusplus_libobj_la-cuddObj.Plo
/bin/sh ./libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o cplusplus/libobj.la  cplusplus/cplusplus_libobj_la-cuddObj.lo  
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddAddAbs.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddAddAbs.Tpo -c -o cudd/cudd_libcudd_la-cuddAddAbs.lo `test -f 'cudd/cuddAddAbs.c' || echo './'`cudd/cuddAddAbs.c
mv -f cudd/.deps/cudd_libcudd_la-cuddAddAbs.Tpo cudd/.deps/cudd_libcudd_la-cuddAddAbs.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddAddApply.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddAddApply.Tpo -c -o cudd/cudd_libcudd_la-cuddAddApply.lo `test -f 'cudd/cuddAddApply.c' || echo './'`cudd/cuddAddApply.c
mv -f cudd/.deps/cudd_libcudd_la-cuddAddApply.Tpo cudd/.deps/cudd_libcudd_la-cuddAddApply.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddAddFind.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddAddFind.Tpo -c -o cudd/cudd_libcudd_la-cuddAddFind.lo `test -f 'cudd/cuddAddFind.c' || echo './'`cudd/cuddAddFind.c
mv -f cudd/.deps/cudd_libcudd_la-cuddAddFind.Tpo cudd/.deps/cudd_libcudd_la-cuddAddFind.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddAddInv.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddAddInv.Tpo -c -o cudd/cudd_libcudd_la-cuddAddInv.lo `test -f 'cudd/cuddAddInv.c' || echo './'`cudd/cuddAddInv.c
mv -f cudd/.deps/cudd_libcudd_la-cuddAddInv.Tpo cudd/.deps/cudd_libcudd_la-cuddAddInv.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddAddIte.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddAddIte.Tpo -c -o cudd/cudd_libcudd_la-cuddAddIte.lo `test -f 'cudd/cuddAddIte.c' || echo './'`cudd/cuddAddIte.c
mv -f cudd/.deps/cudd_libcudd_la-cuddAddIte.Tpo cudd/.deps/cudd_libcudd_la-cuddAddIte.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddAddNeg.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddAddNeg.Tpo -c -o cudd/cudd_libcudd_la-cuddAddNeg.lo `test -f 'cudd/cuddAddNeg.c' || echo './'`cudd/cuddAddNeg.c
mv -f cudd/.deps/cudd_libcudd_la-cuddAddNeg.Tpo cudd/.deps/cudd_libcudd_la-cuddAddNeg.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddAddWalsh.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddAddWalsh.Tpo -c -o cudd/cudd_libcudd_la-cuddAddWalsh.lo `test -f 'cudd/cuddAddWalsh.c' || echo './'`cudd/cuddAddWalsh.c
mv -f cudd/.deps/cudd_libcudd_la-cuddAddWalsh.Tpo cudd/.deps/cudd_libcudd_la-cuddAddWalsh.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddAndAbs.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddAndAbs.Tpo -c -o cudd/cudd_libcudd_la-cuddAndAbs.lo `test -f 'cudd/cuddAndAbs.c' || echo './'`cudd/cuddAndAbs.c
mv -f cudd/.deps/cudd_libcudd_la-cuddAndAbs.Tpo cudd/.deps/cudd_libcudd_la-cuddAndAbs.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddAnneal.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddAnneal.Tpo -c -o cudd/cudd_libcudd_la-cuddAnneal.lo `test -f 'cudd/cuddAnneal.c' || echo './'`cudd/cuddAnneal.c
mv -f cudd/.deps/cudd_libcudd_la-cuddAnneal.Tpo cudd/.deps/cudd_libcudd_la-cuddAnneal.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddApa.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddApa.Tpo -c -o cudd/cudd_libcudd_la-cuddApa.lo `test -f 'cudd/cuddApa.c' || echo './'`cudd/cuddApa.c
mv -f cudd/.deps/cudd_libcudd_la-cuddApa.Tpo cudd/.deps/cudd_libcudd_la-cuddApa.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddAPI.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddAPI.Tpo -c -o cudd/cudd_libcudd_la-cuddAPI.lo `test -f 'cudd/cuddAPI.c' || echo './'`cudd/cuddAPI.c
mv -f cudd/.deps/cudd_libcudd_la-cuddAPI.Tpo cudd/.deps/cudd_libcudd_la-cuddAPI.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddApprox.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddApprox.Tpo -c -o cudd/cudd_libcudd_la-cuddApprox.lo `test -f 'cudd/cuddApprox.c' || echo './'`cudd/cuddApprox.c
mv -f cudd/.deps/cudd_libcudd_la-cuddApprox.Tpo cudd/.deps/cudd_libcudd_la-cuddApprox.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddBddAbs.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddBddAbs.Tpo -c -o cudd/cudd_libcudd_la-cuddBddAbs.lo `test -f 'cudd/cuddBddAbs.c' || echo './'`cudd/cuddBddAbs.c
mv -f cudd/.deps/cudd_libcudd_la-cuddBddAbs.Tpo cudd/.deps/cudd_libcudd_la-cuddBddAbs.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddBddCorr.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddBddCorr.Tpo -c -o cudd/cudd_libcudd_la-cuddBddCorr.lo `test -f 'cudd/cuddBddCorr.c' || echo './'`cudd/cuddBddCorr.c
mv -f cudd/.deps/cudd_libcudd_la-cuddBddCorr.Tpo cudd/.deps/cudd_libcudd_la-cuddBddCorr.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddBddIte.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddBddIte.Tpo -c -o cudd/cudd_libcudd_la-cuddBddIte.lo `test -f 'cudd/cuddBddIte.c' || echo './'`cudd/cuddBddIte.c
mv -f cudd/.deps/cudd_libcudd_la-cuddBddIte.Tpo cudd/.deps/cudd_libcudd_la-cuddBddIte.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddBridge.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddBridge.Tpo -c -o cudd/cudd_libcudd_la-cuddBridge.lo `test -f 'cudd/cuddBridge.c' || echo './'`cudd/cuddBridge.c
mv -f cudd/.deps/cudd_libcudd_la-cuddBridge.Tpo cudd/.deps/cudd_libcudd_la-cuddBridge.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddCache.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddCache.Tpo -c -o cudd/cudd_libcudd_la-cuddCache.lo `test -f 'cudd/cuddCache.c' || echo './'`cudd/cuddCache.c
mv -f cudd/.deps/cudd_libcudd_la-cuddCache.Tpo cudd/.deps/cudd_libcudd_la-cuddCache.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddCheck.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddCheck.Tpo -c -o cudd/cudd_libcudd_la-cuddCheck.lo `test -f 'cudd/cuddCheck.c' || echo './'`cudd/cuddCheck.c
mv -f cudd/.deps/cudd_libcudd_la-cuddCheck.Tpo cudd/.deps/cudd_libcudd_la-cuddCheck.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddClip.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddClip.Tpo -c -o cudd/cudd_libcudd_la-cuddClip.lo `test -f 'cudd/cuddClip.c' || echo './'`cudd/cuddClip.c
mv -f cudd/.deps/cudd_libcudd_la-cuddClip.Tpo cudd/.deps/cudd_libcudd_la-cuddClip.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddCof.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddCof.Tpo -c -o cudd/cudd_libcudd_la-cuddCof.lo `test -f 'cudd/cuddCof.c' || echo './'`cudd/cuddCof.c
mv -f cudd/.deps/cudd_libcudd_la-cuddCof.Tpo cudd/.deps/cudd_libcudd_la-cuddCof.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddCompose.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddCompose.Tpo -c -o cudd/cudd_libcudd_la-cuddCompose.lo `test -f 'cudd/cuddCompose.c' || echo './'`cudd/cuddCompose.c
mv -f cudd/.deps/cudd_libcudd_la-cuddCompose.Tpo cudd/.deps/cudd_libcudd_la-cuddCompose.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddDecomp.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddDecomp.Tpo -c -o cudd/cudd_libcudd_la-cuddDecomp.lo `test -f 'cudd/cuddDecomp.c' || echo './'`cudd/cuddDecomp.c
mv -f cudd/.deps/cudd_libcudd_la-cuddDecomp.Tpo cudd/.deps/cudd_libcudd_la-cuddDecomp.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddEssent.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddEssent.Tpo -c -o cudd/cudd_libcudd_la-cuddEssent.lo `test -f 'cudd/cuddEssent.c' || echo './'`cudd/cuddEssent.c
mv -f cudd/.deps/cudd_libcudd_la-cuddEssent.Tpo cudd/.deps/cudd_libcudd_la-cuddEssent.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddExact.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddExact.Tpo -c -o cudd/cudd_libcudd_la-cuddExact.lo `test -f 'cudd/cuddExact.c' || echo './'`cudd/cuddExact.c
mv -f cudd/.deps/cudd_libcudd_la-cuddExact.Tpo cudd/.deps/cudd_libcudd_la-cuddExact.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddExport.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddExport.Tpo -c -o cudd/cudd_libcudd_la-cuddExport.lo `test -f 'cudd/cuddExport.c' || echo './'`cudd/cuddExport.c
mv -f cudd/.deps/cudd_libcudd_la-cuddExport.Tpo cudd/.deps/cudd_libcudd_la-cuddExport.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddGenCof.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddGenCof.Tpo -c -o cudd/cudd_libcudd_la-cuddGenCof.lo `test -f 'cudd/cuddGenCof.c' || echo './'`cudd/cuddGenCof.c
mv -f cudd/.deps/cudd_libcudd_la-cuddGenCof.Tpo cudd/.deps/cudd_libcudd_la-cuddGenCof.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddGenetic.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddGenetic.Tpo -c -o cudd/cudd_libcudd_la-cuddGenetic.lo `test -f 'cudd/cuddGenetic.c' || echo './'`cudd/cuddGenetic.c
mv -f cudd/.deps/cudd_libcudd_la-cuddGenetic.Tpo cudd/.deps/cudd_libcudd_la-cuddGenetic.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddGroup.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddGroup.Tpo -c -o cudd/cudd_libcudd_la-cuddGroup.lo `test -f 'cudd/cuddGroup.c' || echo './'`cudd/cuddGroup.c
mv -f cudd/.deps/cudd_libcudd_la-cuddGroup.Tpo cudd/.deps/cudd_libcudd_la-cuddGroup.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddHarwell.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddHarwell.Tpo -c -o cudd/cudd_libcudd_la-cuddHarwell.lo `test -f 'cudd/cuddHarwell.c' || echo './'`cudd/cuddHarwell.c
mv -f cudd/.deps/cudd_libcudd_la-cuddHarwell.Tpo cudd/.deps/cudd_libcudd_la-cuddHarwell.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddInit.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddInit.Tpo -c -o cudd/cudd_libcudd_la-cuddInit.lo `test -f 'cudd/cuddInit.c' || echo './'`cudd/cuddInit.c
mv -f cudd/.deps/cudd_libcudd_la-cuddInit.Tpo cudd/.deps/cudd_libcudd_la-cuddInit.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddInteract.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddInteract.Tpo -c -o cudd/cudd_libcudd_la-cuddInteract.lo `test -f 'cudd/cuddInteract.c' || echo './'`cudd/cuddInteract.c
mv -f cudd/.deps/cudd_libcudd_la-cuddInteract.Tpo cudd/.deps/cudd_libcudd_la-cuddInteract.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddLCache.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddLCache.Tpo -c -o cudd/cudd_libcudd_la-cuddLCache.lo `test -f 'cudd/cuddLCache.c' || echo './'`cudd/cuddLCache.c
mv -f cudd/.deps/cudd_libcudd_la-cuddLCache.Tpo cudd/.deps/cudd_libcudd_la-cuddLCache.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddLevelQ.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddLevelQ.Tpo -c -o cudd/cudd_libcudd_la-cuddLevelQ.lo `test -f 'cudd/cuddLevelQ.c' || echo './'`cudd/cuddLevelQ.c
mv -f cudd/.deps/cudd_libcudd_la-cuddLevelQ.Tpo cudd/.deps/cudd_libcudd_la-cuddLevelQ.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddLinear.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddLinear.Tpo -c -o cudd/cudd_libcudd_la-cuddLinear.lo `test -f 'cudd/cuddLinear.c' || echo './'`cudd/cuddLinear.c
mv -f cudd/.deps/cudd_libcudd_la-cuddLinear.Tpo cudd/.deps/cudd_libcudd_la-cuddLinear.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddLiteral.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddLiteral.Tpo -c -o cudd/cudd_libcudd_la-cuddLiteral.lo `test -f 'cudd/cuddLiteral.c' || echo './'`cudd/cuddLiteral.c
mv -f cudd/.deps/cudd_libcudd_la-cuddLiteral.Tpo cudd/.deps/cudd_libcudd_la-cuddLiteral.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddMatMult.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddMatMult.Tpo -c -o cudd/cudd_libcudd_la-cuddMatMult.lo `test -f 'cudd/cuddMatMult.c' || echo './'`cudd/cuddMatMult.c
mv -f cudd/.deps/cudd_libcudd_la-cuddMatMult.Tpo cudd/.deps/cudd_libcudd_la-cuddMatMult.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddPriority.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddPriority.Tpo -c -o cudd/cudd_libcudd_la-cuddPriority.lo `test -f 'cudd/cuddPriority.c' || echo './'`cudd/cuddPriority.c
mv -f cudd/.deps/cudd_libcudd_la-cuddPriority.Tpo cudd/.deps/cudd_libcudd_la-cuddPriority.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddRead.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddRead.Tpo -c -o cudd/cudd_libcudd_la-cuddRead.lo `test -f 'cudd/cuddRead.c' || echo './'`cudd/cuddRead.c
mv -f cudd/.deps/cudd_libcudd_la-cuddRead.Tpo cudd/.deps/cudd_libcudd_la-cuddRead.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddRef.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddRef.Tpo -c -o cudd/cudd_libcudd_la-cuddRef.lo `test -f 'cudd/cuddRef.c' || echo './'`cudd/cuddRef.c
mv -f cudd/.deps/cudd_libcudd_la-cuddRef.Tpo cudd/.deps/cudd_libcudd_la-cuddRef.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddReorder.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddReorder.Tpo -c -o cudd/cudd_libcudd_la-cuddReorder.lo `test -f 'cudd/cuddReorder.c' || echo './'`cudd/cuddReorder.c
mv -f cudd/.deps/cudd_libcudd_la-cuddReorder.Tpo cudd/.deps/cudd_libcudd_la-cuddReorder.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddSat.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddSat.Tpo -c -o cudd/cudd_libcudd_la-cuddSat.lo `test -f 'cudd/cuddSat.c' || echo './'`cudd/cuddSat.c
mv -f cudd/.deps/cudd_libcudd_la-cuddSat.Tpo cudd/.deps/cudd_libcudd_la-cuddSat.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddSign.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddSign.Tpo -c -o cudd/cudd_libcudd_la-cuddSign.lo `test -f 'cudd/cuddSign.c' || echo './'`cudd/cuddSign.c
mv -f cudd/.deps/cudd_libcudd_la-cuddSign.Tpo cudd/.deps/cudd_libcudd_la-cuddSign.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddSolve.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddSolve.Tpo -c -o cudd/cudd_libcudd_la-cuddSolve.lo `test -f 'cudd/cuddSolve.c' || echo './'`cudd/cuddSolve.c
mv -f cudd/.deps/cudd_libcudd_la-cuddSolve.Tpo cudd/.deps/cudd_libcudd_la-cuddSolve.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddSplit.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddSplit.Tpo -c -o cudd/cudd_libcudd_la-cuddSplit.lo `test -f 'cudd/cuddSplit.c' || echo './'`cudd/cuddSplit.c
mv -f cudd/.deps/cudd_libcudd_la-cuddSplit.Tpo cudd/.deps/cudd_libcudd_la-cuddSplit.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddSubsetHB.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddSubsetHB.Tpo -c -o cudd/cudd_libcudd_la-cuddSubsetHB.lo `test -f 'cudd/cuddSubsetHB.c' || echo './'`cudd/cuddSubsetHB.c
mv -f cudd/.deps/cudd_libcudd_la-cuddSubsetHB.Tpo cudd/.deps/cudd_libcudd_la-cuddSubsetHB.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddSubsetSP.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddSubsetSP.Tpo -c -o cudd/cudd_libcudd_la-cuddSubsetSP.lo `test -f 'cudd/cuddSubsetSP.c' || echo './'`cudd/cuddSubsetSP.c
mv -f cudd/.deps/cudd_libcudd_la-cuddSubsetSP.Tpo cudd/.deps/cudd_libcudd_la-cuddSubsetSP.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddSymmetry.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddSymmetry.Tpo -c -o cudd/cudd_libcudd_la-cuddSymmetry.lo `test -f 'cudd/cuddSymmetry.c' || echo './'`cudd/cuddSymmetry.c
mv -f cudd/.deps/cudd_libcudd_la-cuddSymmetry.Tpo cudd/.deps/cudd_libcudd_la-cuddSymmetry.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddTable.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddTable.Tpo -c -o cudd/cudd_libcudd_la-cuddTable.lo `test -f 'cudd/cuddTable.c' || echo './'`cudd/cuddTable.c
mv -f cudd/.deps/cudd_libcudd_la-cuddTable.Tpo cudd/.deps/cudd_libcudd_la-cuddTable.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddUtil.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddUtil.Tpo -c -o cudd/cudd_libcudd_la-cuddUtil.lo `test -f 'cudd/cuddUtil.c' || echo './'`cudd/cuddUtil.c
mv -f cudd/.deps/cudd_libcudd_la-cuddUtil.Tpo cudd/.deps/cudd_libcudd_la-cuddUtil.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddWindow.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddWindow.Tpo -c -o cudd/cudd_libcudd_la-cuddWindow.lo `test -f 'cudd/cuddWindow.c' || echo './'`cudd/cuddWindow.c
mv -f cudd/.deps/cudd_libcudd_la-cuddWindow.Tpo cudd/.deps/cudd_libcudd_la-cuddWindow.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddCount.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddCount.Tpo -c -o cudd/cudd_libcudd_la-cuddZddCount.lo `test -f 'cudd/cuddZddCount.c' || echo './'`cudd/cuddZddCount.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddCount.Tpo cudd/.deps/cudd_libcudd_la-cuddZddCount.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddFuncs.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddFuncs.Tpo -c -o cudd/cudd_libcudd_la-cuddZddFuncs.lo `test -f 'cudd/cuddZddFuncs.c' || echo './'`cudd/cuddZddFuncs.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddFuncs.Tpo cudd/.deps/cudd_libcudd_la-cuddZddFuncs.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddGroup.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddGroup.Tpo -c -o cudd/cudd_libcudd_la-cuddZddGroup.lo `test -f 'cudd/cuddZddGroup.c' || echo './'`cudd/cuddZddGroup.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddGroup.Tpo cudd/.deps/cudd_libcudd_la-cuddZddGroup.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddIsop.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddIsop.Tpo -c -o cudd/cudd_libcudd_la-cuddZddIsop.lo `test -f 'cudd/cuddZddIsop.c' || echo './'`cudd/cuddZddIsop.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddIsop.Tpo cudd/.deps/cudd_libcudd_la-cuddZddIsop.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddLin.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddLin.Tpo -c -o cudd/cudd_libcudd_la-cuddZddLin.lo `test -f 'cudd/cuddZddLin.c' || echo './'`cudd/cuddZddLin.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddLin.Tpo cudd/.deps/cudd_libcudd_la-cuddZddLin.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddMisc.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddMisc.Tpo -c -o cudd/cudd_libcudd_la-cuddZddMisc.lo `test -f 'cudd/cuddZddMisc.c' || echo './'`cudd/cuddZddMisc.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddMisc.Tpo cudd/.deps/cudd_libcudd_la-cuddZddMisc.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddPort.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddPort.Tpo -c -o cudd/cudd_libcudd_la-cuddZddPort.lo `test -f 'cudd/cuddZddPort.c' || echo './'`cudd/cuddZddPort.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddPort.Tpo cudd/.deps/cudd_libcudd_la-cuddZddPort.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddReord.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddReord.Tpo -c -o cudd/cudd_libcudd_la-cuddZddReord.lo `test -f 'cudd/cuddZddReord.c' || echo './'`cudd/cuddZddReord.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddReord.Tpo cudd/.deps/cudd_libcudd_la-cuddZddReord.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddSetop.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddSetop.Tpo -c -o cudd/cudd_libcudd_la-cuddZddSetop.lo `test -f 'cudd/cuddZddSetop.c' || echo './'`cudd/cuddZddSetop.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddSetop.Tpo cudd/.deps/cudd_libcudd_la-cuddZddSetop.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddSymm.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddSymm.Tpo -c -o cudd/cudd_libcudd_la-cuddZddSymm.lo `test -f 'cudd/cuddZddSymm.c' || echo './'`cudd/cuddZddSymm.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddSymm.Tpo cudd/.deps/cudd_libcudd_la-cuddZddSymm.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT cudd/cudd_libcudd_la-cuddZddUtil.lo -MD -MP -MF cudd/.deps/cudd_libcudd_la-cuddZddUtil.Tpo -c -o cudd/cudd_libcudd_la-cuddZddUtil.lo `test -f 'cudd/cuddZddUtil.c' || echo './'`cudd/cuddZddUtil.c
mv -f cudd/.deps/cudd_libcudd_la-cuddZddUtil.Tpo cudd/.deps/cudd_libcudd_la-cuddZddUtil.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-cpu_stats.lo -MD -MP -MF util/.deps/cudd_libcudd_la-cpu_stats.Tpo -c -o util/cudd_libcudd_la-cpu_stats.lo `test -f 'util/cpu_stats.c' || echo './'`util/cpu_stats.c
mv -f util/.deps/cudd_libcudd_la-cpu_stats.Tpo util/.deps/cudd_libcudd_la-cpu_stats.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-cpu_time.lo -MD -MP -MF util/.deps/cudd_libcudd_la-cpu_time.Tpo -c -o util/cudd_libcudd_la-cpu_time.lo `test -f 'util/cpu_time.c' || echo './'`util/cpu_time.c
mv -f util/.deps/cudd_libcudd_la-cpu_time.Tpo util/.deps/cudd_libcudd_la-cpu_time.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-cstringstream.lo -MD -MP -MF util/.deps/cudd_libcudd_la-cstringstream.Tpo -c -o util/cudd_libcudd_la-cstringstream.lo `test -f 'util/cstringstream.c' || echo './'`util/cstringstream.c
mv -f util/.deps/cudd_libcudd_la-cstringstream.Tpo util/.deps/cudd_libcudd_la-cstringstream.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-datalimit.lo -MD -MP -MF util/.deps/cudd_libcudd_la-datalimit.Tpo -c -o util/cudd_libcudd_la-datalimit.lo `test -f 'util/datalimit.c' || echo './'`util/datalimit.c
mv -f util/.deps/cudd_libcudd_la-datalimit.Tpo util/.deps/cudd_libcudd_la-datalimit.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-pathsearch.lo -MD -MP -MF util/.deps/cudd_libcudd_la-pathsearch.Tpo -c -o util/cudd_libcudd_la-pathsearch.lo `test -f 'util/pathsearch.c' || echo './'`util/pathsearch.c
mv -f util/.deps/cudd_libcudd_la-pathsearch.Tpo util/.deps/cudd_libcudd_la-pathsearch.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-pipefork.lo -MD -MP -MF util/.deps/cudd_libcudd_la-pipefork.Tpo -c -o util/cudd_libcudd_la-pipefork.lo `test -f 'util/pipefork.c' || echo './'`util/pipefork.c
mv -f util/.deps/cudd_libcudd_la-pipefork.Tpo util/.deps/cudd_libcudd_la-pipefork.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-prtime.lo -MD -MP -MF util/.deps/cudd_libcudd_la-prtime.Tpo -c -o util/cudd_libcudd_la-prtime.lo `test -f 'util/prtime.c' || echo './'`util/prtime.c
mv -f util/.deps/cudd_libcudd_la-prtime.Tpo util/.deps/cudd_libcudd_la-prtime.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-safe_mem.lo -MD -MP -MF util/.deps/cudd_libcudd_la-safe_mem.Tpo -c -o util/cudd_libcudd_la-safe_mem.lo `test -f 'util/safe_mem.c' || echo './'`util/safe_mem.c
mv -f util/.deps/cudd_libcudd_la-safe_mem.Tpo util/.deps/cudd_libcudd_la-safe_mem.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-strsav.lo -MD -MP -MF util/.deps/cudd_libcudd_la-strsav.Tpo -c -o util/cudd_libcudd_la-strsav.lo `test -f 'util/strsav.c' || echo './'`util/strsav.c
mv -f util/.deps/cudd_libcudd_la-strsav.Tpo util/.deps/cudd_libcudd_la-strsav.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-texpand.lo -MD -MP -MF util/.deps/cudd_libcudd_la-texpand.Tpo -c -o util/cudd_libcudd_la-texpand.lo `test -f 'util/texpand.c' || echo './'`util/texpand.c
mv -f util/.deps/cudd_libcudd_la-texpand.Tpo util/.deps/cudd_libcudd_la-texpand.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT util/cudd_libcudd_la-ucbqsort.lo -MD -MP -MF util/.deps/cudd_libcudd_la-ucbqsort.Tpo -c -o util/cudd_libcudd_la-ucbqsort.lo `test -f 'util/ucbqsort.c' || echo './'`util/ucbqsort.c
mv -f util/.deps/cudd_libcudd_la-ucbqsort.Tpo util/.deps/cudd_libcudd_la-ucbqsort.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT st/cudd_libcudd_la-st.lo -MD -MP -MF st/.deps/cudd_libcudd_la-st.Tpo -c -o st/cudd_libcudd_la-st.lo `test -f 'st/st.c' || echo './'`st/st.c
mv -f st/.deps/cudd_libcudd_la-st.Tpo st/.deps/cudd_libcudd_la-st.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT epd/cudd_libcudd_la-epd.lo -MD -MP -MF epd/.deps/cudd_libcudd_la-epd.Tpo -c -o epd/cudd_libcudd_la-epd.lo `test -f 'epd/epd.c' || echo './'`epd/epd.c
mv -f epd/.deps/cudd_libcudd_la-epd.Tpo epd/.deps/cudd_libcudd_la-epd.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT mtr/cudd_libcudd_la-mtrBasic.lo -MD -MP -MF mtr/.deps/cudd_libcudd_la-mtrBasic.Tpo -c -o mtr/cudd_libcudd_la-mtrBasic.lo `test -f 'mtr/mtrBasic.c' || echo './'`mtr/mtrBasic.c
mv -f mtr/.deps/cudd_libcudd_la-mtrBasic.Tpo mtr/.deps/cudd_libcudd_la-mtrBasic.Plo
/bin/sh ./libtool  --tag=CC   --mode=compile gcc -DHAVE_CONFIG_H -I.  -I./cudd -I./st -I./epd -I./mtr -I./util   -Wall -Wextra -g -O3 -MT mtr/cudd_libcudd_la-mtrGroup.lo -MD -MP -MF mtr/.deps/cudd_libcudd_la-mtrGroup.Tpo -c -o mtr/cudd_libcudd_la-mtrGroup.lo `test -f 'mtr/mtrGroup.c' || echo './'`mtr/mtrGroup.c
mv -f mtr/.deps/cudd_libcudd_la-mtrGroup.Tpo mtr/.deps/cudd_libcudd_la-mtrGroup.Plo
/bin/sh ./libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3 -release 3.0.0 -version-info 0:0:0 -no-undefined  -o cudd/libcudd.la -rpath /usr/local/lib cudd/cudd_libcudd_la-cuddAddAbs.lo cudd/cudd_libcudd_la-cuddAddApply.lo cudd/cudd_libcudd_la-cuddAddFind.lo cudd/cudd_libcudd_la-cuddAddInv.lo cudd/cudd_libcudd_la-cuddAddIte.lo cudd/cudd_libcudd_la-cuddAddNeg.lo cudd/cudd_libcudd_la-cuddAddWalsh.lo cudd/cudd_libcudd_la-cuddAndAbs.lo cudd/cudd_libcudd_la-cuddAnneal.lo cudd/cudd_libcudd_la-cuddApa.lo cudd/cudd_libcudd_la-cuddAPI.lo cudd/cudd_libcudd_la-cuddApprox.lo cudd/cudd_libcudd_la-cuddBddAbs.lo cudd/cudd_libcudd_la-cuddBddCorr.lo cudd/cudd_libcudd_la-cuddBddIte.lo cudd/cudd_libcudd_la-cuddBridge.lo cudd/cudd_libcudd_la-cuddCache.lo cudd/cudd_libcudd_la-cuddCheck.lo cudd/cudd_libcudd_la-cuddClip.lo cudd/cudd_libcudd_la-cuddCof.lo cudd/cudd_libcudd_la-cuddCompose.lo cudd/cudd_libcudd_la-cuddDecomp.lo cudd/cudd_libcudd_la-cuddEssent.lo cudd/cudd_libcudd_la-cuddExact.lo cudd/cudd_libcudd_la-cuddExport.lo cudd/cudd_libcudd_la-cuddGenCof.lo cudd/cudd_libcudd_la-cuddGenetic.lo cudd/cudd_libcudd_la-cuddGroup.lo cudd/cudd_libcudd_la-cuddHarwell.lo cudd/cudd_libcudd_la-cuddInit.lo cudd/cudd_libcudd_la-cuddInteract.lo cudd/cudd_libcudd_la-cuddLCache.lo cudd/cudd_libcudd_la-cuddLevelQ.lo cudd/cudd_libcudd_la-cuddLinear.lo cudd/cudd_libcudd_la-cuddLiteral.lo cudd/cudd_libcudd_la-cuddMatMult.lo cudd/cudd_libcudd_la-cuddPriority.lo cudd/cudd_libcudd_la-cuddRead.lo cudd/cudd_libcudd_la-cuddRef.lo cudd/cudd_libcudd_la-cuddReorder.lo cudd/cudd_libcudd_la-cuddSat.lo cudd/cudd_libcudd_la-cuddSign.lo cudd/cudd_libcudd_la-cuddSolve.lo cudd/cudd_libcudd_la-cuddSplit.lo cudd/cudd_libcudd_la-cuddSubsetHB.lo cudd/cudd_libcudd_la-cuddSubsetSP.lo cudd/cudd_libcudd_la-cuddSymmetry.lo cudd/cudd_libcudd_la-cuddTable.lo cudd/cudd_libcudd_la-cuddUtil.lo cudd/cudd_libcudd_la-cuddWindow.lo cudd/cudd_libcudd_la-cuddZddCount.lo cudd/cudd_libcudd_la-cuddZddFuncs.lo cudd/cudd_libcudd_la-cuddZddGroup.lo cudd/cudd_libcudd_la-cuddZddIsop.lo cudd/cudd_libcudd_la-cuddZddLin.lo cudd/cudd_libcudd_la-cuddZddMisc.lo cudd/cudd_libcudd_la-cuddZddPort.lo cudd/cudd_libcudd_la-cuddZddReord.lo cudd/cudd_libcudd_la-cuddZddSetop.lo cudd/cudd_libcudd_la-cuddZddSymm.lo cudd/cudd_libcudd_la-cuddZddUtil.lo util/cudd_libcudd_la-cpu_stats.lo util/cudd_libcudd_la-cpu_time.lo util/cudd_libcudd_la-cstringstream.lo util/cudd_libcudd_la-datalimit.lo util/cudd_libcudd_la-pathsearch.lo util/cudd_libcudd_la-pipefork.lo util/cudd_libcudd_la-prtime.lo util/cudd_libcudd_la-safe_mem.lo util/cudd_libcudd_la-strsav.lo util/cudd_libcudd_la-texpand.lo util/cudd_libcudd_la-ucbqsort.lo st/cudd_libcudd_la-st.lo epd/cudd_libcudd_la-epd.lo mtr/cudd_libcudd_la-mtrBasic.lo mtr/cudd_libcudd_la-mtrGroup.lo    
/bin/sh ./libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o cplusplus/testobj cplusplus/cplusplus_testobj-testobj.o cplusplus/libobj.la cudd/libcudd.la 