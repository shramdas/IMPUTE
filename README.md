=======
///////    Project of Biostatistics 830  /////////////

C++ code to do imputation of hyplotypes with a HMM...
Each job has its own head file *.h

LoadVCF.h contains functions to load study/reference haplotypes.

CalcGL.h contains three functions (InitializeFirstVector, Transpose, Condition) to be used in HMM LeftWalk & RightWalk.
Please include CalcGL.h if need to call those three functions by:  #include "CalcGL.h"

HMM.h contains functions for HMM walk from left/right/combine.

Impute.h contains functions to do imputation.


compile cpp code :

g++ -I ../libStatGen/include/ -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS main.cpp ../libStatGen/libStatGen.a -lm -lz -o impute

