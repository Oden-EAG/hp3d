NOTE: This may not compile due to uweak_module having matrices that are too 
large. This coud be avoided with dynamic allocation but requires more 
substantial modifications. The simple fix is to change the parameter 
MAXNORD_ADD from 3 to 2 in the module parametersDPG.F90 
(src/modules/parametersDPG.F90). After making this change remember to
recompile the library, since this is not a local file. 

