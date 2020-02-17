#!/bin/bash

INFILE=${CSEISDIR_SRCROOT}/src/include/cseis_modules.txt
TMPFILE=${LIBDIR}/include/tmpfile
TMPFILE_SORT=${LIBDIR}/include/tmpfile.sort
OUTFILE=${LIBDIR}/include/cseis_modules.h
ALLFILE=${LIBDIR}/include/cseis_modules_all.h

echo "...creating module header file $OUTFILE"

modulesAll=$(cat $INFILE | awk '{if(NF>=1 && substr($0,1,1) != "#")print $1}')
numAll=$(cat $INFILE | awk '{if(NF>=1 && substr($0,1,1) != "#") N=N+1}END{print N}' | bc -l)

versionsAll=$(cat $INFILE | awk '{if(NF>=1&&substr($0,1,1)!="#"){ if(NF>=2){print $2" "}else{print "1.0 "}}}')

counter=$(echo "0"|bc -l)
for module in ${modulesAll}
do
    counter=$(echo "$counter + 1"|bc -l)
    moduleUpper=$( echo "$module" | tr '[:lower:]' '[:upper:]' )
    moduleArray[$counter]=$moduleUpper
done
counter=$(echo "0"|bc -l)
for version in ${versionsAll}
do
    counter=$(echo "$counter + 1"|bc -l)
    versionArray[$counter]=$version
done
counter=$(echo "0"|bc -l)

# Write all modules and associated versions into temporary file (purpose: Sorting by name)
if [ -e $TMPFILE ]; then
    rm -f $TMPFILE
fi
touch $TMPFILE
for (( i = 1; i <= $numAll; i = i+1 ))
do
    echo "${moduleArray[i]} ${versionArray[i]}" >> $TMPFILE
done
sort $TMPFILE > ${TMPFILE_SORT}

#-------------------------------------------------------------
# Create header file defining all standard module names
#

echo "#ifndef CS_MODULES_H" > ${OUTFILE}
echo "#define CS_MODULES_H" >> ${OUTFILE}
echo "" >> ${OUTFILE}
echo "" >> ${OUTFILE}
echo "const int N_METHODS  = ${numAll};" >> ${OUTFILE}
echo "" >> ${OUTFILE}
echo "std::string NAMES[N_METHODS] = {" >> ${OUTFILE}

awk -v numModules=$numAll '{printf"  std::string(\""$1"\")";if(NR<numModules){print","}else{print ""}}' ${TMPFILE_SORT} >> ${OUTFILE}

echo "};" >> ${OUTFILE}
echo "" >> ${OUTFILE}
echo "#endif" >> ${OUTFILE}

#*************************************************************
# Create second header file defining all standard modules and associated methods
# ...required for Windows monolitic executable build
#

echo "...creating module header file $ALLFILE"

echo "#ifndef CS_MODULES_ALL_H" > ${ALLFILE}
echo  "#define CS_MODULES_ALL_H" >> ${ALLFILE}
echo "" >> ${ALLFILE}
echo "namespace cseis_system {" >> ${ALLFILE}
echo "  class csTrace;" >> ${ALLFILE}
echo "  class csParamManager;" >> ${ALLFILE}
echo "  class csTraceGather;" >> ${ALLFILE}
echo "  class csParamDef;" >> ${ALLFILE}
echo "  class csInitPhaseEnv;" >> ${ALLFILE}
echo "  class csExecPhaseEnv;" >> ${ALLFILE}
echo "  class csLogWriter;" >> ${ALLFILE}
echo "}" >> ${ALLFILE}
echo "" >> ${ALLFILE}
echo "const int N_METHODS  = $numAll;" >> ${ALLFILE}
echo "" >> ${ALLFILE}

awk '{printf"void init_mod_"tolower($1)"_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);\n"}' ${TMPFILE_SORT} >> ${ALLFILE}
echo "" >> ${ALLFILE}


echo "bool start_exec_mod_dummy_(cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) { return true; }" >> ${ALLFILE}
awk '{printf("bool start_exec_mod_"tolower($1)"_(cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);\n")}' ${TMPFILE_SORT} >> ${ALLFILE}
echo "" >> ${ALLFILE}

echo "void exec_mod_dummy_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) {}" >> ${ALLFILE}
awk '{printf("void exec_mod_"tolower($1)"_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);\n")}' ${TMPFILE_SORT} >> ${ALLFILE}
echo "" >> ${ALLFILE}

echo "void cleanup_mod_dummy_(cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) {}" >> ${ALLFILE}
awk '{printf("void cleanup_mod_"tolower($1)"_(cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);\n")}' ${TMPFILE_SORT} >> ${ALLFILE}
echo "" >> ${ALLFILE}


echo "void(*METHODS_INIT[N_METHODS])( cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter* ) = {" >> ${ALLFILE}
awk '{if(NR>1){print ","};     printf("  init_mod_"tolower($1)"_")}' ${TMPFILE_SORT} >> ${ALLFILE}
echo ""  >> ${ALLFILE}
echo "};" >> ${ALLFILE}
echo "" >> ${ALLFILE}

echo "bool(*METHODS_EXEC_START[N_METHODS])(cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) = {" >> ${ALLFILE}
awk '{if(NR>1){print ","};     printf("  start_exec_mod_"tolower($1)"_")}' ${TMPFILE_SORT} >> ${ALLFILE}
echo ""  >> ${ALLFILE}
echo "};" >> ${ALLFILE}
echo "" >> ${ALLFILE}

echo "void(*METHODS_EXEC[N_METHODS])(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) = {" >> ${ALLFILE}
awk '{if(NR>1){print ","};     printf("  exec_mod_"tolower($1)"_")}' ${TMPFILE_SORT} >> ${ALLFILE}
echo ""  >> ${ALLFILE}
echo "};" >> ${ALLFILE}
echo "" >> ${ALLFILE}

echo "void(*METHODS_CLEANUP[N_METHODS])(cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) = {" >> ${ALLFILE}
awk '{if(NR>1){print ","};     printf("  cleanup_mod_"tolower($1)"_")}' ${TMPFILE_SORT} >> ${ALLFILE}
echo ""  >> ${ALLFILE}
echo "};" >> ${ALLFILE}
echo "" >> ${ALLFILE}
echo "std::string NAMES[N_METHODS] = {" >> ${ALLFILE}


awk '{if(NR>1){print ","};N=1;printf("  std::string(\""toupper($1)"\")")}' ${TMPFILE_SORT} >> ${ALLFILE}
echo ""  >> ${ALLFILE}
echo "};" >> ${ALLFILE}
echo "" >> ${ALLFILE}


awk '{printf("void params_mod_"tolower($1)"_(cseis_system::csParamDef*);\n")}' ${TMPFILE_SORT} >> ${ALLFILE}
echo "" >> ${ALLFILE}
echo "void (*METHODS_PARAM[N_METHODS])(cseis_system::csParamDef*) = {" >> ${ALLFILE}

awk '{if(NR>1){print ","};N=1;printf("  params_mod_"tolower($1)"_")}' ${TMPFILE_SORT} >> ${ALLFILE}
echo ""  >> ${ALLFILE}
echo "};" >> ${ALLFILE}
echo ""  >> ${ALLFILE}
echo "#endif"  >> ${ALLFILE}

#*************************************************************
# Create soft links to all standard modules and system libraries
#

echo "...creating soft links for all libraries in directory $LIBDIR"

#for module in ${modulesAll}
for (( i = 1; i <= $numAll; i = i+1 ))
do
    module=${moduleArray[i]}
    echo "..."$module
    version=${versionArray[i]}
    moduleLower=$( echo "$module" | tr '[:upper:]' '[:lower:]' )
    libname1=${LIBDIR}/libmod_$moduleLower.so
    libname2=${libname1}.${version}
    rm -f $libname1
#    rm -f $libname2
    ln -s $libname2 $libname1
done

echo "...DONE creating soft links for all libraries in directory $LIBDIR"


#rm -f $TMPFILE
#rm -f ${TMPFILE_SORT}
