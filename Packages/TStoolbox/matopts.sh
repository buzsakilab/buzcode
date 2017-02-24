#
# matopts.sh	Shell script for configuring MAT-file standalone applications.
#               These options were tested with the specified compiler.
#
# usage:        Do not call this file directly; it is sourced by the
#               mbuild shell script.  Modify only if you don't like the
#               defaults after running mbuild.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_SA_OPT: Template Options file for building standalone MAT applications
#
# Copyright 1984-2004 The MathWorks, Inc.
# $Revision: 1.30.4.5 $  $Date: 2004/04/25 21:30:54 $
#----------------------------------------------------------------------------
#
    if [ "$TMW_ROOT" = "" ]; then
	TMW_ROOT="$MATLAB"
    fi
    MFLAGS="-I$TMW_ROOT/extern/include"
    MLIBS="-L$TMW_ROOT/bin/$Arch -lmat -lmx"
    LDEXTENSION=''
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The mex script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        hpux)
#----------------------------------------------------------------------------
#           what `which cc`
#           HP92453-01 B.11.11.06 HP C Compiler
            CC='cc'
            COMPFLAGS='+DA2.0'
            CFLAGS="-Ae $COMPFLAGS $MFLAGS -Wp,-H65535"
            CLIBS="$MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           what `which aCC`
#           HP aC++ B3910B A.03.37
#           HP aC++ B3910B A.03.30 Language Support Library
            CXX='aCC'
            CXXFLAGS="$MCXXFLAGS -AA -D_HPUX_SOURCE $COMPFLAGS"
            CXXLIBS="$MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG +Oconservative'
            CXXDEBUGFLAGS='-g'
#
#           what `which f90`
#           HP-UX f90 20020606 (083554)  B3907DB/B3909DB B.11.01.60
#           HP F90 v2.6
#            $ PATCH/11.00:PHCO_95167  Oct  1 1998 13:46:32 $
            FC='f90'
            FFLAGS="$MFLAGS $COMPFLAGS"
            FLIBS="$MLIBS -lm"
            FOPTIMFLAGS='-O +Oconservative'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="$COMPFLAGS"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        glnx86)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
#           gcc -v
#           gcc version 3.2.3
            CC='gcc'
            CFLAGS="$MFLAGS -ansi -D_GNU_SOURCE -fexceptions"
            CLIBS="$RPATH $MLIBS -lm -lstdc++"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           g77 -fversion
#           GNU Fortran (GCC 3.2.3) 3.2.3 20030422 (release)
#           NOTE: g77 is not thread safe
            FC='g77'
            FFLAGS="$MFLAGS -fexceptions"
            FLIBS="$RPATH $MLIBS -lm -lstdc++"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS=''
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        glnxi64)
#----------------------------------------------------------------------------
echo "Error: Did not imbed 'options.sh' code"; exit 1 #imbed options.sh glnxi64 12
#----------------------------------------------------------------------------
            ;;
        glnxa64)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
#           gcc -v
#           gcc version 3.2.3
            CC='gcc'
            CFLAGS="$MFLAGS -ansi -D_GNU_SOURCE -fexceptions"
            CLIBS="$RPATH $MLIBS -lm -lstdc++"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           g77 -fversion
#           GNU Fortran (GCC 3.2.3) 3.2.3 20030422 (release)
#           NOTE: g77 is not thread safe
            FC='g77'
            FFLAGS="$MFLAGS -fexceptions"
            FLIBS="$RPATH $MLIBS -lm -lstdc++"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS=''
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sol2)
#----------------------------------------------------------------------------
#           cc -V
#           Sun C 5.5 Patch 112760-06 2004/01/13
            CC='cc'
            CFLAGS="$MFLAGS -dalign -xlibmieee -D__EXTENSIONS__"
            CLIBS="$MLIBS -lm"
            COPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CDEBUGFLAGS='-g'  
            CXXDEBUGFLAGS='-g'
#
#           f90 -V
#           Sun Fortran 95 7.1 Patch 112762-09 2004/01/26
            FC='f90'
            FFLAGS="$MFLAGS -dalign"
            FLIBS="$MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS=''
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'  
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        mac)
#----------------------------------------------------------------------------
#           gcc-3.3 -v
#           gcc version 3.3 20030304 (Apple Computer, Inc. build 1435)
            CC='gcc-3.3'
            CFLAGS="$MFLAGS -fno-common -no-cpp-precomp -fexceptions"
            CLIBS="$MLIBS -lstdc++"
            COPTIMFLAGS='-O3 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           f77 -V
#           FORTRAN 77 Compiler 8.2a
            FC='f77'
            FFLAGS="-f -N15 -N11 -s -Q51 -W $MFLAGS"
            ABSOFTLIBDIR=`which $FC | sed -n -e '1s|bin/'$FC'|lib|p'`
            FLIBS="-L$ABSOFTLIBDIR -lfio -lf77math"
            FLIBS="$MLIBS -L$ABSOFTLIBDIR -lfio -lf77math -lstdc++"
            FOPTIMFLAGS='-O -cpu:g4'
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDFLAGS='-Wl,-flat_namespace -undefined suppress'
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#           LDEXTENSION="$LDEXTENSION"
#----------------------------------------------------------------------------
#############################################################################
