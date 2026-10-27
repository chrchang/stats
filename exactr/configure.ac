AC_INIT([exactr],[0.1.0],[chrchang@alumni.caltech.edu])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([could not determine R_HOME])
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

dnl The double-double arithmetic in plink2_highprec requires that the compiler
dnl not contract x*y+z into fma(x,y,z) on its own; see the comments in
dnl include/plink2_float.h.  -ffp-contract=off isn't portable, so only use it
dnl when the compiler accepts it.
AC_MSG_CHECKING([whether ${CXX} accepts -ffp-contract=off])
fp_contract_save_CXXFLAGS="${CXXFLAGS}"
CXXFLAGS="${CXXFLAGS} -ffp-contract=off -Werror"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [])],
                  [FP_CONTRACT_FLAGS="-ffp-contract=off"
                   AC_MSG_RESULT([yes])],
                  [FP_CONTRACT_FLAGS=""
                   AC_MSG_RESULT([no])])
CXXFLAGS="${fp_contract_save_CXXFLAGS}"
AC_SUBST([FP_CONTRACT_FLAGS])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
