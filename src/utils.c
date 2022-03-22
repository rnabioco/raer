#include <Rinternals.h>

// https://stackoverflow.com/questions/26666614/how-do-i-check-if-an-externalptr-is-null-from-within-r
SEXP isnull(SEXP pointer) {
  return ScalarLogical(!R_ExternalPtrAddr(pointer));
}
