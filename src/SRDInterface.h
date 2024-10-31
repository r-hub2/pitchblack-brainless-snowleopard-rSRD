#ifndef SRD_INTERFACE_H_
#define SRD_INTERFACE_H_

#include <Rcpp.h>

#if defined(R_BUILD)
 #define STRICT_R_HEADERS
 #include "R.h"
 // textual substitution
 #define std::cout Rcpp::Rcout
#endif

#endif /* SRD_INTERFACE_H_ */
