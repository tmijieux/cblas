#ifndef TDP_CBLAS_H
#define TDP_CBLAS_H

#include "ddot.h"
#include "dgemm.h"
#include "daxpy.h"

#define ALIGNED(X_) __attribute__(( aligned(X_) ))
#define IS_ALIGNED(ADDR_, ALIGN_) ((((long)(ADDR_)) & ((ALIGN_)-1)) == 0)

#endif // TDP_CBLAS_H
