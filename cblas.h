#ifndef TDP_CBLAS_H
#define TDP_CBLAS_H

#include "ddot.h"
#include "dgemm.h"
#include "dgemv.h"
#include "daxpy.h"

#define ALIGNED(X_) __attribute__(( aligned(X_) ))
#define IS_ALIGNED(ADDR_, ALIGN_) ((((long)(ADDR_)) & ((ALIGN_)-1)) == 0)
#define UPPER_LD(VALUE_, ALIGN_) (((VALUE_) + ((ALIGN_)-1)) & ~((ALIGN_)-1))

#endif // TDP_CBLAS_H
