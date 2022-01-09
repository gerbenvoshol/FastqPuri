#ifndef _FNV_1A_H
#define _FNV_1A_H

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

uint64_t fnv_1a(const char *key, int len, int seed);

#ifdef __cplusplus
}
#endif

#endif