/* Copyright (c) 2014-2018 The Bitcoin Core developers
 * Copyright (c) 2022 Mark Friedenbach
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SHA2__SHA256_INTERNAL_H
#define SHA2__SHA256_INTERNAL_H

#include <sha2/sha256.h>

#if defined(__x86_64__) || defined(__amd64__) || defined(__i386__)
extern void transform_sha256_sse4(uint32_t* s, const unsigned char* chunk, size_t blocks);

extern void transform_sha256multi_sse41_4way(struct sha256* out, const uint32_t* s, const unsigned char* in);
extern void transform_sha256d64_sse41_4way(struct sha256 out[4], const struct sha256 in[8]);

extern void transform_sha256multi_avx2_8way(struct sha256* out, const uint32_t* s, const unsigned char* in);
extern void transform_sha256d64_avx2_8way(struct sha256 out[8], const struct sha256 in[16]);

extern void transform_sha256_shani(uint32_t* s, const unsigned char* chunk, size_t blocks);
extern void transform_sha256d64_shani_2way(struct sha256 out[2], const struct sha256 in[4]);
#endif
#if defined(__arm__) || defined(__aarch32__) || defined(__arm64__) || defined(__aarch64__) || defined(_M_ARM)
extern void transform_sha256_armv8(uint32_t* s, const unsigned char* chunk, size_t blocks);
extern void transform_sha256d64_armv8_2way(struct sha256 out[2], const struct sha256 in[4]);
#endif

#endif /* SHA2__SHA256_INTERNAL_H */

/* End of File
 */
