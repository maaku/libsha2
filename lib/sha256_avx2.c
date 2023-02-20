/* Copyright (c) 2017-2019 The Bitcoin Core developers
 * Copyright (c) 2022 Mark Friedenbach
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#if defined(__x86_64__) || defined(__amd64__)

#include <sha2/sha256.h>
#include "sha256_internal.h"

#include <stdint.h> /* for uint32_t */
#include <immintrin.h> /* for assembly intrinsics */

#include "common.h"

/* In order to keep the code compatible with the C89 standard, inline functions
 * are to be avoided.  Therefore most of the following inlines are implemented
 * as preprocessor macros.  To prevent bugs due to multiple expansion, macros
 * which reference a parameter more than once are kept as inlines, so if using a
 * compiler which does not understand the C99 automatic inline attribute, the
 * `inline` keyword will need to be removed with a macro, like so:.
 */

/* #define inline */

/* This may change performance slightly as automatic inlining will no longer be
 * performed for these functions.
 */

#define K(x) _mm256_set1_epi32(x)

#define Add(x, y) _mm256_add_epi32((x), (y))
#define Add3(x, y, z) Add(Add((x), (y)), (z))
#define Add4(x, y, z, w) Add(Add((x), (y)), Add((z), (w)))
#define Add5(x, y, z, w, v) Add(Add3((x), (y), (z)), Add((w), (v)))
static inline __attribute__((always_inline)) __m256i Inc_avx2(__m256i *x, __m256i y) { *x = Add(*x, y); return *x; }
static inline __attribute__((always_inline)) __m256i Inc3_avx2(__m256i *x, __m256i y, __m256i z) { *x = Add3(*x, y, z); return *x; }
static inline __attribute__((always_inline)) __m256i Inc4_avx2(__m256i *x, __m256i y, __m256i z, __m256i w) { *x = Add4(*x, y, z, w); return *x; }
#define Xor(x, y) _mm256_xor_si256((x), (y))
#define Xor3(x, y, z) Xor(Xor((x), (y)), (z))
#define Or(x, y) _mm256_or_si256((x), (y))
#define And(x, y) _mm256_and_si256((x), (y))
#define ShR(x, n) _mm256_srli_epi32((x), (n))
#define ShL(x, n) _mm256_slli_epi32((x), (n))

static inline __attribute__((always_inline)) __m256i Ch_avx2(__m256i x, __m256i y, __m256i z) { return Xor(z, And(x, Xor(y, z))); }
static inline __attribute__((always_inline)) __m256i Maj_avx2(__m256i x, __m256i y, __m256i z) { return Or(And(x, y), And(z, Or(x, y))); }
static inline __attribute__((always_inline)) __m256i Sigma0_avx2(__m256i x) { return Xor3(Or(ShR(x, 2), ShL(x, 30)), Or(ShR(x, 13), ShL(x, 19)), Or(ShR(x, 22), ShL(x, 10))); }
static inline __attribute__((always_inline)) __m256i Sigma1_avx2(__m256i x) { return Xor3(Or(ShR(x, 6), ShL(x, 26)), Or(ShR(x, 11), ShL(x, 21)), Or(ShR(x, 25), ShL(x, 7))); }
static inline __attribute__((always_inline)) __m256i sigma0_avx2(__m256i x) { return Xor3(Or(ShR(x, 7), ShL(x, 25)), Or(ShR(x, 18), ShL(x, 14)), ShR(x, 3)); }
static inline __attribute__((always_inline)) __m256i sigma1_avx2(__m256i x) { return Xor3(Or(ShR(x, 17), ShL(x, 15)), Or(ShR(x, 19), ShL(x, 13)), ShR(x, 10)); }

/** One round of SHA-256. */
static inline __attribute__((always_inline)) void Round_avx2(__m256i a, __m256i b, __m256i c, __m256i *d, __m256i e, __m256i f, __m256i g, __m256i *h, __m256i k)
{
        __m256i t1 = Add4(*h, Sigma1_avx2(e), Ch_avx2(e, f, g), k);
        __m256i t2 = Add(Sigma0_avx2(a), Maj_avx2(a, b, c));
        *d = Add(*d, t1);
        *h = Add(t1, t2);
}

static inline __attribute__((always_inline)) __m256i Read8_avx2(const unsigned char* chunk)
{
        return _mm256_shuffle_epi8(
                _mm256_set_epi32(
                        ReadLE32(chunk + 0),
                        ReadLE32(chunk + 64),
                        ReadLE32(chunk + 128),
                        ReadLE32(chunk + 192),
                        ReadLE32(chunk + 256),
                        ReadLE32(chunk + 320),
                        ReadLE32(chunk + 384),
                        ReadLE32(chunk + 448)),
                _mm256_set_epi32(
                        0x0C0D0E0FUL, 0x08090A0BUL, 0x04050607UL, 0x00010203UL,
                        0x0C0D0E0FUL, 0x08090A0BUL, 0x04050607UL, 0x00010203UL));
}

static inline __attribute__((always_inline)) void Write8_avx2(unsigned char *out, __m256i v)
{
        v = _mm256_shuffle_epi8(v, _mm256_set_epi32(
                0x0C0D0E0FUL, 0x08090A0BUL, 0x04050607UL, 0x00010203UL,
                0x0C0D0E0FUL, 0x08090A0BUL, 0x04050607UL, 0x00010203UL));
        WriteLE32(out + 0, _mm256_extract_epi32(v, 7));
        WriteLE32(out + 32, _mm256_extract_epi32(v, 6));
        WriteLE32(out + 64, _mm256_extract_epi32(v, 5));
        WriteLE32(out + 96, _mm256_extract_epi32(v, 4));
        WriteLE32(out + 128, _mm256_extract_epi32(v, 3));
        WriteLE32(out + 160, _mm256_extract_epi32(v, 2));
        WriteLE32(out + 192, _mm256_extract_epi32(v, 1));
        WriteLE32(out + 224, _mm256_extract_epi32(v, 0));
}

void transform_sha256multi_avx2_8way(struct sha256* out, const uint32_t* s, const unsigned char* in)
{
        /* Transform 1 */
        __m256i a = K(s[0]);
        __m256i b = K(s[1]);
        __m256i c = K(s[2]);
        __m256i d = K(s[3]);
        __m256i e = K(s[4]);
        __m256i f = K(s[5]);
        __m256i g = K(s[6]);
        __m256i h = K(s[7]);

        __m256i w0 = Read8_avx2(&in[0]),
                w1 = Read8_avx2(&in[4]),
                w2 = Read8_avx2(&in[8]),
                w3 = Read8_avx2(&in[12]),
                w4 = Read8_avx2(&in[16]),
                w5 = Read8_avx2(&in[20]),
                w6 = Read8_avx2(&in[24]),
                w7 = Read8_avx2(&in[28]),
                w8 = Read8_avx2(&in[32]),
                w9 = Read8_avx2(&in[36]),
                w10 = Read8_avx2(&in[40]),
                w11 = Read8_avx2(&in[44]),
                w12 = Read8_avx2(&in[48]),
                w13 = Read8_avx2(&in[52]),
                w14 = Read8_avx2(&in[56]),
                w15 = Read8_avx2(&in[60]);

        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x428a2f98ul), w0));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x71374491ul), w1));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0xb5c0fbcful), w2));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0xe9b5dba5ul), w3));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x3956c25bul), w4));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x59f111f1ul), w5));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x923f82a4ul), w6));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0xab1c5ed5ul), w7));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0xd807aa98ul), w8));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x12835b01ul), w9));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x243185beul), w10));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x550c7dc3ul), w11));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x72be5d74ul), w12));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x80deb1feul), w13));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x9bdc06a7ul), w14));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0xc19bf174ul), w15));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0xe49b69c1ul), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0xefbe4786ul), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x0fc19dc6ul), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x240ca1ccul), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x2de92c6ful), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x4a7484aaul), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x5cb0a9dcul), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x76f988daul), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x983e5152ul), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0xa831c66dul), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0xb00327c8ul), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0xbf597fc7ul), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0xc6e00bf3ul), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0xd5a79147ul), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x06ca6351ul), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x14292967ul), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x27b70a85ul), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x2e1b2138ul), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x4d2c6dfcul), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x53380d13ul), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x650a7354ul), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x766a0abbul), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x81c2c92eul), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x92722c85ul), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0xa2bfe8a1ul), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0xa81a664bul), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0xc24b8b70ul), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0xc76c51a3ul), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0xd192e819ul), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0xd6990624ul), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0xf40e3585ul), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x106aa070ul), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x19a4c116ul), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x1e376c08ul), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x2748774cul), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x34b0bcb5ul), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x391c0cb3ul), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x4ed8aa4aul), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x5b9cca4ful), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x682e6ff3ul), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x748f82eeul), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x78a5636ful), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x84c87814ul), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x8cc70208ul), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x90befffaul), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0xa4506cebul), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0xbef9a3f7ul), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0xc67178f2ul), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));

        /* Output */
        Write8_avx2(&out->u8[0], Add(a, K(s[0])));
        Write8_avx2(&out->u8[4], Add(b, K(s[1])));
        Write8_avx2(&out->u8[8], Add(c, K(s[2])));
        Write8_avx2(&out->u8[12], Add(d, K(s[3])));
        Write8_avx2(&out->u8[16], Add(e, K(s[4])));
        Write8_avx2(&out->u8[20], Add(f, K(s[5])));
        Write8_avx2(&out->u8[24], Add(g, K(s[6])));
        Write8_avx2(&out->u8[28], Add(h, K(s[7])));
}

void transform_sha256d64_avx2_8way(struct sha256 out[8], const struct sha256 in[16])
{
        /* Transform 1 */
        __m256i a = K(0x6a09e667ul);
        __m256i b = K(0xbb67ae85ul);
        __m256i c = K(0x3c6ef372ul);
        __m256i d = K(0xa54ff53aul);
        __m256i e = K(0x510e527ful);
        __m256i f = K(0x9b05688cul);
        __m256i g = K(0x1f83d9abul);
        __m256i h = K(0x5be0cd19ul);

        __m256i w0 = Read8_avx2(&in[0].u8[0]),
                w1 = Read8_avx2(&in[0].u8[4]),
                w2 = Read8_avx2(&in[0].u8[8]),
                w3 = Read8_avx2(&in[0].u8[12]),
                w4 = Read8_avx2(&in[0].u8[16]),
                w5 = Read8_avx2(&in[0].u8[20]),
                w6 = Read8_avx2(&in[0].u8[24]),
                w7 = Read8_avx2(&in[0].u8[28]),
                w8 = Read8_avx2(&in[1].u8[0]),
                w9 = Read8_avx2(&in[1].u8[4]),
                w10 = Read8_avx2(&in[1].u8[8]),
                w11 = Read8_avx2(&in[1].u8[12]),
                w12 = Read8_avx2(&in[1].u8[16]),
                w13 = Read8_avx2(&in[1].u8[20]),
                w14 = Read8_avx2(&in[1].u8[24]),
                w15 = Read8_avx2(&in[1].u8[28]);

        __m256i t0, t1, t2, t3, t4, t5, t6, t7;

        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x428a2f98ul), w0));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x71374491ul), w1));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0xb5c0fbcful), w2));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0xe9b5dba5ul), w3));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x3956c25bul), w4));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x59f111f1ul), w5));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x923f82a4ul), w6));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0xab1c5ed5ul), w7));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0xd807aa98ul), w8));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x12835b01ul), w9));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x243185beul), w10));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x550c7dc3ul), w11));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x72be5d74ul), w12));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x80deb1feul), w13));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x9bdc06a7ul), w14));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0xc19bf174ul), w15));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0xe49b69c1ul), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0xefbe4786ul), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x0fc19dc6ul), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x240ca1ccul), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x2de92c6ful), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x4a7484aaul), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x5cb0a9dcul), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x76f988daul), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x983e5152ul), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0xa831c66dul), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0xb00327c8ul), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0xbf597fc7ul), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0xc6e00bf3ul), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0xd5a79147ul), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x06ca6351ul), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x14292967ul), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x27b70a85ul), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x2e1b2138ul), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x4d2c6dfcul), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x53380d13ul), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x650a7354ul), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x766a0abbul), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x81c2c92eul), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x92722c85ul), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0xa2bfe8a1ul), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0xa81a664bul), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0xc24b8b70ul), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0xc76c51a3ul), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0xd192e819ul), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0xd6990624ul), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0xf40e3585ul), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x106aa070ul), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x19a4c116ul), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x1e376c08ul), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x2748774cul), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x34b0bcb5ul), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x391c0cb3ul), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x4ed8aa4aul), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x5b9cca4ful), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x682e6ff3ul), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x748f82eeul), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x78a5636ful), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x84c87814ul), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x8cc70208ul), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x90befffaul), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0xa4506cebul), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0xbef9a3f7ul), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0xc67178f2ul), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));

        t0 = a = Add(a, K(0x6a09e667ul));
        t1 = b = Add(b, K(0xbb67ae85ul));
        t2 = c = Add(c, K(0x3c6ef372ul));
        t3 = d = Add(d, K(0xa54ff53aul));
        t4 = e = Add(e, K(0x510e527ful));
        t5 = f = Add(f, K(0x9b05688cul));
        t6 = g = Add(g, K(0x1f83d9abul));
        t7 = h = Add(h, K(0x5be0cd19ul));

        /* Transform 2 */
        Round_avx2(a, b, c, &d, e, f, g, &h, K(0xc28a2f98ul));
        Round_avx2(h, a, b, &c, d, e, f, &g, K(0x71374491ul));
        Round_avx2(g, h, a, &b, c, d, e, &f, K(0xb5c0fbcful));
        Round_avx2(f, g, h, &a, b, c, d, &e, K(0xe9b5dba5ul));
        Round_avx2(e, f, g, &h, a, b, c, &d, K(0x3956c25bul));
        Round_avx2(d, e, f, &g, h, a, b, &c, K(0x59f111f1ul));
        Round_avx2(c, d, e, &f, g, h, a, &b, K(0x923f82a4ul));
        Round_avx2(b, c, d, &e, f, g, h, &a, K(0xab1c5ed5ul));
        Round_avx2(a, b, c, &d, e, f, g, &h, K(0xd807aa98ul));
        Round_avx2(h, a, b, &c, d, e, f, &g, K(0x12835b01ul));
        Round_avx2(g, h, a, &b, c, d, e, &f, K(0x243185beul));
        Round_avx2(f, g, h, &a, b, c, d, &e, K(0x550c7dc3ul));
        Round_avx2(e, f, g, &h, a, b, c, &d, K(0x72be5d74ul));
        Round_avx2(d, e, f, &g, h, a, b, &c, K(0x80deb1feul));
        Round_avx2(c, d, e, &f, g, h, a, &b, K(0x9bdc06a7ul));
        Round_avx2(b, c, d, &e, f, g, h, &a, K(0xc19bf374ul));
        Round_avx2(a, b, c, &d, e, f, g, &h, K(0x649b69c1ul));
        Round_avx2(h, a, b, &c, d, e, f, &g, K(0xf0fe4786ul));
        Round_avx2(g, h, a, &b, c, d, e, &f, K(0x0fe1edc6ul));
        Round_avx2(f, g, h, &a, b, c, d, &e, K(0x240cf254ul));
        Round_avx2(e, f, g, &h, a, b, c, &d, K(0x4fe9346ful));
        Round_avx2(d, e, f, &g, h, a, b, &c, K(0x6cc984beul));
        Round_avx2(c, d, e, &f, g, h, a, &b, K(0x61b9411eul));
        Round_avx2(b, c, d, &e, f, g, h, &a, K(0x16f988faul));
        Round_avx2(a, b, c, &d, e, f, g, &h, K(0xf2c65152ul));
        Round_avx2(h, a, b, &c, d, e, f, &g, K(0xa88e5a6dul));
        Round_avx2(g, h, a, &b, c, d, e, &f, K(0xb019fc65ul));
        Round_avx2(f, g, h, &a, b, c, d, &e, K(0xb9d99ec7ul));
        Round_avx2(e, f, g, &h, a, b, c, &d, K(0x9a1231c3ul));
        Round_avx2(d, e, f, &g, h, a, b, &c, K(0xe70eeaa0ul));
        Round_avx2(c, d, e, &f, g, h, a, &b, K(0xfdb1232bul));
        Round_avx2(b, c, d, &e, f, g, h, &a, K(0xc7353eb0ul));
        Round_avx2(a, b, c, &d, e, f, g, &h, K(0x3069bad5ul));
        Round_avx2(h, a, b, &c, d, e, f, &g, K(0xcb976d5ful));
        Round_avx2(g, h, a, &b, c, d, e, &f, K(0x5a0f118ful));
        Round_avx2(f, g, h, &a, b, c, d, &e, K(0xdc1eeefdul));
        Round_avx2(e, f, g, &h, a, b, c, &d, K(0x0a35b689ul));
        Round_avx2(d, e, f, &g, h, a, b, &c, K(0xde0b7a04ul));
        Round_avx2(c, d, e, &f, g, h, a, &b, K(0x58f4ca9dul));
        Round_avx2(b, c, d, &e, f, g, h, &a, K(0xe15d5b16ul));
        Round_avx2(a, b, c, &d, e, f, g, &h, K(0x007f3e86ul));
        Round_avx2(h, a, b, &c, d, e, f, &g, K(0x37088980ul));
        Round_avx2(g, h, a, &b, c, d, e, &f, K(0xa507ea32ul));
        Round_avx2(f, g, h, &a, b, c, d, &e, K(0x6fab9537ul));
        Round_avx2(e, f, g, &h, a, b, c, &d, K(0x17406110ul));
        Round_avx2(d, e, f, &g, h, a, b, &c, K(0x0d8cd6f1ul));
        Round_avx2(c, d, e, &f, g, h, a, &b, K(0xcdaa3b6dul));
        Round_avx2(b, c, d, &e, f, g, h, &a, K(0xc0bbbe37ul));
        Round_avx2(a, b, c, &d, e, f, g, &h, K(0x83613bdaul));
        Round_avx2(h, a, b, &c, d, e, f, &g, K(0xdb48a363ul));
        Round_avx2(g, h, a, &b, c, d, e, &f, K(0x0b02e931ul));
        Round_avx2(f, g, h, &a, b, c, d, &e, K(0x6fd15ca7ul));
        Round_avx2(e, f, g, &h, a, b, c, &d, K(0x521afacaul));
        Round_avx2(d, e, f, &g, h, a, b, &c, K(0x31338431ul));
        Round_avx2(c, d, e, &f, g, h, a, &b, K(0x6ed41a95ul));
        Round_avx2(b, c, d, &e, f, g, h, &a, K(0x6d437890ul));
        Round_avx2(a, b, c, &d, e, f, g, &h, K(0xc39c91f2ul));
        Round_avx2(h, a, b, &c, d, e, f, &g, K(0x9eccabbdul));
        Round_avx2(g, h, a, &b, c, d, e, &f, K(0xb5c9a0e6ul));
        Round_avx2(f, g, h, &a, b, c, d, &e, K(0x532fb63cul));
        Round_avx2(e, f, g, &h, a, b, c, &d, K(0xd2c741c6ul));
        Round_avx2(d, e, f, &g, h, a, b, &c, K(0x07237ea3ul));
        Round_avx2(c, d, e, &f, g, h, a, &b, K(0xa4954b68ul));
        Round_avx2(b, c, d, &e, f, g, h, &a, K(0x4c191d76ul));

        w0 = Add(t0, a);
        w1 = Add(t1, b);
        w2 = Add(t2, c);
        w3 = Add(t3, d);
        w4 = Add(t4, e);
        w5 = Add(t5, f);
        w6 = Add(t6, g);
        w7 = Add(t7, h);

        /* Transform 3 */
        a = K(0x6a09e667ul);
        b = K(0xbb67ae85ul);
        c = K(0x3c6ef372ul);
        d = K(0xa54ff53aul);
        e = K(0x510e527ful);
        f = K(0x9b05688cul);
        g = K(0x1f83d9abul);
        h = K(0x5be0cd19ul);

        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x428a2f98ul), w0));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x71374491ul), w1));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0xb5c0fbcful), w2));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0xe9b5dba5ul), w3));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x3956c25bul), w4));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x59f111f1ul), w5));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x923f82a4ul), w6));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0xab1c5ed5ul), w7));
        Round_avx2(a, b, c, &d, e, f, g, &h, K(0x5807aa98ul));
        Round_avx2(h, a, b, &c, d, e, f, &g, K(0x12835b01ul));
        Round_avx2(g, h, a, &b, c, d, e, &f, K(0x243185beul));
        Round_avx2(f, g, h, &a, b, c, d, &e, K(0x550c7dc3ul));
        Round_avx2(e, f, g, &h, a, b, c, &d, K(0x72be5d74ul));
        Round_avx2(d, e, f, &g, h, a, b, &c, K(0x80deb1feul));
        Round_avx2(c, d, e, &f, g, h, a, &b, K(0x9bdc06a7ul));
        Round_avx2(b, c, d, &e, f, g, h, &a, K(0xc19bf274ul));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0xe49b69c1ul), Inc_avx2(&w0, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0xefbe4786ul), Inc3_avx2(&w1, K(0xa00000ul), sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x0fc19dc6ul), Inc3_avx2(&w2, sigma1_avx2(w0), sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x240ca1ccul), Inc3_avx2(&w3, sigma1_avx2(w1), sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x2de92c6ful), Inc3_avx2(&w4, sigma1_avx2(w2), sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x4a7484aaul), Inc3_avx2(&w5, sigma1_avx2(w3), sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x5cb0a9dcul), Inc4_avx2(&w6, sigma1_avx2(w4), K(0x100ul), sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x76f988daul), Inc4_avx2(&w7, sigma1_avx2(w5), w0, K(0x11002000ul))));
        w8 = Add3(K(0x80000000ul), sigma1_avx2(w6), w1);
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x983e5152ul), w8));
        w9 = Add(sigma1_avx2(w7), w2);
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0xa831c66dul), w9));
        w10 = Add(sigma1_avx2(w8), w3);
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0xb00327c8ul), w10));
        w11 = Add(sigma1_avx2(w9), w4);
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0xbf597fc7ul), w11));
        w12 = Add(sigma1_avx2(w10), w5);
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0xc6e00bf3ul), w12));
        w13 = Add(sigma1_avx2(w11), w6);
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0xd5a79147ul), w13));
        w14 = Add3(sigma1_avx2(w12), w7, K(0x400022ul));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x06ca6351ul), w14));
        w15 = Add4(K(0x100ul), sigma1_avx2(w13), w8, sigma0_avx2(w0));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x14292967ul), w15));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x27b70a85ul), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x2e1b2138ul), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x4d2c6dfcul), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x53380d13ul), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x650a7354ul), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x766a0abbul), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x81c2c92eul), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x92722c85ul), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0xa2bfe8a1ul), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0xa81a664bul), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0xc24b8b70ul), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0xc76c51a3ul), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0xd192e819ul), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0xd6990624ul), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0xf40e3585ul), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x106aa070ul), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x19a4c116ul), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x1e376c08ul), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x2748774cul), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x34b0bcb5ul), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x391c0cb3ul), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0x4ed8aa4aul), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add(K(0x5b9cca4ful), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add(K(0x682e6ff3ul), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add(K(0x748f82eeul), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add(K(0x78a5636ful), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add(K(0x84c87814ul), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add(K(0x8cc70208ul), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add(K(0x90befffaul), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add(K(0xa4506cebul), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add5(K(0xbef9a3f7ul), w14, sigma1_avx2(w12), w7, sigma0_avx2(w15)));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add5(K(0xc67178f2ul), w15, sigma1_avx2(w13), w8, sigma0_avx2(w0)));

        /* Output */
        Write8_avx2(&out->u8[0], Add(a, K(0x6a09e667ul)));
        Write8_avx2(&out->u8[4], Add(b, K(0xbb67ae85ul)));
        Write8_avx2(&out->u8[8], Add(c, K(0x3c6ef372ul)));
        Write8_avx2(&out->u8[12], Add(d, K(0xa54ff53aul)));
        Write8_avx2(&out->u8[16], Add(e, K(0x510e527ful)));
        Write8_avx2(&out->u8[20], Add(f, K(0x9b05688cul)));
        Write8_avx2(&out->u8[24], Add(g, K(0x1f83d9abul)));
        Write8_avx2(&out->u8[28], Add(h, K(0x5be0cd19ul)));
}

#else
/* -Wempty-translation-unit
 * ISO C requires a translation unit to contain at least one declaration
 */
typedef int make_iso_compilers_happy;
#endif

/* End of File
 */
