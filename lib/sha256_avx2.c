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

#define K_avx2(x) _mm256_set1_epi32(x)

#define Add_avx2(x, y) _mm256_add_epi32((x), (y))
#define Add3_avx2(x, y, z) Add_avx2(Add_avx2((x), (y)), (z))
#define Add4_avx2(x, y, z, w) Add_avx2(Add_avx2((x), (y)), Add_avx2((z), (w)))
#define Add5_avx2(x, y, z, w, v) Add_avx2(Add3_avx2((x), (y), (z)), Add_avx2((w), (v)))
static inline __attribute__((always_inline)) __m256i Inc_avx2(__m256i *x, __m256i y) { *x = Add_avx2(*x, y); return *x; }
static inline __attribute__((always_inline)) __m256i Inc3_avx2(__m256i *x, __m256i y, __m256i z) { *x = Add3_avx2(*x, y, z); return *x; }
static inline __attribute__((always_inline)) __m256i Inc4_avx2(__m256i *x, __m256i y, __m256i z, __m256i w) { *x = Add4_avx2(*x, y, z, w); return *x; }
#define Xor_avx2(x, y) _mm256_xor_si256((x), (y))
#define Xor3_avx2(x, y, z) Xor_avx2(Xor_avx2((x), (y)), (z))
#define Or_avx2(x, y) _mm256_or_si256((x), (y))
#define And_avx2(x, y) _mm256_and_si256((x), (y))
#define ShR_avx2(x, n) _mm256_srli_epi32((x), (n))
#define ShL_avx2(x, n) _mm256_slli_epi32((x), (n))

static inline __attribute__((always_inline)) __m256i Ch_avx2(__m256i x, __m256i y, __m256i z) { return Xor_avx2(z, And_avx2(x, Xor_avx2(y, z))); }
static inline __attribute__((always_inline)) __m256i Maj_avx2(__m256i x, __m256i y, __m256i z) { return Or_avx2(And_avx2(x, y), And_avx2(z, Or_avx2(x, y))); }
static inline __attribute__((always_inline)) __m256i Sigma0_avx2(__m256i x) { return Xor3_avx2(Or_avx2(ShR_avx2(x, 2), ShL_avx2(x, 30)), Or_avx2(ShR_avx2(x, 13), ShL_avx2(x, 19)), Or_avx2(ShR_avx2(x, 22), ShL_avx2(x, 10))); }
static inline __attribute__((always_inline)) __m256i Sigma1_avx2(__m256i x) { return Xor3_avx2(Or_avx2(ShR_avx2(x, 6), ShL_avx2(x, 26)), Or_avx2(ShR_avx2(x, 11), ShL_avx2(x, 21)), Or_avx2(ShR_avx2(x, 25), ShL_avx2(x, 7))); }
static inline __attribute__((always_inline)) __m256i sigma0_avx2(__m256i x) { return Xor3_avx2(Or_avx2(ShR_avx2(x, 7), ShL_avx2(x, 25)), Or_avx2(ShR_avx2(x, 18), ShL_avx2(x, 14)), ShR_avx2(x, 3)); }
static inline __attribute__((always_inline)) __m256i sigma1_avx2(__m256i x) { return Xor3_avx2(Or_avx2(ShR_avx2(x, 17), ShL_avx2(x, 15)), Or_avx2(ShR_avx2(x, 19), ShL_avx2(x, 13)), ShR_avx2(x, 10)); }

/** One round of SHA-256. */
static inline __attribute__((always_inline)) void Round_avx2(__m256i a, __m256i b, __m256i c, __m256i *d, __m256i e, __m256i f, __m256i g, __m256i *h, __m256i k)
{
        __m256i t1 = Add4_avx2(*h, Sigma1_avx2(e), Ch_avx2(e, f, g), k);
        __m256i t2 = Add_avx2(Sigma0_avx2(a), Maj_avx2(a, b, c));
        *d = Add_avx2(*d, t1);
        *h = Add_avx2(t1, t2);
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
                        202182159, 134810123, 67438087, 66051,
                        202182159, 134810123, 67438087, 66051));
}

static inline __attribute__((always_inline)) void Write8_avx2(unsigned char *out, __m256i v)
{
        v = _mm256_shuffle_epi8(v, _mm256_set_epi32(
                202182159, 134810123, 67438087, 66051,
                202182159, 134810123, 67438087, 66051));
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
        __m256i a = K_avx2(*(int32_t*)&s[0]);
        __m256i b = K_avx2(*(int32_t*)&s[1]);
        __m256i c = K_avx2(*(int32_t*)&s[2]);
        __m256i d = K_avx2(*(int32_t*)&s[3]);
        __m256i e = K_avx2(*(int32_t*)&s[4]);
        __m256i f = K_avx2(*(int32_t*)&s[5]);
        __m256i g = K_avx2(*(int32_t*)&s[6]);
        __m256i h = K_avx2(*(int32_t*)&s[7]);

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

        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(1116352408), w0));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(1899447441), w1));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-1245643825), w2));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-373957723), w3));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(961987163), w4));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1508970993), w5));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-1841331548), w6));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(-1424204075), w7));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-670586216), w8));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(310598401), w9));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(607225278), w10));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(1426881987), w11));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(1925078388), w12));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-2132889090), w13));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-1680079193), w14));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(-1046744716), w15));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-459576895), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(-272742522), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(264347078), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(604807628), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(770255983), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1249150122), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(1555081692), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(1996064986), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-1740746414), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(-1473132947), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-1341970488), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-1084653625), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(-958395405), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-710438585), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(113926993), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(338241895), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(666307205), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(773529912), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(1294757372), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(1396182291), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(1695183700), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1986661051), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-2117940946), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(-1838011259), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-1564481375), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(-1474664885), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-1035236496), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-949202525), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(-778901479), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-694614492), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-200395387), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(275423344), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(430227734), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(506948616), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(659060556), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(883997877), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(958139571), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1322822218), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(1537002063), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(1747873779), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(1955562222), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(2024104815), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-2067236844), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-1933114872), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(-1866530822), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-1538233109), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-1090935817), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(-965641998), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));

        /* Output */
        Write8_avx2(&out->u8[0], Add_avx2(a, K_avx2(s[0])));
        Write8_avx2(&out->u8[4], Add_avx2(b, K_avx2(s[1])));
        Write8_avx2(&out->u8[8], Add_avx2(c, K_avx2(s[2])));
        Write8_avx2(&out->u8[12], Add_avx2(d, K_avx2(s[3])));
        Write8_avx2(&out->u8[16], Add_avx2(e, K_avx2(s[4])));
        Write8_avx2(&out->u8[20], Add_avx2(f, K_avx2(s[5])));
        Write8_avx2(&out->u8[24], Add_avx2(g, K_avx2(s[6])));
        Write8_avx2(&out->u8[28], Add_avx2(h, K_avx2(s[7])));
}

void transform_sha256d64_avx2_8way(struct sha256 out[8], const struct sha256 in[16])
{
        /* Transform 1 */
        __m256i a = K_avx2(1779033703);
        __m256i b = K_avx2(-1150833019);
        __m256i c = K_avx2(1013904242);
        __m256i d = K_avx2(-1521486534);
        __m256i e = K_avx2(1359893119);
        __m256i f = K_avx2(-1694144372);
        __m256i g = K_avx2(528734635);
        __m256i h = K_avx2(1541459225);

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

        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(1116352408), w0));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(1899447441), w1));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-1245643825), w2));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-373957723), w3));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(961987163), w4));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1508970993), w5));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-1841331548), w6));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(-1424204075), w7));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-670586216), w8));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(310598401), w9));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(607225278), w10));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(1426881987), w11));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(1925078388), w12));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-2132889090), w13));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-1680079193), w14));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(-1046744716), w15));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-459576895), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(-272742522), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(264347078), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(604807628), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(770255983), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1249150122), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(1555081692), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(1996064986), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-1740746414), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(-1473132947), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-1341970488), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-1084653625), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(-958395405), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-710438585), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(113926993), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(338241895), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(666307205), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(773529912), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(1294757372), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(1396182291), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(1695183700), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1986661051), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-2117940946), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(-1838011259), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-1564481375), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(-1474664885), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-1035236496), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-949202525), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(-778901479), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-694614492), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-200395387), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(275423344), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(430227734), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(506948616), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(659060556), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(883997877), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(958139571), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1322822218), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(1537002063), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(1747873779), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(1955562222), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(2024104815), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-2067236844), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-1933114872), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(-1866530822), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-1538233109), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-1090935817), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(-965641998), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));

        t0 = a = Add_avx2(a, K_avx2(1779033703));
        t1 = b = Add_avx2(b, K_avx2(-1150833019));
        t2 = c = Add_avx2(c, K_avx2(1013904242));
        t3 = d = Add_avx2(d, K_avx2(-1521486534));
        t4 = e = Add_avx2(e, K_avx2(1359893119));
        t5 = f = Add_avx2(f, K_avx2(-1694144372));
        t6 = g = Add_avx2(g, K_avx2(528734635));
        t7 = h = Add_avx2(h, K_avx2(1541459225));

        /* Transform 2 */
        Round_avx2(a, b, c, &d, e, f, g, &h, K_avx2(-1031131240));
        Round_avx2(h, a, b, &c, d, e, f, &g, K_avx2(1899447441));
        Round_avx2(g, h, a, &b, c, d, e, &f, K_avx2(-1245643825));
        Round_avx2(f, g, h, &a, b, c, d, &e, K_avx2(-373957723));
        Round_avx2(e, f, g, &h, a, b, c, &d, K_avx2(961987163));
        Round_avx2(d, e, f, &g, h, a, b, &c, K_avx2(1508970993));
        Round_avx2(c, d, e, &f, g, h, a, &b, K_avx2(-1841331548));
        Round_avx2(b, c, d, &e, f, g, h, &a, K_avx2(-1424204075));
        Round_avx2(a, b, c, &d, e, f, g, &h, K_avx2(-670586216));
        Round_avx2(h, a, b, &c, d, e, f, &g, K_avx2(310598401));
        Round_avx2(g, h, a, &b, c, d, e, &f, K_avx2(607225278));
        Round_avx2(f, g, h, &a, b, c, d, &e, K_avx2(1426881987));
        Round_avx2(e, f, g, &h, a, b, c, &d, K_avx2(1925078388));
        Round_avx2(d, e, f, &g, h, a, b, &c, K_avx2(-2132889090));
        Round_avx2(c, d, e, &f, g, h, a, &b, K_avx2(-1680079193));
        Round_avx2(b, c, d, &e, f, g, h, &a, K_avx2(-1046744204));
        Round_avx2(a, b, c, &d, e, f, g, &h, K_avx2(1687906753));
        Round_avx2(h, a, b, &c, d, e, f, &g, K_avx2(-251771002));
        Round_avx2(g, h, a, &b, c, d, e, &f, K_avx2(266464710));
        Round_avx2(f, g, h, &a, b, c, d, &e, K_avx2(604828244));
        Round_avx2(e, f, g, &h, a, b, c, &d, K_avx2(1340683375));
        Round_avx2(d, e, f, &g, h, a, b, &c, K_avx2(1825146046));
        Round_avx2(c, d, e, &f, g, h, a, &b, K_avx2(1639530782));
        Round_avx2(b, c, d, &e, f, g, h, &a, K_avx2(385452282));
        Round_avx2(a, b, c, &d, e, f, g, &h, K_avx2(-221884078));
        Round_avx2(h, a, b, &c, d, e, f, &g, K_avx2(-1467065747));
        Round_avx2(g, h, a, &b, c, d, e, &f, K_avx2(-1340474267));
        Round_avx2(f, g, h, &a, b, c, d, &e, K_avx2(-1176920377));
        Round_avx2(e, f, g, &h, a, b, c, &d, K_avx2(-1710083645));
        Round_avx2(d, e, f, &g, h, a, b, &c, K_avx2(-418452832));
        Round_avx2(c, d, e, &f, g, h, a, &b, K_avx2(-38722773));
        Round_avx2(b, c, d, &e, f, g, h, &a, K_avx2(-952811856));
        Round_avx2(a, b, c, &d, e, f, g, &h, K_avx2(812235477));
        Round_avx2(h, a, b, &c, d, e, f, &g, K_avx2(-879268513));
        Round_avx2(g, h, a, &b, c, d, e, &f, K_avx2(1510936975));
        Round_avx2(f, g, h, &a, b, c, d, &e, K_avx2(-601952515));
        Round_avx2(e, f, g, &h, a, b, c, &d, K_avx2(171292297));
        Round_avx2(d, e, f, &g, h, a, b, &c, K_avx2(-569673212));
        Round_avx2(c, d, e, &f, g, h, a, &b, K_avx2(1492437661));
        Round_avx2(b, c, d, &e, f, g, h, &a, K_avx2(-513975530));
        Round_avx2(a, b, c, &d, e, f, g, &h, K_avx2(8339078));
        Round_avx2(h, a, b, &c, d, e, f, &g, K_avx2(923306368));
        Round_avx2(g, h, a, &b, c, d, e, &f, K_avx2(-1526207950));
        Round_avx2(f, g, h, &a, b, c, d, &e, K_avx2(1873515831));
        Round_avx2(e, f, g, &h, a, b, c, &d, K_avx2(390095120));
        Round_avx2(d, e, f, &g, h, a, b, &c, K_avx2(227333873));
        Round_avx2(c, d, e, &f, g, h, a, &b, K_avx2(-844481683));
        Round_avx2(b, c, d, &e, f, g, h, &a, K_avx2(-1061437897));
        Round_avx2(a, b, c, &d, e, f, g, &h, K_avx2(-2090779686));
        Round_avx2(h, a, b, &c, d, e, f, &g, K_avx2(-615996573));
        Round_avx2(g, h, a, &b, c, d, e, &f, K_avx2(184740145));
        Round_avx2(f, g, h, &a, b, c, d, &e, K_avx2(1875991719));
        Round_avx2(e, f, g, &h, a, b, c, &d, K_avx2(1377499850));
        Round_avx2(d, e, f, &g, h, a, b, &c, K_avx2(825459761));
        Round_avx2(c, d, e, &f, g, h, a, &b, K_avx2(1859394197));
        Round_avx2(b, c, d, &e, f, g, h, &a, K_avx2(1833138320));
        Round_avx2(a, b, c, &d, e, f, g, &h, K_avx2(-1013149198));
        Round_avx2(h, a, b, &c, d, e, f, &g, K_avx2(-1630753859));
        Round_avx2(g, h, a, &b, c, d, e, &f, K_avx2(-1245077274));
        Round_avx2(f, g, h, &a, b, c, d, &e, K_avx2(1395635772));
        Round_avx2(e, f, g, &h, a, b, c, &d, K_avx2(-758693434));
        Round_avx2(d, e, f, &g, h, a, b, &c, K_avx2(119766691));
        Round_avx2(c, d, e, &f, g, h, a, &b, K_avx2(-1533719704));
        Round_avx2(b, c, d, &e, f, g, h, &a, K_avx2(1276714358));

        w0 = Add_avx2(t0, a);
        w1 = Add_avx2(t1, b);
        w2 = Add_avx2(t2, c);
        w3 = Add_avx2(t3, d);
        w4 = Add_avx2(t4, e);
        w5 = Add_avx2(t5, f);
        w6 = Add_avx2(t6, g);
        w7 = Add_avx2(t7, h);

        /* Transform 3 */
        a = K_avx2(1779033703);
        b = K_avx2(-1150833019);
        c = K_avx2(1013904242);
        d = K_avx2(-1521486534);
        e = K_avx2(1359893119);
        f = K_avx2(-1694144372);
        g = K_avx2(528734635);
        h = K_avx2(1541459225);

        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(1116352408), w0));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(1899447441), w1));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-1245643825), w2));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-373957723), w3));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(961987163), w4));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1508970993), w5));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-1841331548), w6));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(-1424204075), w7));
        Round_avx2(a, b, c, &d, e, f, g, &h, K_avx2(1476897432));
        Round_avx2(h, a, b, &c, d, e, f, &g, K_avx2(310598401));
        Round_avx2(g, h, a, &b, c, d, e, &f, K_avx2(607225278));
        Round_avx2(f, g, h, &a, b, c, d, &e, K_avx2(1426881987));
        Round_avx2(e, f, g, &h, a, b, c, &d, K_avx2(1925078388));
        Round_avx2(d, e, f, &g, h, a, b, &c, K_avx2(-2132889090));
        Round_avx2(c, d, e, &f, g, h, a, &b, K_avx2(-1680079193));
        Round_avx2(b, c, d, &e, f, g, h, &a, K_avx2(-1046744460));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-459576895), Inc_avx2(&w0, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(-272742522), Inc3_avx2(&w1, K_avx2(10485760), sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(264347078), Inc3_avx2(&w2, sigma1_avx2(w0), sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(604807628), Inc3_avx2(&w3, sigma1_avx2(w1), sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(770255983), Inc3_avx2(&w4, sigma1_avx2(w2), sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1249150122), Inc3_avx2(&w5, sigma1_avx2(w3), sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(1555081692), Inc4_avx2(&w6, sigma1_avx2(w4), K_avx2(256), sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(1996064986), Inc4_avx2(&w7, sigma1_avx2(w5), w0, K_avx2(285220864))));
        w8 = Add3_avx2(K_avx2(-2147483648), sigma1_avx2(w6), w1);
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-1740746414), w8));
        w9 = Add_avx2(sigma1_avx2(w7), w2);
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(-1473132947), w9));
        w10 = Add_avx2(sigma1_avx2(w8), w3);
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-1341970488), w10));
        w11 = Add_avx2(sigma1_avx2(w9), w4);
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-1084653625), w11));
        w12 = Add_avx2(sigma1_avx2(w10), w5);
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(-958395405), w12));
        w13 = Add_avx2(sigma1_avx2(w11), w6);
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-710438585), w13));
        w14 = Add3_avx2(sigma1_avx2(w12), w7, K_avx2(4194338));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(113926993), w14));
        w15 = Add4_avx2(K_avx2(256), sigma1_avx2(w13), w8, sigma0_avx2(w0));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(338241895), w15));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(666307205), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(773529912), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(1294757372), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(1396182291), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(1695183700), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1986661051), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-2117940946), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(-1838011259), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(-1564481375), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(-1474664885), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-1035236496), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-949202525), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(-778901479), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-694614492), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(-200395387), Inc4_avx2(&w14, sigma1_avx2(w12), w7, sigma0_avx2(w15))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(275423344), Inc4_avx2(&w15, sigma1_avx2(w13), w8, sigma0_avx2(w0))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(430227734), Inc4_avx2(&w0, sigma1_avx2(w14), w9, sigma0_avx2(w1))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(506948616), Inc4_avx2(&w1, sigma1_avx2(w15), w10, sigma0_avx2(w2))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(659060556), Inc4_avx2(&w2, sigma1_avx2(w0), w11, sigma0_avx2(w3))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(883997877), Inc4_avx2(&w3, sigma1_avx2(w1), w12, sigma0_avx2(w4))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(958139571), Inc4_avx2(&w4, sigma1_avx2(w2), w13, sigma0_avx2(w5))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(1322822218), Inc4_avx2(&w5, sigma1_avx2(w3), w14, sigma0_avx2(w6))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add_avx2(K_avx2(1537002063), Inc4_avx2(&w6, sigma1_avx2(w4), w15, sigma0_avx2(w7))));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add_avx2(K_avx2(1747873779), Inc4_avx2(&w7, sigma1_avx2(w5), w0, sigma0_avx2(w8))));
        Round_avx2(a, b, c, &d, e, f, g, &h, Add_avx2(K_avx2(1955562222), Inc4_avx2(&w8, sigma1_avx2(w6), w1, sigma0_avx2(w9))));
        Round_avx2(h, a, b, &c, d, e, f, &g, Add_avx2(K_avx2(2024104815), Inc4_avx2(&w9, sigma1_avx2(w7), w2, sigma0_avx2(w10))));
        Round_avx2(g, h, a, &b, c, d, e, &f, Add_avx2(K_avx2(-2067236844), Inc4_avx2(&w10, sigma1_avx2(w8), w3, sigma0_avx2(w11))));
        Round_avx2(f, g, h, &a, b, c, d, &e, Add_avx2(K_avx2(-1933114872), Inc4_avx2(&w11, sigma1_avx2(w9), w4, sigma0_avx2(w12))));
        Round_avx2(e, f, g, &h, a, b, c, &d, Add_avx2(K_avx2(-1866530822), Inc4_avx2(&w12, sigma1_avx2(w10), w5, sigma0_avx2(w13))));
        Round_avx2(d, e, f, &g, h, a, b, &c, Add_avx2(K_avx2(-1538233109), Inc4_avx2(&w13, sigma1_avx2(w11), w6, sigma0_avx2(w14))));
        Round_avx2(c, d, e, &f, g, h, a, &b, Add5_avx2(K_avx2(-1090935817), w14, sigma1_avx2(w12), w7, sigma0_avx2(w15)));
        Round_avx2(b, c, d, &e, f, g, h, &a, Add5_avx2(K_avx2(-965641998), w15, sigma1_avx2(w13), w8, sigma0_avx2(w0)));

        /* Output */
        Write8_avx2(&out->u8[0], Add_avx2(a, K_avx2(1779033703)));
        Write8_avx2(&out->u8[4], Add_avx2(b, K_avx2(-1150833019)));
        Write8_avx2(&out->u8[8], Add_avx2(c, K_avx2(1013904242)));
        Write8_avx2(&out->u8[12], Add_avx2(d, K_avx2(-1521486534)));
        Write8_avx2(&out->u8[16], Add_avx2(e, K_avx2(1359893119)));
        Write8_avx2(&out->u8[20], Add_avx2(f, K_avx2(-1694144372)));
        Write8_avx2(&out->u8[24], Add_avx2(g, K_avx2(528734635)));
        Write8_avx2(&out->u8[28], Add_avx2(h, K_avx2(1541459225)));
}

#else
/* -Wempty-translation-unit
 * ISO C requires a translation unit to contain at least one declaration
 */
typedef int avx2_make_iso_compilers_happy;
#endif

/* End of File
 */
