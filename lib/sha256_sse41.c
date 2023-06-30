/* Copyright (c) 2018-2019 The Bitcoin Core developers
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

#define K_sse41(x) _mm_set1_epi32(x)

#define Add_sse41(x, y) _mm_add_epi32((x), (y))
#define Add3_sse41(x, y, z) Add_sse41(Add_sse41((x), (y)), (z))
#define Add4_sse41(x, y, z, w) Add_sse41(Add_sse41((x), (y)), Add_sse41((z), (w)))
#define Add5_sse41(x, y, z, w, v) Add_sse41(Add3_sse41((x), (y), (z)), Add_sse41((w), (v)))
static inline __attribute__((always_inline)) __m128i Inc_sse41(__m128i *x, __m128i y) { *x = Add_sse41(*x, y); return *x; }
static inline __attribute__((always_inline)) __m128i Inc3_sse41(__m128i *x, __m128i y, __m128i z) { *x = Add3_sse41(*x, y, z); return *x; }
static inline __attribute__((always_inline)) __m128i Inc4_sse41(__m128i *x, __m128i y, __m128i z, __m128i w) { *x = Add4_sse41(*x, y, z, w); return *x; }
#define Xor_sse41(x, y) _mm_xor_si128((x), (y))
#define Xor3_sse41(x, y, z) Xor_sse41(Xor_sse41((x), (y)), (z))
#define Or_sse41(x, y) _mm_or_si128((x), (y))
#define And_sse41(x, y) _mm_and_si128((x), (y))
#define ShR_sse41(x, n) _mm_srli_epi32((x), (n))
#define ShL_sse41(x, n) _mm_slli_epi32((x), (n))

static inline __attribute__((always_inline)) __m128i Ch_sse41(__m128i x, __m128i y, __m128i z) { return Xor_sse41(z, And_sse41(x, Xor_sse41(y, z))); }
static inline __attribute__((always_inline)) __m128i Maj_sse41(__m128i x, __m128i y, __m128i z) { return Or_sse41(And_sse41(x, y), And_sse41(z, Or_sse41(x, y))); }
static inline __attribute__((always_inline)) __m128i Sigma0_sse41(__m128i x) { return Xor3_sse41(Or_sse41(ShR_sse41(x, 2), ShL_sse41(x, 30)), Or_sse41(ShR_sse41(x, 13), ShL_sse41(x, 19)), Or_sse41(ShR_sse41(x, 22), ShL_sse41(x, 10))); }
static inline __attribute__((always_inline)) __m128i Sigma1_sse41(__m128i x) { return Xor3_sse41(Or_sse41(ShR_sse41(x, 6), ShL_sse41(x, 26)), Or_sse41(ShR_sse41(x, 11), ShL_sse41(x, 21)), Or_sse41(ShR_sse41(x, 25), ShL_sse41(x, 7))); }
static inline __attribute__((always_inline)) __m128i sigma0_sse41(__m128i x) { return Xor3_sse41(Or_sse41(ShR_sse41(x, 7), ShL_sse41(x, 25)), Or_sse41(ShR_sse41(x, 18), ShL_sse41(x, 14)), ShR_sse41(x, 3)); }
static inline __attribute__((always_inline)) __m128i sigma1_sse41(__m128i x) { return Xor3_sse41(Or_sse41(ShR_sse41(x, 17), ShL_sse41(x, 15)), Or_sse41(ShR_sse41(x, 19), ShL_sse41(x, 13)), ShR_sse41(x, 10)); }

/** One round of SHA-256. */
static inline __attribute__((always_inline)) void Round_sse41(__m128i a, __m128i b, __m128i c, __m128i *d, __m128i e, __m128i f, __m128i g, __m128i *h, __m128i k)
{
        __m128i t1 = Add4_sse41(*h, Sigma1_sse41(e), Ch_sse41(e, f, g), k);
        __m128i t2 = Add_sse41(Sigma0_sse41(a), Maj_sse41(a, b, c));
        *d = Add_sse41(*d, t1);
        *h = Add_sse41(t1, t2);
}

static inline __attribute__((always_inline)) __m128i Read4_sse41(const unsigned char* chunk) {
        return _mm_shuffle_epi8(
                _mm_set_epi32(
                        ReadLE32(chunk + 0),
                        ReadLE32(chunk + 64),
                        ReadLE32(chunk + 128),
                        ReadLE32(chunk + 192)),
                _mm_set_epi32(202182159, 134810123, 67438087, 66051));
}

static inline __attribute__((always_inline)) void Write4_sse41(unsigned char *out, __m128i v) {
        v = _mm_shuffle_epi8(v, _mm_set_epi32(202182159, 134810123, 67438087, 66051));
        WriteLE32(out + 0, _mm_extract_epi32(v, 3));
        WriteLE32(out + 32, _mm_extract_epi32(v, 2));
        WriteLE32(out + 64, _mm_extract_epi32(v, 1));
        WriteLE32(out + 96, _mm_extract_epi32(v, 0));
}

void transform_sha256multi_sse41_4way(struct sha256* out, const uint32_t* s, const unsigned char* in)
{
        /* Transform 1 */
        __m128i a = K_sse41(s[0]);
        __m128i b = K_sse41(s[1]);
        __m128i c = K_sse41(s[2]);
        __m128i d = K_sse41(s[3]);
        __m128i e = K_sse41(s[4]);
        __m128i f = K_sse41(s[5]);
        __m128i g = K_sse41(s[6]);
        __m128i h = K_sse41(s[7]);

        __m128i w0 = Read4_sse41(in + 0),
                w1 = Read4_sse41(in + 4),
                w2 = Read4_sse41(in + 8),
                w3 = Read4_sse41(in + 12),
                w4 = Read4_sse41(in + 16),
                w5 = Read4_sse41(in + 20),
                w6 = Read4_sse41(in + 24),
                w7 = Read4_sse41(in + 28),
                w8 = Read4_sse41(in + 32),
                w9 = Read4_sse41(in + 36),
                w10 = Read4_sse41(in + 40),
                w11 = Read4_sse41(in + 44),
                w12 = Read4_sse41(in + 48),
                w13 = Read4_sse41(in + 52),
                w14 = Read4_sse41(in + 56),
                w15 = Read4_sse41(in + 60);

        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(1116352408), w0));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(1899447441), w1));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-1245643825), w2));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-373957723), w3));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(961987163), w4));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1508970993), w5));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-1841331548), w6));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(-1424204075), w7));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-670586216), w8));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(310598401), w9));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(607225278), w10));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(1426881987), w11));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(1925078388), w12));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-2132889090), w13));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-1680079193), w14));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(-1046744716), w15));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-459576895), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(-272742522), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(264347078), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(604807628), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(770255983), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1249150122), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(1555081692), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(1996064986), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-1740746414), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(-1473132947), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-1341970488), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-1084653625), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(-958395405), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-710438585), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(113926993), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(338241895), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(666307205), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(773529912), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(1294757372), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(1396182291), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(1695183700), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1986661051), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-2117940946), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(-1838011259), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-1564481375), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(-1474664885), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-1035236496), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-949202525), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(-778901479), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-694614492), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-200395387), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(275423344), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(430227734), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(506948616), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(659060556), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(883997877), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(958139571), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1322822218), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(1537002063), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(1747873779), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(1955562222), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(2024104815), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-2067236844), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-1933114872), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(-1866530822), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-1538233109), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-1090935817), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(-965641998), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));

        /* Output */
        Write4_sse41(&out->u8[0], Add_sse41(a, K_sse41(s[0])));
        Write4_sse41(&out->u8[4], Add_sse41(b, K_sse41(s[1])));
        Write4_sse41(&out->u8[8], Add_sse41(c, K_sse41(s[2])));
        Write4_sse41(&out->u8[12], Add_sse41(d, K_sse41(s[3])));
        Write4_sse41(&out->u8[16], Add_sse41(e, K_sse41(s[4])));
        Write4_sse41(&out->u8[20], Add_sse41(f, K_sse41(s[5])));
        Write4_sse41(&out->u8[24], Add_sse41(g, K_sse41(s[6])));
        Write4_sse41(&out->u8[28], Add_sse41(h, K_sse41(s[7])));
}

void transform_sha256d64_sse41_4way(struct sha256 out[4], const struct sha256 in[8])
{
        /* Transform 1 */
        __m128i a = K_sse41(1779033703);
        __m128i b = K_sse41(-1150833019);
        __m128i c = K_sse41(1013904242);
        __m128i d = K_sse41(-1521486534);
        __m128i e = K_sse41(1359893119);
        __m128i f = K_sse41(-1694144372);
        __m128i g = K_sse41(528734635);
        __m128i h = K_sse41(1541459225);

        __m128i w0 = Read4_sse41(&in[0].u8[0]),
                w1 = Read4_sse41(&in[0].u8[4]),
                w2 = Read4_sse41(&in[0].u8[8]),
                w3 = Read4_sse41(&in[0].u8[12]),
                w4 = Read4_sse41(&in[0].u8[16]),
                w5 = Read4_sse41(&in[0].u8[20]),
                w6 = Read4_sse41(&in[0].u8[24]),
                w7 = Read4_sse41(&in[0].u8[28]),
                w8 = Read4_sse41(&in[1].u8[0]),
                w9 = Read4_sse41(&in[1].u8[4]),
                w10 = Read4_sse41(&in[1].u8[8]),
                w11 = Read4_sse41(&in[1].u8[12]),
                w12 = Read4_sse41(&in[1].u8[16]),
                w13 = Read4_sse41(&in[1].u8[20]),
                w14 = Read4_sse41(&in[1].u8[24]),
                w15 = Read4_sse41(&in[1].u8[28]);

        __m128i t0, t1, t2, t3, t4, t5, t6, t7;

        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(1116352408), w0));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(1899447441), w1));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-1245643825), w2));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-373957723), w3));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(961987163), w4));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1508970993), w5));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-1841331548), w6));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(-1424204075), w7));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-670586216), w8));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(310598401), w9));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(607225278), w10));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(1426881987), w11));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(1925078388), w12));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-2132889090), w13));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-1680079193), w14));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(-1046744716), w15));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-459576895), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(-272742522), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(264347078), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(604807628), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(770255983), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1249150122), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(1555081692), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(1996064986), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-1740746414), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(-1473132947), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-1341970488), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-1084653625), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(-958395405), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-710438585), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(113926993), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(338241895), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(666307205), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(773529912), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(1294757372), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(1396182291), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(1695183700), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1986661051), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-2117940946), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(-1838011259), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-1564481375), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(-1474664885), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-1035236496), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-949202525), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(-778901479), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-694614492), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-200395387), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(275423344), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(430227734), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(506948616), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(659060556), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(883997877), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(958139571), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1322822218), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(1537002063), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(1747873779), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(1955562222), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(2024104815), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-2067236844), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-1933114872), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(-1866530822), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-1538233109), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-1090935817), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(-965641998), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));

        t0 = a = Add_sse41(a, K_sse41(1779033703));
        t1 = b = Add_sse41(b, K_sse41(-1150833019));
        t2 = c = Add_sse41(c, K_sse41(1013904242));
        t3 = d = Add_sse41(d, K_sse41(-1521486534));
        t4 = e = Add_sse41(e, K_sse41(1359893119));
        t5 = f = Add_sse41(f, K_sse41(-1694144372));
        t6 = g = Add_sse41(g, K_sse41(528734635));
        t7 = h = Add_sse41(h, K_sse41(1541459225));

        /* Transform 2 */
        Round_sse41(a, b, c, &d, e, f, g, &h, K_sse41(-1031131240));
        Round_sse41(h, a, b, &c, d, e, f, &g, K_sse41(1899447441));
        Round_sse41(g, h, a, &b, c, d, e, &f, K_sse41(-1245643825));
        Round_sse41(f, g, h, &a, b, c, d, &e, K_sse41(-373957723));
        Round_sse41(e, f, g, &h, a, b, c, &d, K_sse41(961987163));
        Round_sse41(d, e, f, &g, h, a, b, &c, K_sse41(1508970993));
        Round_sse41(c, d, e, &f, g, h, a, &b, K_sse41(-1841331548));
        Round_sse41(b, c, d, &e, f, g, h, &a, K_sse41(-1424204075));
        Round_sse41(a, b, c, &d, e, f, g, &h, K_sse41(-670586216));
        Round_sse41(h, a, b, &c, d, e, f, &g, K_sse41(310598401));
        Round_sse41(g, h, a, &b, c, d, e, &f, K_sse41(607225278));
        Round_sse41(f, g, h, &a, b, c, d, &e, K_sse41(1426881987));
        Round_sse41(e, f, g, &h, a, b, c, &d, K_sse41(1925078388));
        Round_sse41(d, e, f, &g, h, a, b, &c, K_sse41(-2132889090));
        Round_sse41(c, d, e, &f, g, h, a, &b, K_sse41(-1680079193));
        Round_sse41(b, c, d, &e, f, g, h, &a, K_sse41(-1046744204));
        Round_sse41(a, b, c, &d, e, f, g, &h, K_sse41(1687906753));
        Round_sse41(h, a, b, &c, d, e, f, &g, K_sse41(-251771002));
        Round_sse41(g, h, a, &b, c, d, e, &f, K_sse41(266464710));
        Round_sse41(f, g, h, &a, b, c, d, &e, K_sse41(604828244));
        Round_sse41(e, f, g, &h, a, b, c, &d, K_sse41(1340683375));
        Round_sse41(d, e, f, &g, h, a, b, &c, K_sse41(1825146046));
        Round_sse41(c, d, e, &f, g, h, a, &b, K_sse41(1639530782));
        Round_sse41(b, c, d, &e, f, g, h, &a, K_sse41(385452282));
        Round_sse41(a, b, c, &d, e, f, g, &h, K_sse41(-221884078));
        Round_sse41(h, a, b, &c, d, e, f, &g, K_sse41(-1467065747));
        Round_sse41(g, h, a, &b, c, d, e, &f, K_sse41(-1340474267));
        Round_sse41(f, g, h, &a, b, c, d, &e, K_sse41(-1176920377));
        Round_sse41(e, f, g, &h, a, b, c, &d, K_sse41(-1710083645));
        Round_sse41(d, e, f, &g, h, a, b, &c, K_sse41(-418452832));
        Round_sse41(c, d, e, &f, g, h, a, &b, K_sse41(-38722773));
        Round_sse41(b, c, d, &e, f, g, h, &a, K_sse41(-952811856));
        Round_sse41(a, b, c, &d, e, f, g, &h, K_sse41(812235477));
        Round_sse41(h, a, b, &c, d, e, f, &g, K_sse41(-879268513));
        Round_sse41(g, h, a, &b, c, d, e, &f, K_sse41(1510936975));
        Round_sse41(f, g, h, &a, b, c, d, &e, K_sse41(-601952515));
        Round_sse41(e, f, g, &h, a, b, c, &d, K_sse41(171292297));
        Round_sse41(d, e, f, &g, h, a, b, &c, K_sse41(-569673212));
        Round_sse41(c, d, e, &f, g, h, a, &b, K_sse41(1492437661));
        Round_sse41(b, c, d, &e, f, g, h, &a, K_sse41(-513975530));
        Round_sse41(a, b, c, &d, e, f, g, &h, K_sse41(8339078));
        Round_sse41(h, a, b, &c, d, e, f, &g, K_sse41(923306368));
        Round_sse41(g, h, a, &b, c, d, e, &f, K_sse41(-1526207950));
        Round_sse41(f, g, h, &a, b, c, d, &e, K_sse41(1873515831));
        Round_sse41(e, f, g, &h, a, b, c, &d, K_sse41(390095120));
        Round_sse41(d, e, f, &g, h, a, b, &c, K_sse41(227333873));
        Round_sse41(c, d, e, &f, g, h, a, &b, K_sse41(-844481683));
        Round_sse41(b, c, d, &e, f, g, h, &a, K_sse41(-1061437897));
        Round_sse41(a, b, c, &d, e, f, g, &h, K_sse41(-2090779686));
        Round_sse41(h, a, b, &c, d, e, f, &g, K_sse41(-615996573));
        Round_sse41(g, h, a, &b, c, d, e, &f, K_sse41(184740145));
        Round_sse41(f, g, h, &a, b, c, d, &e, K_sse41(1875991719));
        Round_sse41(e, f, g, &h, a, b, c, &d, K_sse41(1377499850));
        Round_sse41(d, e, f, &g, h, a, b, &c, K_sse41(825459761));
        Round_sse41(c, d, e, &f, g, h, a, &b, K_sse41(1859394197));
        Round_sse41(b, c, d, &e, f, g, h, &a, K_sse41(1833138320));
        Round_sse41(a, b, c, &d, e, f, g, &h, K_sse41(-1013149198));
        Round_sse41(h, a, b, &c, d, e, f, &g, K_sse41(-1630753859));
        Round_sse41(g, h, a, &b, c, d, e, &f, K_sse41(-1245077274));
        Round_sse41(f, g, h, &a, b, c, d, &e, K_sse41(1395635772));
        Round_sse41(e, f, g, &h, a, b, c, &d, K_sse41(-758693434));
        Round_sse41(d, e, f, &g, h, a, b, &c, K_sse41(119766691));
        Round_sse41(c, d, e, &f, g, h, a, &b, K_sse41(-1533719704));
        Round_sse41(b, c, d, &e, f, g, h, &a, K_sse41(1276714358));

        w0 = Add_sse41(t0, a);
        w1 = Add_sse41(t1, b);
        w2 = Add_sse41(t2, c);
        w3 = Add_sse41(t3, d);
        w4 = Add_sse41(t4, e);
        w5 = Add_sse41(t5, f);
        w6 = Add_sse41(t6, g);
        w7 = Add_sse41(t7, h);

        /* Transform 3 */
        a = K_sse41(1779033703);
        b = K_sse41(-1150833019);
        c = K_sse41(1013904242);
        d = K_sse41(-1521486534);
        e = K_sse41(1359893119);
        f = K_sse41(-1694144372);
        g = K_sse41(528734635);
        h = K_sse41(1541459225);

        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(1116352408), w0));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(1899447441), w1));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-1245643825), w2));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-373957723), w3));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(961987163), w4));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1508970993), w5));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-1841331548), w6));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(-1424204075), w7));
        Round_sse41(a, b, c, &d, e, f, g, &h, K_sse41(1476897432));
        Round_sse41(h, a, b, &c, d, e, f, &g, K_sse41(310598401));
        Round_sse41(g, h, a, &b, c, d, e, &f, K_sse41(607225278));
        Round_sse41(f, g, h, &a, b, c, d, &e, K_sse41(1426881987));
        Round_sse41(e, f, g, &h, a, b, c, &d, K_sse41(1925078388));
        Round_sse41(d, e, f, &g, h, a, b, &c, K_sse41(-2132889090));
        Round_sse41(c, d, e, &f, g, h, a, &b, K_sse41(-1680079193));
        Round_sse41(b, c, d, &e, f, g, h, &a, K_sse41(-1046744460));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-459576895), Inc_sse41(&w0, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(-272742522), Inc3_sse41(&w1, K_sse41(10485760), sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(264347078), Inc3_sse41(&w2, sigma1_sse41(w0), sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(604807628), Inc3_sse41(&w3, sigma1_sse41(w1), sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(770255983), Inc3_sse41(&w4, sigma1_sse41(w2), sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1249150122), Inc3_sse41(&w5, sigma1_sse41(w3), sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(1555081692), Inc4_sse41(&w6, sigma1_sse41(w4), K_sse41(256), sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(1996064986), Inc4_sse41(&w7, sigma1_sse41(w5), w0, K_sse41(285220864))));
        w8 = Add3_sse41(K_sse41(-2147483648), sigma1_sse41(w6), w1);
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-1740746414), w8));
        w9 = Add_sse41(sigma1_sse41(w7), w2);
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(-1473132947), w9));
        w10 = Add_sse41(sigma1_sse41(w8), w3);
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-1341970488), w10));
        w11 = Add_sse41(sigma1_sse41(w9), w4);
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-1084653625), w11));
        w12 = Add_sse41(sigma1_sse41(w10), w5);
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(-958395405), w12));
        w13 = Add_sse41(sigma1_sse41(w11), w6);
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-710438585), w13));
        w14 = Add3_sse41(sigma1_sse41(w12), w7, K_sse41(4194338));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(113926993), w14));
        w15 = Add4_sse41(K_sse41(256), sigma1_sse41(w13), w8, sigma0_sse41(w0));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(338241895), w15));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(666307205), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(773529912), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(1294757372), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(1396182291), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(1695183700), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1986661051), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-2117940946), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(-1838011259), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(-1564481375), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(-1474664885), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-1035236496), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-949202525), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(-778901479), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-694614492), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(-200395387), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(275423344), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(430227734), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(506948616), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(659060556), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(883997877), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(958139571), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(1322822218), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add_sse41(K_sse41(1537002063), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add_sse41(K_sse41(1747873779), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add_sse41(K_sse41(1955562222), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add_sse41(K_sse41(2024104815), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add_sse41(K_sse41(-2067236844), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add_sse41(K_sse41(-1933114872), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add_sse41(K_sse41(-1866530822), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add_sse41(K_sse41(-1538233109), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add5_sse41(K_sse41(-1090935817), w14, sigma1_sse41(w12), w7, sigma0_sse41(w15)));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add5_sse41(K_sse41(-965641998), w15, sigma1_sse41(w13), w8, sigma0_sse41(w0)));

        /* Output */
        Write4_sse41(&out->u8[0], Add_sse41(a, K_sse41(1779033703)));
        Write4_sse41(&out->u8[4], Add_sse41(b, K_sse41(-1150833019)));
        Write4_sse41(&out->u8[8], Add_sse41(c, K_sse41(1013904242)));
        Write4_sse41(&out->u8[12], Add_sse41(d, K_sse41(-1521486534)));
        Write4_sse41(&out->u8[16], Add_sse41(e, K_sse41(1359893119)));
        Write4_sse41(&out->u8[20], Add_sse41(f, K_sse41(-1694144372)));
        Write4_sse41(&out->u8[24], Add_sse41(g, K_sse41(528734635)));
        Write4_sse41(&out->u8[28], Add_sse41(h, K_sse41(1541459225)));
}

#else
/* -Wempty-translation-unit
 * ISO C requires a translation unit to contain at least one declaration
 */
typedef int sse41_make_iso_compilers_happy;
#endif

/* End of File
 */
