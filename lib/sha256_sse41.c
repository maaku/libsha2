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

#define K(x) _mm_set1_epi32(x)

#define Add(x, y) _mm_add_epi32((x), (y))
#define Add3(x, y, z) Add(Add((x), (y)), (z))
#define Add4(x, y, z, w) Add(Add((x), (y)), Add((z), (w)))
#define Add5(x, y, z, w, v) Add(Add3((x), (y), (z)), Add((w), (v)))
static inline __attribute__((always_inline)) __m128i Inc_sse41(__m128i *x, __m128i y) { *x = Add(*x, y); return *x; }
static inline __attribute__((always_inline)) __m128i Inc3_sse41(__m128i *x, __m128i y, __m128i z) { *x = Add3(*x, y, z); return *x; }
static inline __attribute__((always_inline)) __m128i Inc4_sse41(__m128i *x, __m128i y, __m128i z, __m128i w) { *x = Add4(*x, y, z, w); return *x; }
#define Xor(x, y) _mm_xor_si128((x), (y))
#define Xor3(x, y, z) Xor(Xor((x), (y)), (z))
#define Or(x, y) _mm_or_si128((x), (y))
#define And(x, y) _mm_and_si128((x), (y))
#define ShR(x, n) _mm_srli_epi32((x), (n))
#define ShL(x, n) _mm_slli_epi32((x), (n))

static inline __attribute__((always_inline)) __m128i Ch_sse41(__m128i x, __m128i y, __m128i z) { return Xor(z, And(x, Xor(y, z))); }
static inline __attribute__((always_inline)) __m128i Maj_sse41(__m128i x, __m128i y, __m128i z) { return Or(And(x, y), And(z, Or(x, y))); }
static inline __attribute__((always_inline)) __m128i Sigma0_sse41(__m128i x) { return Xor3(Or(ShR(x, 2), ShL(x, 30)), Or(ShR(x, 13), ShL(x, 19)), Or(ShR(x, 22), ShL(x, 10))); }
static inline __attribute__((always_inline)) __m128i Sigma1_sse41(__m128i x) { return Xor3(Or(ShR(x, 6), ShL(x, 26)), Or(ShR(x, 11), ShL(x, 21)), Or(ShR(x, 25), ShL(x, 7))); }
static inline __attribute__((always_inline)) __m128i sigma0_sse41(__m128i x) { return Xor3(Or(ShR(x, 7), ShL(x, 25)), Or(ShR(x, 18), ShL(x, 14)), ShR(x, 3)); }
static inline __attribute__((always_inline)) __m128i sigma1_sse41(__m128i x) { return Xor3(Or(ShR(x, 17), ShL(x, 15)), Or(ShR(x, 19), ShL(x, 13)), ShR(x, 10)); }

/** One round of SHA-256. */
static inline __attribute__((always_inline)) void Round_sse41(__m128i a, __m128i b, __m128i c, __m128i *d, __m128i e, __m128i f, __m128i g, __m128i *h, __m128i k)
{
        __m128i t1 = Add4(*h, Sigma1_sse41(e), Ch_sse41(e, f, g), k);
        __m128i t2 = Add(Sigma0_sse41(a), Maj_sse41(a, b, c));
        *d = Add(*d, t1);
        *h = Add(t1, t2);
}

static inline __attribute__((always_inline)) __m128i Read4_sse41(const unsigned char* chunk) {
        return _mm_shuffle_epi8(
                _mm_set_epi32(
                        ReadLE32(chunk + 0),
                        ReadLE32(chunk + 64),
                        ReadLE32(chunk + 128),
                        ReadLE32(chunk + 192)),
                _mm_set_epi32(
                        0x0C0D0E0FUL,
                        0x08090A0BUL,
                        0x04050607UL,
                        0x00010203UL));
}

static inline __attribute__((always_inline)) void Write4_sse41(unsigned char *out, __m128i v) {
        v = _mm_shuffle_epi8(v, _mm_set_epi32(
                0x0C0D0E0FUL,
                0x08090A0BUL,
                0x04050607UL,
                0x00010203UL));
        WriteLE32(out + 0, _mm_extract_epi32(v, 3));
        WriteLE32(out + 32, _mm_extract_epi32(v, 2));
        WriteLE32(out + 64, _mm_extract_epi32(v, 1));
        WriteLE32(out + 96, _mm_extract_epi32(v, 0));
}

void transform_sha256multi_sse41_4way(struct sha256* out, const uint32_t* s, const unsigned char* in)
{
        /* Transform 1 */
        __m128i a = K(s[0]);
        __m128i b = K(s[1]);
        __m128i c = K(s[2]);
        __m128i d = K(s[3]);
        __m128i e = K(s[4]);
        __m128i f = K(s[5]);
        __m128i g = K(s[6]);
        __m128i h = K(s[7]);

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

        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x428a2f98ul), w0));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x71374491ul), w1));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0xb5c0fbcful), w2));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0xe9b5dba5ul), w3));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x3956c25bul), w4));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x59f111f1ul), w5));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x923f82a4ul), w6));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0xab1c5ed5ul), w7));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0xd807aa98ul), w8));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x12835b01ul), w9));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x243185beul), w10));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x550c7dc3ul), w11));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x72be5d74ul), w12));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x80deb1feul), w13));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x9bdc06a7ul), w14));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0xc19bf174ul), w15));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0xe49b69c1ul), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0xefbe4786ul), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x0fc19dc6ul), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x240ca1ccul), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x2de92c6ful), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x4a7484aaul), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x5cb0a9dcul), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x76f988daul), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x983e5152ul), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0xa831c66dul), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0xb00327c8ul), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0xbf597fc7ul), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0xc6e00bf3ul), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0xd5a79147ul), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x06ca6351ul), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x14292967ul), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x27b70a85ul), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x2e1b2138ul), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x4d2c6dfcul), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x53380d13ul), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x650a7354ul), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x766a0abbul), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x81c2c92eul), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x92722c85ul), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0xa2bfe8a1ul), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0xa81a664bul), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0xc24b8b70ul), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0xc76c51a3ul), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0xd192e819ul), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0xd6990624ul), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0xf40e3585ul), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x106aa070ul), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x19a4c116ul), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x1e376c08ul), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x2748774cul), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x34b0bcb5ul), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x391c0cb3ul), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x4ed8aa4aul), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x5b9cca4ful), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x682e6ff3ul), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x748f82eeul), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x78a5636ful), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x84c87814ul), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x8cc70208ul), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x90befffaul), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0xa4506cebul), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0xbef9a3f7ul), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0xc67178f2ul), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));

        /* Output */
        Write4_sse41(&out->u8[0], Add(a, K(s[0])));
        Write4_sse41(&out->u8[4], Add(b, K(s[1])));
        Write4_sse41(&out->u8[8], Add(c, K(s[2])));
        Write4_sse41(&out->u8[12], Add(d, K(s[3])));
        Write4_sse41(&out->u8[16], Add(e, K(s[4])));
        Write4_sse41(&out->u8[20], Add(f, K(s[5])));
        Write4_sse41(&out->u8[24], Add(g, K(s[6])));
        Write4_sse41(&out->u8[28], Add(h, K(s[7])));
}

void transform_sha256d64_sse41_4way(struct sha256 out[4], const struct sha256 in[8])
{
        /* Transform 1 */
        __m128i a = K(0x6a09e667ul);
        __m128i b = K(0xbb67ae85ul);
        __m128i c = K(0x3c6ef372ul);
        __m128i d = K(0xa54ff53aul);
        __m128i e = K(0x510e527ful);
        __m128i f = K(0x9b05688cul);
        __m128i g = K(0x1f83d9abul);
        __m128i h = K(0x5be0cd19ul);

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

        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x428a2f98ul), w0));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x71374491ul), w1));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0xb5c0fbcful), w2));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0xe9b5dba5ul), w3));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x3956c25bul), w4));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x59f111f1ul), w5));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x923f82a4ul), w6));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0xab1c5ed5ul), w7));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0xd807aa98ul), w8));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x12835b01ul), w9));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x243185beul), w10));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x550c7dc3ul), w11));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x72be5d74ul), w12));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x80deb1feul), w13));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x9bdc06a7ul), w14));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0xc19bf174ul), w15));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0xe49b69c1ul), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0xefbe4786ul), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x0fc19dc6ul), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x240ca1ccul), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x2de92c6ful), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x4a7484aaul), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x5cb0a9dcul), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x76f988daul), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x983e5152ul), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0xa831c66dul), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0xb00327c8ul), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0xbf597fc7ul), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0xc6e00bf3ul), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0xd5a79147ul), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x06ca6351ul), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x14292967ul), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x27b70a85ul), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x2e1b2138ul), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x4d2c6dfcul), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x53380d13ul), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x650a7354ul), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x766a0abbul), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x81c2c92eul), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x92722c85ul), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0xa2bfe8a1ul), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0xa81a664bul), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0xc24b8b70ul), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0xc76c51a3ul), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0xd192e819ul), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0xd6990624ul), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0xf40e3585ul), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x106aa070ul), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x19a4c116ul), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x1e376c08ul), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x2748774cul), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x34b0bcb5ul), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x391c0cb3ul), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x4ed8aa4aul), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x5b9cca4ful), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x682e6ff3ul), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x748f82eeul), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x78a5636ful), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x84c87814ul), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x8cc70208ul), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x90befffaul), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0xa4506cebul), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0xbef9a3f7ul), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0xc67178f2ul), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));

        t0 = a = Add(a, K(0x6a09e667ul));
        t1 = b = Add(b, K(0xbb67ae85ul));
        t2 = c = Add(c, K(0x3c6ef372ul));
        t3 = d = Add(d, K(0xa54ff53aul));
        t4 = e = Add(e, K(0x510e527ful));
        t5 = f = Add(f, K(0x9b05688cul));
        t6 = g = Add(g, K(0x1f83d9abul));
        t7 = h = Add(h, K(0x5be0cd19ul));

        /* Transform 2 */
        Round_sse41(a, b, c, &d, e, f, g, &h, K(0xc28a2f98ul));
        Round_sse41(h, a, b, &c, d, e, f, &g, K(0x71374491ul));
        Round_sse41(g, h, a, &b, c, d, e, &f, K(0xb5c0fbcful));
        Round_sse41(f, g, h, &a, b, c, d, &e, K(0xe9b5dba5ul));
        Round_sse41(e, f, g, &h, a, b, c, &d, K(0x3956c25bul));
        Round_sse41(d, e, f, &g, h, a, b, &c, K(0x59f111f1ul));
        Round_sse41(c, d, e, &f, g, h, a, &b, K(0x923f82a4ul));
        Round_sse41(b, c, d, &e, f, g, h, &a, K(0xab1c5ed5ul));
        Round_sse41(a, b, c, &d, e, f, g, &h, K(0xd807aa98ul));
        Round_sse41(h, a, b, &c, d, e, f, &g, K(0x12835b01ul));
        Round_sse41(g, h, a, &b, c, d, e, &f, K(0x243185beul));
        Round_sse41(f, g, h, &a, b, c, d, &e, K(0x550c7dc3ul));
        Round_sse41(e, f, g, &h, a, b, c, &d, K(0x72be5d74ul));
        Round_sse41(d, e, f, &g, h, a, b, &c, K(0x80deb1feul));
        Round_sse41(c, d, e, &f, g, h, a, &b, K(0x9bdc06a7ul));
        Round_sse41(b, c, d, &e, f, g, h, &a, K(0xc19bf374ul));
        Round_sse41(a, b, c, &d, e, f, g, &h, K(0x649b69c1ul));
        Round_sse41(h, a, b, &c, d, e, f, &g, K(0xf0fe4786ul));
        Round_sse41(g, h, a, &b, c, d, e, &f, K(0x0fe1edc6ul));
        Round_sse41(f, g, h, &a, b, c, d, &e, K(0x240cf254ul));
        Round_sse41(e, f, g, &h, a, b, c, &d, K(0x4fe9346ful));
        Round_sse41(d, e, f, &g, h, a, b, &c, K(0x6cc984beul));
        Round_sse41(c, d, e, &f, g, h, a, &b, K(0x61b9411eul));
        Round_sse41(b, c, d, &e, f, g, h, &a, K(0x16f988faul));
        Round_sse41(a, b, c, &d, e, f, g, &h, K(0xf2c65152ul));
        Round_sse41(h, a, b, &c, d, e, f, &g, K(0xa88e5a6dul));
        Round_sse41(g, h, a, &b, c, d, e, &f, K(0xb019fc65ul));
        Round_sse41(f, g, h, &a, b, c, d, &e, K(0xb9d99ec7ul));
        Round_sse41(e, f, g, &h, a, b, c, &d, K(0x9a1231c3ul));
        Round_sse41(d, e, f, &g, h, a, b, &c, K(0xe70eeaa0ul));
        Round_sse41(c, d, e, &f, g, h, a, &b, K(0xfdb1232bul));
        Round_sse41(b, c, d, &e, f, g, h, &a, K(0xc7353eb0ul));
        Round_sse41(a, b, c, &d, e, f, g, &h, K(0x3069bad5ul));
        Round_sse41(h, a, b, &c, d, e, f, &g, K(0xcb976d5ful));
        Round_sse41(g, h, a, &b, c, d, e, &f, K(0x5a0f118ful));
        Round_sse41(f, g, h, &a, b, c, d, &e, K(0xdc1eeefdul));
        Round_sse41(e, f, g, &h, a, b, c, &d, K(0x0a35b689ul));
        Round_sse41(d, e, f, &g, h, a, b, &c, K(0xde0b7a04ul));
        Round_sse41(c, d, e, &f, g, h, a, &b, K(0x58f4ca9dul));
        Round_sse41(b, c, d, &e, f, g, h, &a, K(0xe15d5b16ul));
        Round_sse41(a, b, c, &d, e, f, g, &h, K(0x007f3e86ul));
        Round_sse41(h, a, b, &c, d, e, f, &g, K(0x37088980ul));
        Round_sse41(g, h, a, &b, c, d, e, &f, K(0xa507ea32ul));
        Round_sse41(f, g, h, &a, b, c, d, &e, K(0x6fab9537ul));
        Round_sse41(e, f, g, &h, a, b, c, &d, K(0x17406110ul));
        Round_sse41(d, e, f, &g, h, a, b, &c, K(0x0d8cd6f1ul));
        Round_sse41(c, d, e, &f, g, h, a, &b, K(0xcdaa3b6dul));
        Round_sse41(b, c, d, &e, f, g, h, &a, K(0xc0bbbe37ul));
        Round_sse41(a, b, c, &d, e, f, g, &h, K(0x83613bdaul));
        Round_sse41(h, a, b, &c, d, e, f, &g, K(0xdb48a363ul));
        Round_sse41(g, h, a, &b, c, d, e, &f, K(0x0b02e931ul));
        Round_sse41(f, g, h, &a, b, c, d, &e, K(0x6fd15ca7ul));
        Round_sse41(e, f, g, &h, a, b, c, &d, K(0x521afacaul));
        Round_sse41(d, e, f, &g, h, a, b, &c, K(0x31338431ul));
        Round_sse41(c, d, e, &f, g, h, a, &b, K(0x6ed41a95ul));
        Round_sse41(b, c, d, &e, f, g, h, &a, K(0x6d437890ul));
        Round_sse41(a, b, c, &d, e, f, g, &h, K(0xc39c91f2ul));
        Round_sse41(h, a, b, &c, d, e, f, &g, K(0x9eccabbdul));
        Round_sse41(g, h, a, &b, c, d, e, &f, K(0xb5c9a0e6ul));
        Round_sse41(f, g, h, &a, b, c, d, &e, K(0x532fb63cul));
        Round_sse41(e, f, g, &h, a, b, c, &d, K(0xd2c741c6ul));
        Round_sse41(d, e, f, &g, h, a, b, &c, K(0x07237ea3ul));
        Round_sse41(c, d, e, &f, g, h, a, &b, K(0xa4954b68ul));
        Round_sse41(b, c, d, &e, f, g, h, &a, K(0x4c191d76ul));

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

        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x428a2f98ul), w0));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x71374491ul), w1));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0xb5c0fbcful), w2));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0xe9b5dba5ul), w3));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x3956c25bul), w4));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x59f111f1ul), w5));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x923f82a4ul), w6));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0xab1c5ed5ul), w7));
        Round_sse41(a, b, c, &d, e, f, g, &h, K(0x5807aa98ul));
        Round_sse41(h, a, b, &c, d, e, f, &g, K(0x12835b01ul));
        Round_sse41(g, h, a, &b, c, d, e, &f, K(0x243185beul));
        Round_sse41(f, g, h, &a, b, c, d, &e, K(0x550c7dc3ul));
        Round_sse41(e, f, g, &h, a, b, c, &d, K(0x72be5d74ul));
        Round_sse41(d, e, f, &g, h, a, b, &c, K(0x80deb1feul));
        Round_sse41(c, d, e, &f, g, h, a, &b, K(0x9bdc06a7ul));
        Round_sse41(b, c, d, &e, f, g, h, &a, K(0xc19bf274ul));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0xe49b69c1ul), Inc_sse41(&w0, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0xefbe4786ul), Inc3_sse41(&w1, K(0xa00000ul), sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x0fc19dc6ul), Inc3_sse41(&w2, sigma1_sse41(w0), sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x240ca1ccul), Inc3_sse41(&w3, sigma1_sse41(w1), sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x2de92c6ful), Inc3_sse41(&w4, sigma1_sse41(w2), sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x4a7484aaul), Inc3_sse41(&w5, sigma1_sse41(w3), sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x5cb0a9dcul), Inc4_sse41(&w6, sigma1_sse41(w4), K(0x100ul), sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x76f988daul), Inc4_sse41(&w7, sigma1_sse41(w5), w0, K(0x11002000ul))));
        w8 = Add3(K(0x80000000ul), sigma1_sse41(w6), w1);
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x983e5152ul), w8));
        w9 = Add(sigma1_sse41(w7), w2);
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0xa831c66dul), w9));
        w10 = Add(sigma1_sse41(w8), w3);
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0xb00327c8ul), w10));
        w11 = Add(sigma1_sse41(w9), w4);
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0xbf597fc7ul), w11));
        w12 = Add(sigma1_sse41(w10), w5);
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0xc6e00bf3ul), w12));
        w13 = Add(sigma1_sse41(w11), w6);
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0xd5a79147ul), w13));
        w14 = Add3(sigma1_sse41(w12), w7, K(0x400022ul));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x06ca6351ul), w14));
        w15 = Add4(K(0x100ul), sigma1_sse41(w13), w8, sigma0_sse41(w0));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x14292967ul), w15));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x27b70a85ul), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x2e1b2138ul), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x4d2c6dfcul), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x53380d13ul), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x650a7354ul), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x766a0abbul), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x81c2c92eul), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x92722c85ul), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0xa2bfe8a1ul), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0xa81a664bul), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0xc24b8b70ul), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0xc76c51a3ul), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0xd192e819ul), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0xd6990624ul), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0xf40e3585ul), Inc4_sse41(&w14, sigma1_sse41(w12), w7, sigma0_sse41(w15))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x106aa070ul), Inc4_sse41(&w15, sigma1_sse41(w13), w8, sigma0_sse41(w0))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x19a4c116ul), Inc4_sse41(&w0, sigma1_sse41(w14), w9, sigma0_sse41(w1))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x1e376c08ul), Inc4_sse41(&w1, sigma1_sse41(w15), w10, sigma0_sse41(w2))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x2748774cul), Inc4_sse41(&w2, sigma1_sse41(w0), w11, sigma0_sse41(w3))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x34b0bcb5ul), Inc4_sse41(&w3, sigma1_sse41(w1), w12, sigma0_sse41(w4))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x391c0cb3ul), Inc4_sse41(&w4, sigma1_sse41(w2), w13, sigma0_sse41(w5))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0x4ed8aa4aul), Inc4_sse41(&w5, sigma1_sse41(w3), w14, sigma0_sse41(w6))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add(K(0x5b9cca4ful), Inc4_sse41(&w6, sigma1_sse41(w4), w15, sigma0_sse41(w7))));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add(K(0x682e6ff3ul), Inc4_sse41(&w7, sigma1_sse41(w5), w0, sigma0_sse41(w8))));
        Round_sse41(a, b, c, &d, e, f, g, &h, Add(K(0x748f82eeul), Inc4_sse41(&w8, sigma1_sse41(w6), w1, sigma0_sse41(w9))));
        Round_sse41(h, a, b, &c, d, e, f, &g, Add(K(0x78a5636ful), Inc4_sse41(&w9, sigma1_sse41(w7), w2, sigma0_sse41(w10))));
        Round_sse41(g, h, a, &b, c, d, e, &f, Add(K(0x84c87814ul), Inc4_sse41(&w10, sigma1_sse41(w8), w3, sigma0_sse41(w11))));
        Round_sse41(f, g, h, &a, b, c, d, &e, Add(K(0x8cc70208ul), Inc4_sse41(&w11, sigma1_sse41(w9), w4, sigma0_sse41(w12))));
        Round_sse41(e, f, g, &h, a, b, c, &d, Add(K(0x90befffaul), Inc4_sse41(&w12, sigma1_sse41(w10), w5, sigma0_sse41(w13))));
        Round_sse41(d, e, f, &g, h, a, b, &c, Add(K(0xa4506cebul), Inc4_sse41(&w13, sigma1_sse41(w11), w6, sigma0_sse41(w14))));
        Round_sse41(c, d, e, &f, g, h, a, &b, Add5(K(0xbef9a3f7ul), w14, sigma1_sse41(w12), w7, sigma0_sse41(w15)));
        Round_sse41(b, c, d, &e, f, g, h, &a, Add5(K(0xc67178f2ul), w15, sigma1_sse41(w13), w8, sigma0_sse41(w0)));

        /* Output */
        Write4_sse41(&out->u8[0], Add(a, K(0x6a09e667ul)));
        Write4_sse41(&out->u8[4], Add(b, K(0xbb67ae85ul)));
        Write4_sse41(&out->u8[8], Add(c, K(0x3c6ef372ul)));
        Write4_sse41(&out->u8[12], Add(d, K(0xa54ff53aul)));
        Write4_sse41(&out->u8[16], Add(e, K(0x510e527ful)));
        Write4_sse41(&out->u8[20], Add(f, K(0x9b05688cul)));
        Write4_sse41(&out->u8[24], Add(g, K(0x1f83d9abul)));
        Write4_sse41(&out->u8[28], Add(h, K(0x5be0cd19ul)));
}

#undef K
#undef Add
#undef Add3
#undef Add4
#undef Add5
#undef Xor
#undef Xor3
#undef Or
#undef And
#undef ShR
#undef ShL

#else
/* -Wempty-translation-unit
 * ISO C requires a translation unit to contain at least one declaration
 */
typedef int make_iso_compilers_happy;
#endif

/* End of File
 */
