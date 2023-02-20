/* Copyright (c) 2014-2019 The Bitcoin Core developers
 * Copyright (c) 2022 Mark Friedenbach
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <sha2/sha256.h>
#include "common.h"

#include <assert.h>
#include <string.h>

#include "compat/cpuid.h"

#if defined(__linux__) && (defined(__arm__) || defined(__aarch64__))
#include <sys/auxv.h>
#include <asm/hwcap.h>
#endif

#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#include "sha256_internal.h"

static inline __attribute__((always_inline)) uint32_t Ch_noasm(uint32_t x, uint32_t y, uint32_t z) { return z ^ (x & (y ^ z)); }
static inline __attribute__((always_inline)) uint32_t Maj_noasm(uint32_t x, uint32_t y, uint32_t z) { return (x & y) | (z & (x | y)); }
static inline __attribute__((always_inline)) uint32_t Sigma0_noasm(uint32_t x) { return (x >> 2 | x << 30) ^ (x >> 13 | x << 19) ^ (x >> 22 | x << 10); }
static inline __attribute__((always_inline)) uint32_t Sigma1_noasm(uint32_t x) { return (x >> 6 | x << 26) ^ (x >> 11 | x << 21) ^ (x >> 25 | x << 7); }
static inline __attribute__((always_inline)) uint32_t sigma0_noasm(uint32_t x) { return (x >> 7 | x << 25) ^ (x >> 18 | x << 14) ^ (x >> 3); }
static inline __attribute__((always_inline)) uint32_t sigma1_noasm(uint32_t x) { return (x >> 17 | x << 15) ^ (x >> 19 | x << 13) ^ (x >> 10); }

/** One round of SHA-256. */
static inline __attribute__((always_inline)) void Round_noasm(uint32_t a, uint32_t b, uint32_t c, uint32_t* d, uint32_t e, uint32_t f, uint32_t g, uint32_t* h, uint32_t k)
{
        uint32_t t1 = *h + Sigma1_noasm(e) + Ch_noasm(e, f, g) + k;
        uint32_t t2 = Sigma0_noasm(a) + Maj_noasm(a, b, c);
        *d += t1;
        *h = t1 + t2;
}

/** Initialize SHA-256 state. */
static inline __attribute__((always_inline)) void Initialize(uint32_t* s)
{
        s[0] = 0x6a09e667ul;
        s[1] = 0xbb67ae85ul;
        s[2] = 0x3c6ef372ul;
        s[3] = 0xa54ff53aul;
        s[4] = 0x510e527ful;
        s[5] = 0x9b05688cul;
        s[6] = 0x1f83d9abul;
        s[7] = 0x5be0cd19ul;
}

/** Perform a number of SHA-256 transformations, processing 64-byte chunks. */
static void transform_noasm(uint32_t* s, const unsigned char* chunk, size_t blocks)
{
        while (blocks--) {
                uint32_t a = s[0], b = s[1], c = s[2], d = s[3], e = s[4], f = s[5], g = s[6], h = s[7];
                uint32_t w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15;

                Round_noasm(a, b, c, &d, e, f, g, &h, 0x428a2f98 + (w0 = ReadBE32(chunk + 0)));
                Round_noasm(h, a, b, &c, d, e, f, &g, 0x71374491 + (w1 = ReadBE32(chunk + 4)));
                Round_noasm(g, h, a, &b, c, d, e, &f, 0xb5c0fbcf + (w2 = ReadBE32(chunk + 8)));
                Round_noasm(f, g, h, &a, b, c, d, &e, 0xe9b5dba5 + (w3 = ReadBE32(chunk + 12)));
                Round_noasm(e, f, g, &h, a, b, c, &d, 0x3956c25b + (w4 = ReadBE32(chunk + 16)));
                Round_noasm(d, e, f, &g, h, a, b, &c, 0x59f111f1 + (w5 = ReadBE32(chunk + 20)));
                Round_noasm(c, d, e, &f, g, h, a, &b, 0x923f82a4 + (w6 = ReadBE32(chunk + 24)));
                Round_noasm(b, c, d, &e, f, g, h, &a, 0xab1c5ed5 + (w7 = ReadBE32(chunk + 28)));
                Round_noasm(a, b, c, &d, e, f, g, &h, 0xd807aa98 + (w8 = ReadBE32(chunk + 32)));
                Round_noasm(h, a, b, &c, d, e, f, &g, 0x12835b01 + (w9 = ReadBE32(chunk + 36)));
                Round_noasm(g, h, a, &b, c, d, e, &f, 0x243185be + (w10 = ReadBE32(chunk + 40)));
                Round_noasm(f, g, h, &a, b, c, d, &e, 0x550c7dc3 + (w11 = ReadBE32(chunk + 44)));
                Round_noasm(e, f, g, &h, a, b, c, &d, 0x72be5d74 + (w12 = ReadBE32(chunk + 48)));
                Round_noasm(d, e, f, &g, h, a, b, &c, 0x80deb1fe + (w13 = ReadBE32(chunk + 52)));
                Round_noasm(c, d, e, &f, g, h, a, &b, 0x9bdc06a7 + (w14 = ReadBE32(chunk + 56)));
                Round_noasm(b, c, d, &e, f, g, h, &a, 0xc19bf174 + (w15 = ReadBE32(chunk + 60)));

                Round_noasm(a, b, c, &d, e, f, g, &h, 0xe49b69c1 + (w0 += sigma1_noasm(w14) + w9 + sigma0_noasm(w1)));
                Round_noasm(h, a, b, &c, d, e, f, &g, 0xefbe4786 + (w1 += sigma1_noasm(w15) + w10 + sigma0_noasm(w2)));
                Round_noasm(g, h, a, &b, c, d, e, &f, 0x0fc19dc6 + (w2 += sigma1_noasm(w0) + w11 + sigma0_noasm(w3)));
                Round_noasm(f, g, h, &a, b, c, d, &e, 0x240ca1cc + (w3 += sigma1_noasm(w1) + w12 + sigma0_noasm(w4)));
                Round_noasm(e, f, g, &h, a, b, c, &d, 0x2de92c6f + (w4 += sigma1_noasm(w2) + w13 + sigma0_noasm(w5)));
                Round_noasm(d, e, f, &g, h, a, b, &c, 0x4a7484aa + (w5 += sigma1_noasm(w3) + w14 + sigma0_noasm(w6)));
                Round_noasm(c, d, e, &f, g, h, a, &b, 0x5cb0a9dc + (w6 += sigma1_noasm(w4) + w15 + sigma0_noasm(w7)));
                Round_noasm(b, c, d, &e, f, g, h, &a, 0x76f988da + (w7 += sigma1_noasm(w5) + w0 + sigma0_noasm(w8)));
                Round_noasm(a, b, c, &d, e, f, g, &h, 0x983e5152 + (w8 += sigma1_noasm(w6) + w1 + sigma0_noasm(w9)));
                Round_noasm(h, a, b, &c, d, e, f, &g, 0xa831c66d + (w9 += sigma1_noasm(w7) + w2 + sigma0_noasm(w10)));
                Round_noasm(g, h, a, &b, c, d, e, &f, 0xb00327c8 + (w10 += sigma1_noasm(w8) + w3 + sigma0_noasm(w11)));
                Round_noasm(f, g, h, &a, b, c, d, &e, 0xbf597fc7 + (w11 += sigma1_noasm(w9) + w4 + sigma0_noasm(w12)));
                Round_noasm(e, f, g, &h, a, b, c, &d, 0xc6e00bf3 + (w12 += sigma1_noasm(w10) + w5 + sigma0_noasm(w13)));
                Round_noasm(d, e, f, &g, h, a, b, &c, 0xd5a79147 + (w13 += sigma1_noasm(w11) + w6 + sigma0_noasm(w14)));
                Round_noasm(c, d, e, &f, g, h, a, &b, 0x06ca6351 + (w14 += sigma1_noasm(w12) + w7 + sigma0_noasm(w15)));
                Round_noasm(b, c, d, &e, f, g, h, &a, 0x14292967 + (w15 += sigma1_noasm(w13) + w8 + sigma0_noasm(w0)));

                Round_noasm(a, b, c, &d, e, f, g, &h, 0x27b70a85 + (w0 += sigma1_noasm(w14) + w9 + sigma0_noasm(w1)));
                Round_noasm(h, a, b, &c, d, e, f, &g, 0x2e1b2138 + (w1 += sigma1_noasm(w15) + w10 + sigma0_noasm(w2)));
                Round_noasm(g, h, a, &b, c, d, e, &f, 0x4d2c6dfc + (w2 += sigma1_noasm(w0) + w11 + sigma0_noasm(w3)));
                Round_noasm(f, g, h, &a, b, c, d, &e, 0x53380d13 + (w3 += sigma1_noasm(w1) + w12 + sigma0_noasm(w4)));
                Round_noasm(e, f, g, &h, a, b, c, &d, 0x650a7354 + (w4 += sigma1_noasm(w2) + w13 + sigma0_noasm(w5)));
                Round_noasm(d, e, f, &g, h, a, b, &c, 0x766a0abb + (w5 += sigma1_noasm(w3) + w14 + sigma0_noasm(w6)));
                Round_noasm(c, d, e, &f, g, h, a, &b, 0x81c2c92e + (w6 += sigma1_noasm(w4) + w15 + sigma0_noasm(w7)));
                Round_noasm(b, c, d, &e, f, g, h, &a, 0x92722c85 + (w7 += sigma1_noasm(w5) + w0 + sigma0_noasm(w8)));
                Round_noasm(a, b, c, &d, e, f, g, &h, 0xa2bfe8a1 + (w8 += sigma1_noasm(w6) + w1 + sigma0_noasm(w9)));
                Round_noasm(h, a, b, &c, d, e, f, &g, 0xa81a664b + (w9 += sigma1_noasm(w7) + w2 + sigma0_noasm(w10)));
                Round_noasm(g, h, a, &b, c, d, e, &f, 0xc24b8b70 + (w10 += sigma1_noasm(w8) + w3 + sigma0_noasm(w11)));
                Round_noasm(f, g, h, &a, b, c, d, &e, 0xc76c51a3 + (w11 += sigma1_noasm(w9) + w4 + sigma0_noasm(w12)));
                Round_noasm(e, f, g, &h, a, b, c, &d, 0xd192e819 + (w12 += sigma1_noasm(w10) + w5 + sigma0_noasm(w13)));
                Round_noasm(d, e, f, &g, h, a, b, &c, 0xd6990624 + (w13 += sigma1_noasm(w11) + w6 + sigma0_noasm(w14)));
                Round_noasm(c, d, e, &f, g, h, a, &b, 0xf40e3585 + (w14 += sigma1_noasm(w12) + w7 + sigma0_noasm(w15)));
                Round_noasm(b, c, d, &e, f, g, h, &a, 0x106aa070 + (w15 += sigma1_noasm(w13) + w8 + sigma0_noasm(w0)));

                Round_noasm(a, b, c, &d, e, f, g, &h, 0x19a4c116 + (w0 += sigma1_noasm(w14) + w9 + sigma0_noasm(w1)));
                Round_noasm(h, a, b, &c, d, e, f, &g, 0x1e376c08 + (w1 += sigma1_noasm(w15) + w10 + sigma0_noasm(w2)));
                Round_noasm(g, h, a, &b, c, d, e, &f, 0x2748774c + (w2 += sigma1_noasm(w0) + w11 + sigma0_noasm(w3)));
                Round_noasm(f, g, h, &a, b, c, d, &e, 0x34b0bcb5 + (w3 += sigma1_noasm(w1) + w12 + sigma0_noasm(w4)));
                Round_noasm(e, f, g, &h, a, b, c, &d, 0x391c0cb3 + (w4 += sigma1_noasm(w2) + w13 + sigma0_noasm(w5)));
                Round_noasm(d, e, f, &g, h, a, b, &c, 0x4ed8aa4a + (w5 += sigma1_noasm(w3) + w14 + sigma0_noasm(w6)));
                Round_noasm(c, d, e, &f, g, h, a, &b, 0x5b9cca4f + (w6 += sigma1_noasm(w4) + w15 + sigma0_noasm(w7)));
                Round_noasm(b, c, d, &e, f, g, h, &a, 0x682e6ff3 + (w7 += sigma1_noasm(w5) + w0 + sigma0_noasm(w8)));
                Round_noasm(a, b, c, &d, e, f, g, &h, 0x748f82ee + (w8 += sigma1_noasm(w6) + w1 + sigma0_noasm(w9)));
                Round_noasm(h, a, b, &c, d, e, f, &g, 0x78a5636f + (w9 += sigma1_noasm(w7) + w2 + sigma0_noasm(w10)));
                Round_noasm(g, h, a, &b, c, d, e, &f, 0x84c87814 + (w10 += sigma1_noasm(w8) + w3 + sigma0_noasm(w11)));
                Round_noasm(f, g, h, &a, b, c, d, &e, 0x8cc70208 + (w11 += sigma1_noasm(w9) + w4 + sigma0_noasm(w12)));
                Round_noasm(e, f, g, &h, a, b, c, &d, 0x90befffa + (w12 += sigma1_noasm(w10) + w5 + sigma0_noasm(w13)));
                Round_noasm(d, e, f, &g, h, a, b, &c, 0xa4506ceb + (w13 += sigma1_noasm(w11) + w6 + sigma0_noasm(w14)));
                Round_noasm(c, d, e, &f, g, h, a, &b, 0xbef9a3f7 + (w14 + sigma1_noasm(w12) + w7 + sigma0_noasm(w15)));
                Round_noasm(b, c, d, &e, f, g, h, &a, 0xc67178f2 + (w15 + sigma1_noasm(w13) + w8 + sigma0_noasm(w0)));

                s[0] += a;
                s[1] += b;
                s[2] += c;
                s[3] += d;
                s[4] += e;
                s[5] += f;
                s[6] += g;
                s[7] += h;
                chunk += 64;
        }
}

static void transform_d64_noasm(struct sha256 out[1], const struct sha256 in[2])
{
        /* Transform 1 */
        uint32_t a = 0x6a09e667ul;
        uint32_t b = 0xbb67ae85ul;
        uint32_t c = 0x3c6ef372ul;
        uint32_t d = 0xa54ff53aul;
        uint32_t e = 0x510e527ful;
        uint32_t f = 0x9b05688cul;
        uint32_t g = 0x1f83d9abul;
        uint32_t h = 0x5be0cd19ul;

        uint32_t w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15;

        uint32_t t0, t1, t2, t3, t4, t5, t6, t7;

        Round_noasm(a, b, c, &d, e, f, g, &h, 0x428a2f98ul + (w0 = ReadBE32(&in[0].u8[0])));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x71374491ul + (w1 = ReadBE32(&in[0].u8[4])));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0xb5c0fbcful + (w2 = ReadBE32(&in[0].u8[8])));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0xe9b5dba5ul + (w3 = ReadBE32(&in[0].u8[12])));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x3956c25bul + (w4 = ReadBE32(&in[0].u8[16])));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x59f111f1ul + (w5 = ReadBE32(&in[0].u8[20])));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x923f82a4ul + (w6 = ReadBE32(&in[0].u8[24])));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xab1c5ed5ul + (w7 = ReadBE32(&in[0].u8[28])));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0xd807aa98ul + (w8 = ReadBE32(&in[1].u8[0])));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x12835b01ul + (w9 = ReadBE32(&in[1].u8[4])));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x243185beul + (w10 = ReadBE32(&in[1].u8[8])));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x550c7dc3ul + (w11 = ReadBE32(&in[1].u8[12])));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x72be5d74ul + (w12 = ReadBE32(&in[1].u8[16])));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x80deb1feul + (w13 = ReadBE32(&in[1].u8[20])));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x9bdc06a7ul + (w14 = ReadBE32(&in[1].u8[24])));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xc19bf174ul + (w15 = ReadBE32(&in[1].u8[28])));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0xe49b69c1ul + (w0 += sigma1_noasm(w14) + w9 + sigma0_noasm(w1)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0xefbe4786ul + (w1 += sigma1_noasm(w15) + w10 + sigma0_noasm(w2)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x0fc19dc6ul + (w2 += sigma1_noasm(w0) + w11 + sigma0_noasm(w3)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x240ca1ccul + (w3 += sigma1_noasm(w1) + w12 + sigma0_noasm(w4)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x2de92c6ful + (w4 += sigma1_noasm(w2) + w13 + sigma0_noasm(w5)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x4a7484aaul + (w5 += sigma1_noasm(w3) + w14 + sigma0_noasm(w6)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x5cb0a9dcul + (w6 += sigma1_noasm(w4) + w15 + sigma0_noasm(w7)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x76f988daul + (w7 += sigma1_noasm(w5) + w0 + sigma0_noasm(w8)));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x983e5152ul + (w8 += sigma1_noasm(w6) + w1 + sigma0_noasm(w9)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0xa831c66dul + (w9 += sigma1_noasm(w7) + w2 + sigma0_noasm(w10)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0xb00327c8ul + (w10 += sigma1_noasm(w8) + w3 + sigma0_noasm(w11)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0xbf597fc7ul + (w11 += sigma1_noasm(w9) + w4 + sigma0_noasm(w12)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0xc6e00bf3ul + (w12 += sigma1_noasm(w10) + w5 + sigma0_noasm(w13)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0xd5a79147ul + (w13 += sigma1_noasm(w11) + w6 + sigma0_noasm(w14)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x06ca6351ul + (w14 += sigma1_noasm(w12) + w7 + sigma0_noasm(w15)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x14292967ul + (w15 += sigma1_noasm(w13) + w8 + sigma0_noasm(w0)));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x27b70a85ul + (w0 += sigma1_noasm(w14) + w9 + sigma0_noasm(w1)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x2e1b2138ul + (w1 += sigma1_noasm(w15) + w10 + sigma0_noasm(w2)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x4d2c6dfcul + (w2 += sigma1_noasm(w0) + w11 + sigma0_noasm(w3)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x53380d13ul + (w3 += sigma1_noasm(w1) + w12 + sigma0_noasm(w4)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x650a7354ul + (w4 += sigma1_noasm(w2) + w13 + sigma0_noasm(w5)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x766a0abbul + (w5 += sigma1_noasm(w3) + w14 + sigma0_noasm(w6)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x81c2c92eul + (w6 += sigma1_noasm(w4) + w15 + sigma0_noasm(w7)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x92722c85ul + (w7 += sigma1_noasm(w5) + w0 + sigma0_noasm(w8)));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0xa2bfe8a1ul + (w8 += sigma1_noasm(w6) + w1 + sigma0_noasm(w9)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0xa81a664bul + (w9 += sigma1_noasm(w7) + w2 + sigma0_noasm(w10)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0xc24b8b70ul + (w10 += sigma1_noasm(w8) + w3 + sigma0_noasm(w11)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0xc76c51a3ul + (w11 += sigma1_noasm(w9) + w4 + sigma0_noasm(w12)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0xd192e819ul + (w12 += sigma1_noasm(w10) + w5 + sigma0_noasm(w13)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0xd6990624ul + (w13 += sigma1_noasm(w11) + w6 + sigma0_noasm(w14)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0xf40e3585ul + (w14 += sigma1_noasm(w12) + w7 + sigma0_noasm(w15)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x106aa070ul + (w15 += sigma1_noasm(w13) + w8 + sigma0_noasm(w0)));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x19a4c116ul + (w0 += sigma1_noasm(w14) + w9 + sigma0_noasm(w1)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x1e376c08ul + (w1 += sigma1_noasm(w15) + w10 + sigma0_noasm(w2)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x2748774cul + (w2 += sigma1_noasm(w0) + w11 + sigma0_noasm(w3)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x34b0bcb5ul + (w3 += sigma1_noasm(w1) + w12 + sigma0_noasm(w4)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x391c0cb3ul + (w4 += sigma1_noasm(w2) + w13 + sigma0_noasm(w5)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x4ed8aa4aul + (w5 += sigma1_noasm(w3) + w14 + sigma0_noasm(w6)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x5b9cca4ful + (w6 += sigma1_noasm(w4) + w15 + sigma0_noasm(w7)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x682e6ff3ul + (w7 += sigma1_noasm(w5) + w0 + sigma0_noasm(w8)));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x748f82eeul + (w8 += sigma1_noasm(w6) + w1 + sigma0_noasm(w9)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x78a5636ful + (w9 += sigma1_noasm(w7) + w2 + sigma0_noasm(w10)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x84c87814ul + (w10 += sigma1_noasm(w8) + w3 + sigma0_noasm(w11)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x8cc70208ul + (w11 += sigma1_noasm(w9) + w4 + sigma0_noasm(w12)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x90befffaul + (w12 += sigma1_noasm(w10) + w5 + sigma0_noasm(w13)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0xa4506cebul + (w13 += sigma1_noasm(w11) + w6 + sigma0_noasm(w14)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0xbef9a3f7ul + (w14 + sigma1_noasm(w12) + w7 + sigma0_noasm(w15)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xc67178f2ul + (w15 + sigma1_noasm(w13) + w8 + sigma0_noasm(w0)));

        a += 0x6a09e667ul; t0 = a;
        b += 0xbb67ae85ul; t1 = b;
        c += 0x3c6ef372ul; t2 = c;
        d += 0xa54ff53aul; t3 = d;
        e += 0x510e527ful; t4 = e;
        f += 0x9b05688cul; t5 = f;
        g += 0x1f83d9abul; t6 = g;
        h += 0x5be0cd19ul; t7 = h;

        /* Transform 2 */
        Round_noasm(a, b, c, &d, e, f, g, &h, 0xc28a2f98ul);
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x71374491ul);
        Round_noasm(g, h, a, &b, c, d, e, &f, 0xb5c0fbcful);
        Round_noasm(f, g, h, &a, b, c, d, &e, 0xe9b5dba5ul);
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x3956c25bul);
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x59f111f1ul);
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x923f82a4ul);
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xab1c5ed5ul);
        Round_noasm(a, b, c, &d, e, f, g, &h, 0xd807aa98ul);
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x12835b01ul);
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x243185beul);
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x550c7dc3ul);
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x72be5d74ul);
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x80deb1feul);
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x9bdc06a7ul);
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xc19bf374ul);
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x649b69c1ul);
        Round_noasm(h, a, b, &c, d, e, f, &g, 0xf0fe4786ul);
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x0fe1edc6ul);
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x240cf254ul);
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x4fe9346ful);
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x6cc984beul);
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x61b9411eul);
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x16f988faul);
        Round_noasm(a, b, c, &d, e, f, g, &h, 0xf2c65152ul);
        Round_noasm(h, a, b, &c, d, e, f, &g, 0xa88e5a6dul);
        Round_noasm(g, h, a, &b, c, d, e, &f, 0xb019fc65ul);
        Round_noasm(f, g, h, &a, b, c, d, &e, 0xb9d99ec7ul);
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x9a1231c3ul);
        Round_noasm(d, e, f, &g, h, a, b, &c, 0xe70eeaa0ul);
        Round_noasm(c, d, e, &f, g, h, a, &b, 0xfdb1232bul);
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xc7353eb0ul);
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x3069bad5ul);
        Round_noasm(h, a, b, &c, d, e, f, &g, 0xcb976d5ful);
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x5a0f118ful);
        Round_noasm(f, g, h, &a, b, c, d, &e, 0xdc1eeefdul);
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x0a35b689ul);
        Round_noasm(d, e, f, &g, h, a, b, &c, 0xde0b7a04ul);
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x58f4ca9dul);
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xe15d5b16ul);
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x007f3e86ul);
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x37088980ul);
        Round_noasm(g, h, a, &b, c, d, e, &f, 0xa507ea32ul);
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x6fab9537ul);
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x17406110ul);
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x0d8cd6f1ul);
        Round_noasm(c, d, e, &f, g, h, a, &b, 0xcdaa3b6dul);
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xc0bbbe37ul);
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x83613bdaul);
        Round_noasm(h, a, b, &c, d, e, f, &g, 0xdb48a363ul);
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x0b02e931ul);
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x6fd15ca7ul);
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x521afacaul);
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x31338431ul);
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x6ed41a95ul);
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x6d437890ul);
        Round_noasm(a, b, c, &d, e, f, g, &h, 0xc39c91f2ul);
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x9eccabbdul);
        Round_noasm(g, h, a, &b, c, d, e, &f, 0xb5c9a0e6ul);
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x532fb63cul);
        Round_noasm(e, f, g, &h, a, b, c, &d, 0xd2c741c6ul);
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x07237ea3ul);
        Round_noasm(c, d, e, &f, g, h, a, &b, 0xa4954b68ul);
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x4c191d76ul);

        w0 = t0 + a;
        w1 = t1 + b;
        w2 = t2 + c;
        w3 = t3 + d;
        w4 = t4 + e;
        w5 = t5 + f;
        w6 = t6 + g;
        w7 = t7 + h;

        /* Transform 3 */
        a = 0x6a09e667ul;
        b = 0xbb67ae85ul;
        c = 0x3c6ef372ul;
        d = 0xa54ff53aul;
        e = 0x510e527ful;
        f = 0x9b05688cul;
        g = 0x1f83d9abul;
        h = 0x5be0cd19ul;

        Round_noasm(a, b, c, &d, e, f, g, &h, 0x428a2f98ul + w0);
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x71374491ul + w1);
        Round_noasm(g, h, a, &b, c, d, e, &f, 0xb5c0fbcful + w2);
        Round_noasm(f, g, h, &a, b, c, d, &e, 0xe9b5dba5ul + w3);
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x3956c25bul + w4);
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x59f111f1ul + w5);
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x923f82a4ul + w6);
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xab1c5ed5ul + w7);
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x5807aa98ul);
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x12835b01ul);
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x243185beul);
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x550c7dc3ul);
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x72be5d74ul);
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x80deb1feul);
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x9bdc06a7ul);
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xc19bf274ul);
        Round_noasm(a, b, c, &d, e, f, g, &h, 0xe49b69c1ul + (w0 += sigma0_noasm(w1)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0xefbe4786ul + (w1 += 0xa00000ul + sigma0_noasm(w2)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x0fc19dc6ul + (w2 += sigma1_noasm(w0) + sigma0_noasm(w3)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x240ca1ccul + (w3 += sigma1_noasm(w1) + sigma0_noasm(w4)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x2de92c6ful + (w4 += sigma1_noasm(w2) + sigma0_noasm(w5)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x4a7484aaul + (w5 += sigma1_noasm(w3) + sigma0_noasm(w6)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x5cb0a9dcul + (w6 += sigma1_noasm(w4) + 0x100ul + sigma0_noasm(w7)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x76f988daul + (w7 += sigma1_noasm(w5) + w0 + 0x11002000ul));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x983e5152ul + (w8 = 0x80000000ul + sigma1_noasm(w6) + w1));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0xa831c66dul + (w9 = sigma1_noasm(w7) + w2));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0xb00327c8ul + (w10 = sigma1_noasm(w8) + w3));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0xbf597fc7ul + (w11 = sigma1_noasm(w9) + w4));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0xc6e00bf3ul + (w12 = sigma1_noasm(w10) + w5));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0xd5a79147ul + (w13 = sigma1_noasm(w11) + w6));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x06ca6351ul + (w14 = sigma1_noasm(w12) + w7 + 0x400022ul));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x14292967ul + (w15 = 0x100ul + sigma1_noasm(w13) + w8 + sigma0_noasm(w0)));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x27b70a85ul + (w0 += sigma1_noasm(w14) + w9 + sigma0_noasm(w1)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x2e1b2138ul + (w1 += sigma1_noasm(w15) + w10 + sigma0_noasm(w2)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x4d2c6dfcul + (w2 += sigma1_noasm(w0) + w11 + sigma0_noasm(w3)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x53380d13ul + (w3 += sigma1_noasm(w1) + w12 + sigma0_noasm(w4)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x650a7354ul + (w4 += sigma1_noasm(w2) + w13 + sigma0_noasm(w5)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x766a0abbul + (w5 += sigma1_noasm(w3) + w14 + sigma0_noasm(w6)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x81c2c92eul + (w6 += sigma1_noasm(w4) + w15 + sigma0_noasm(w7)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x92722c85ul + (w7 += sigma1_noasm(w5) + w0 + sigma0_noasm(w8)));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0xa2bfe8a1ul + (w8 += sigma1_noasm(w6) + w1 + sigma0_noasm(w9)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0xa81a664bul + (w9 += sigma1_noasm(w7) + w2 + sigma0_noasm(w10)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0xc24b8b70ul + (w10 += sigma1_noasm(w8) + w3 + sigma0_noasm(w11)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0xc76c51a3ul + (w11 += sigma1_noasm(w9) + w4 + sigma0_noasm(w12)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0xd192e819ul + (w12 += sigma1_noasm(w10) + w5 + sigma0_noasm(w13)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0xd6990624ul + (w13 += sigma1_noasm(w11) + w6 + sigma0_noasm(w14)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0xf40e3585ul + (w14 += sigma1_noasm(w12) + w7 + sigma0_noasm(w15)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x106aa070ul + (w15 += sigma1_noasm(w13) + w8 + sigma0_noasm(w0)));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x19a4c116ul + (w0 += sigma1_noasm(w14) + w9 + sigma0_noasm(w1)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x1e376c08ul + (w1 += sigma1_noasm(w15) + w10 + sigma0_noasm(w2)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x2748774cul + (w2 += sigma1_noasm(w0) + w11 + sigma0_noasm(w3)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x34b0bcb5ul + (w3 += sigma1_noasm(w1) + w12 + sigma0_noasm(w4)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x391c0cb3ul + (w4 += sigma1_noasm(w2) + w13 + sigma0_noasm(w5)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0x4ed8aa4aul + (w5 += sigma1_noasm(w3) + w14 + sigma0_noasm(w6)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0x5b9cca4ful + (w6 += sigma1_noasm(w4) + w15 + sigma0_noasm(w7)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0x682e6ff3ul + (w7 += sigma1_noasm(w5) + w0 + sigma0_noasm(w8)));
        Round_noasm(a, b, c, &d, e, f, g, &h, 0x748f82eeul + (w8 += sigma1_noasm(w6) + w1 + sigma0_noasm(w9)));
        Round_noasm(h, a, b, &c, d, e, f, &g, 0x78a5636ful + (w9 += sigma1_noasm(w7) + w2 + sigma0_noasm(w10)));
        Round_noasm(g, h, a, &b, c, d, e, &f, 0x84c87814ul + (w10 += sigma1_noasm(w8) + w3 + sigma0_noasm(w11)));
        Round_noasm(f, g, h, &a, b, c, d, &e, 0x8cc70208ul + (w11 += sigma1_noasm(w9) + w4 + sigma0_noasm(w12)));
        Round_noasm(e, f, g, &h, a, b, c, &d, 0x90befffaul + (w12 += sigma1_noasm(w10) + w5 + sigma0_noasm(w13)));
        Round_noasm(d, e, f, &g, h, a, b, &c, 0xa4506cebul + (w13 += sigma1_noasm(w11) + w6 + sigma0_noasm(w14)));
        Round_noasm(c, d, e, &f, g, h, a, &b, 0xbef9a3f7ul + (w14 + sigma1_noasm(w12) + w7 + sigma0_noasm(w15)));
        Round_noasm(b, c, d, &e, f, g, h, &a, 0xc67178f2ul + (w15 + sigma1_noasm(w13) + w8 + sigma0_noasm(w0)));

        /* Output */
        WriteBE32(&out->u8[0], a + 0x6a09e667ul);
        WriteBE32(&out->u8[4], b + 0xbb67ae85ul);
        WriteBE32(&out->u8[8], c + 0x3c6ef372ul);
        WriteBE32(&out->u8[12], d + 0xa54ff53aul);
        WriteBE32(&out->u8[16], e + 0x510e527ful);
        WriteBE32(&out->u8[20], f + 0x9b05688cul);
        WriteBE32(&out->u8[24], g + 0x1f83d9abul);
        WriteBE32(&out->u8[28], h + 0x5be0cd19ul);
}

typedef void (*transform_t)(uint32_t*, const unsigned char*, size_t);
typedef void (*transform_multi_t)(struct sha256*, const uint32_t*, const unsigned char*);
typedef void (*transform_d64_t)(struct sha256[], const struct sha256[]);

void transform_d64_wrapper(struct sha256 out[1], const struct sha256 in[2], transform_t tr)
{
        uint32_t s[8];
        static const unsigned char padding1[64] = {
                0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0
        };
        unsigned char buffer2[64] = {
                0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0
        };
        Initialize(s);
        tr(s, in->u8, 1);
        tr(s, padding1, 1);
        WriteBE32(buffer2 + 0, s[0]);
        WriteBE32(buffer2 + 4, s[1]);
        WriteBE32(buffer2 + 8, s[2]);
        WriteBE32(buffer2 + 12, s[3]);
        WriteBE32(buffer2 + 16, s[4]);
        WriteBE32(buffer2 + 20, s[5]);
        WriteBE32(buffer2 + 24, s[6]);
        WriteBE32(buffer2 + 28, s[7]);
        Initialize(s);
        tr(s, buffer2, 1);
        WriteBE32(&out->u8[0], s[0]);
        WriteBE32(&out->u8[4], s[1]);
        WriteBE32(&out->u8[8], s[2]);
        WriteBE32(&out->u8[12], s[3]);
        WriteBE32(&out->u8[16], s[4]);
        WriteBE32(&out->u8[20], s[5]);
        WriteBE32(&out->u8[24], s[6]);
        WriteBE32(&out->u8[28], s[7]);
}
#if defined(__x86_64__) || defined(__amd64__)
void transform_sha256d64_shani(struct sha256 out[1], const struct sha256 in[2])
{
        transform_d64_wrapper(out, in, transform_sha256_shani);
}
void transform_sha256d64_sse4(struct sha256 out[1], const struct sha256 in[2])
{
        transform_d64_wrapper(out, in, transform_sha256_sse4);
}
#endif /* defined(__x86_64__) || defined(__amd64__) || defined(__i386__) */
#if defined(__arm__) || defined(__aarch32__) || defined(__arm64__) || defined(__aarch64__) || defined(_M_ARM)
void transform_sha256d64_armv8(struct sha256 out[1], const struct sha256 in[2])
{
        transform_d64_wrapper(out, in, transform_sha256_armv8);
}
#endif /* defined(__arm__) || defined(__aarch32__) || defined(__arm64__) || defined(__aarch64__) || defined(_M_ARM) */

transform_t transform = transform_noasm;
transform_multi_t transform_2way = NULL;
transform_multi_t transform_4way = NULL;
transform_multi_t transform_8way = NULL;
transform_d64_t transform_d64 = transform_d64_noasm;
transform_d64_t transform_d64_2way = NULL;
transform_d64_t transform_d64_4way = NULL;
transform_d64_t transform_d64_8way = NULL;

#ifndef NDEBUG
static int self_test() {
        /* Input state (equal to the initial SHA256 state) */
        static const uint32_t init[8] = {
                0x6a09e667ul, 0xbb67ae85ul, 0x3c6ef372ul, 0xa54ff53aul, 0x510e527ful, 0x9b05688cul, 0x1f83d9abul, 0x5be0cd19ul
        };
        /* Some random input data to test with */
        static const unsigned char data[641] = "-" /* Intentionally not aligned */
                "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do "
                "eiusmod tempor incididunt ut labore et dolore magna aliqua. Et m"
                "olestie ac feugiat sed lectus vestibulum mattis ullamcorper. Mor"
                "bi blandit cursus risus at ultrices mi tempus imperdiet nulla. N"
                "unc congue nisi vita suscipit tellus mauris. Imperdiet proin fer"
                "mentum leo vel orci. Massa tempor nec feugiat nisl pretium fusce"
                " id velit. Telus in metus vulputate eu scelerisque felis. Mi tem"
                "pus imperdiet nulla malesuada pellentesque. Tristique magna sit.";
        static const struct sha256* data_d64 = (const struct sha256*)&data[1];
        /* Expected output state for hashing the i*64 first input bytes above (excluding SHA256 padding). */
        static const uint32_t result[9][8] = {
                {0x6a09e667ul, 0xbb67ae85ul, 0x3c6ef372ul, 0xa54ff53aul, 0x510e527ful, 0x9b05688cul, 0x1f83d9abul, 0x5be0cd19ul},
                {0x91f8ec6bul, 0x4da10fe3ul, 0x1c9c292cul, 0x45e18185ul, 0x435cc111ul, 0x3ca26f09ul, 0xeb954caeul, 0x402a7069ul},
                {0xcabea5acul, 0x374fb97cul, 0x182ad996ul, 0x7bd69cbful, 0x450ff900ul, 0xc1d2be8aul, 0x6a41d505ul, 0xe6212dc3ul},
                {0xbcff09d6ul, 0x3e76f36eul, 0x3ecb2501ul, 0x78866e97ul, 0xe1c1e2fdul, 0x32f4eafful, 0x8aa6c4e5ul, 0xdfc024bcul},
                {0xa08c5d94ul, 0x0a862f93ul, 0x6b7f2f40ul, 0x8f9fae76ul, 0x6d40439ful, 0x79dcee0cul, 0x3e39ff3aul, 0xdc3bdbb1ul},
                {0x216a0895ul, 0x9f1a3662ul, 0xe99946f9ul, 0x87ba4364ul, 0x0fb5db2cul, 0x12bed3d3ul, 0x6689c0c7ul, 0x292f1b04ul},
                {0xca3067f8ul, 0xbc8c2656ul, 0x37cb7e0dul, 0x9b6b8b0ful, 0x46dc380bul, 0xf1287f57ul, 0xc42e4b23ul, 0x3fefe94dul},
                {0x3e4c4039ul, 0xbb6fca8cul, 0x6f27d2f7ul, 0x301e44a4ul, 0x8352ba14ul, 0x5769ce37ul, 0x48a1155ful, 0xc0e1c4c6ul},
                {0xfe2fa9ddul, 0x69d0862bul, 0x1ae0db23ul, 0x471f9244ul, 0xf55c0145ul, 0xc30f9c3bul, 0x40a84ea0ul, 0x5b8a266cul},
        };
        /* Expected output for each of the individual 8 64-byte messages under full double SHA256 (including padding). */
        static const unsigned char result_d64[256] = {
                0x09, 0x3a, 0xc4, 0xd0, 0x0f, 0xf7, 0x57, 0xe1, 0x72, 0x85, 0x79, 0x42, 0xfe, 0xe7, 0xe0, 0xa0,
                0xfc, 0x52, 0xd7, 0xdb, 0x07, 0x63, 0x45, 0xfb, 0x53, 0x14, 0x7d, 0x17, 0x22, 0x86, 0xf0, 0x52,
                0x48, 0xb6, 0x11, 0x9e, 0x6e, 0x48, 0x81, 0x6d, 0xcc, 0x57, 0x1f, 0xb2, 0x97, 0xa8, 0xd5, 0x25,
                0x9b, 0x82, 0xaa, 0x89, 0xe2, 0xfd, 0x2d, 0x56, 0xe8, 0x28, 0x83, 0x0b, 0xe2, 0xfa, 0x53, 0xb7,
                0xd6, 0x6b, 0x07, 0x85, 0x83, 0xb0, 0x10, 0xa2, 0xf5, 0x51, 0x3c, 0xf9, 0x60, 0x03, 0xab, 0x45,
                0x6c, 0x15, 0x6e, 0xef, 0xb5, 0xac, 0x3e, 0x6c, 0xdf, 0xb4, 0x92, 0x22, 0x2d, 0xce, 0xbf, 0x3e,
                0xe9, 0xe5, 0xf6, 0x29, 0x0e, 0x01, 0x4f, 0xd2, 0xd4, 0x45, 0x65, 0xb3, 0xbb, 0xf2, 0x4c, 0x16,
                0x37, 0x50, 0x3c, 0x6e, 0x49, 0x8c, 0x5a, 0x89, 0x2b, 0x1b, 0xab, 0xc4, 0x37, 0xd1, 0x46, 0xe9,
                0x3d, 0x0e, 0x85, 0xa2, 0x50, 0x73, 0xa1, 0x5e, 0x54, 0x37, 0xd7, 0x94, 0x17, 0x56, 0xc2, 0xd8,
                0xe5, 0x9f, 0xed, 0x4e, 0xae, 0x15, 0x42, 0x06, 0x0d, 0x74, 0x74, 0x5e, 0x24, 0x30, 0xce, 0xd1,
                0x9e, 0x50, 0xa3, 0x9a, 0xb8, 0xf0, 0x4a, 0x57, 0x69, 0x78, 0x67, 0x12, 0x84, 0x58, 0xbe, 0xc7,
                0x36, 0xaa, 0xee, 0x7c, 0x64, 0xa3, 0x76, 0xec, 0xff, 0x55, 0x41, 0x00, 0x2a, 0x44, 0x68, 0x4d,
                0xb6, 0x53, 0x9e, 0x1c, 0x95, 0xb7, 0xca, 0xdc, 0x7f, 0x7d, 0x74, 0x27, 0x5c, 0x8e, 0xa6, 0x84,
                0xb5, 0xac, 0x87, 0xa9, 0xf3, 0xff, 0x75, 0xf2, 0x34, 0xcd, 0x1a, 0x3b, 0x82, 0x2c, 0x2b, 0x4e,
                0x6a, 0x46, 0x30, 0xa6, 0x89, 0x86, 0x23, 0xac, 0xf8, 0xa5, 0x15, 0xe9, 0x0a, 0xaa, 0x1e, 0x9a,
                0xd7, 0x93, 0x6b, 0x28, 0xe4, 0x3b, 0xfd, 0x59, 0xc6, 0xed, 0x7c, 0x5f, 0xa5, 0x41, 0xcb, 0x51
        };

        int i;

        /* Test transform() for 0 through 8 transformations. */
        for (i = 0; i <= 8; ++i) {
                uint32_t state[8];
                memcpy(state, init, 8 * sizeof(uint32_t));
                transform(state, data + 1, i);
                if (memcmp(state, result[i], 8 * sizeof(uint32_t))) return 0;
        }

        /* Test transform_d64 */
        {
            struct sha256 out[1];
            transform_d64(out, data_d64);
            if (memcmp(out, result_d64, 32)) return 0;
        }

        /* Test transform_d64_2way, if available. */
        if (transform_d64_2way) {
                struct sha256 out[2];
                transform_d64_2way(out, data_d64);
                if (memcmp(out, result_d64, 64)) return 0;
        }

        /* Test transform_d64_4way, if available. */
        if (transform_d64_4way) {
                struct sha256 out[4];
                transform_d64_4way(out, data_d64);
                if (memcmp(out, result_d64, 128)) return 0;
        }

        /* Test transform_d64_8way, if available. */
        if (transform_d64_8way) {
                struct sha256 out[8];
                transform_d64_8way(out, data_d64);
                if (memcmp(out, result_d64, 256)) return 0;
        }

        return !0;
}
#endif /* NDEBUG */

#if (defined(__x86_64__) || defined(__amd64__) || defined(__i386__))
/** Check whether the OS has enabled AVX registers. */
static int AVXEnabled()
{
        uint32_t a, d;
        __asm__("xgetbv" : "=a"(a), "=d"(d) : "c"(0));
        return (a & 6) == 6;
}
#endif

const char* sha256_auto_detect()
{
        static char ret[255] = "standard";
#if defined(HAVE_GETCPUID)
        int have_sse4 = 0;
        int have_xsave = 0;
        int have_avx = 0;
        int have_avx2 = 0;
        int have_shani = 0;
        int enabled_avx = 0;

        uint32_t eax=0, ebx=0, ecx=0, edx=0;

        (void)AVXEnabled;
        (void)have_sse4;
        (void)have_avx;
        (void)have_xsave;
        (void)have_avx2;
        (void)have_shani;
        (void)enabled_avx;

        GetCPUID(1, 0, &eax, &ebx, &ecx, &edx);
        have_sse4 = (ecx >> 19) & 1;
        have_xsave = (ecx >> 27) & 1;
        have_avx = (ecx >> 28) & 1;
        if (have_xsave && have_avx) {
                enabled_avx = AVXEnabled();
        }
        if (have_sse4) {
                GetCPUID(7, 0, &eax, &ebx, &ecx, &edx);
                have_avx2 = (ebx >> 5) & 1;
                have_shani = (ebx >> 29) & 1;
        }

#if !defined(BUILD_BITCOIN_INTERNAL)
        if (have_shani) {
                transform = transform_sha256_shani;
                transform_d64 = transform_sha256d64_shani;
                transform_d64_2way = transform_sha256d64_shani_2way;
                strcpy(ret, "shani(1way,2way)");
                have_sse4 = 0; /* Disable SSE4/AVX2; */
                have_avx2 = 0;
        }
#endif

        if (have_sse4) {
#if defined(__x86_64__) || defined(__amd64__)
                transform = transform_sha256_sse4;
                transform_d64 = transform_sha256d64_sse4;
                strcpy(ret, "sse4(1way)");
#endif
#if !defined(BUILD_BITCOIN_INTERNAL)
                transform_4way = transform_sha256multi_sse41_4way;
                transform_d64_4way = transform_sha256d64_sse41_4way;
                strcat(ret, ",sse41(4way)");
#endif
    }

#if !defined(BUILD_BITCOIN_INTERNAL)
        if (have_avx2 && have_avx && enabled_avx) {
                transform_8way = transform_sha256multi_avx2_8way;
                transform_d64_8way = transform_sha256d64_avx2_8way;
                strcat(ret, ",avx2(8way)");
        }
#endif

#elif defined(__aarch64__)
        int have_arm_shani = 0;

#if defined(__linux__)
#if defined(__arm__) /* 32-bit */
        if (getauxval(AT_HWCAP2) & HWCAP2_SHA2) {
                have_arm_shani = !0;
        }
#endif
#if defined(__aarch64__) /* 64-bit */
        if (getauxval(AT_HWCAP) & HWCAP_SHA2) {
                have_arm_shani = !0;
        }
#endif
#endif

#if defined(__APPLE__)
        int val = 0;
        size_t len = sizeof(val);
        if (sysctlbyname("hw.optional.arm.FEAT_SHA256", &val, &len, NULL, 0) == 0) {
                have_arm_shani = !0;
        }
#endif

        if (have_arm_shani) {
                transform = transform_sha256_armv8;
                transform_d64 = transform_sha256d64_armv8;
                transform_d64_2way = transform_sha256d64_armv8_2way;
                strcpy(ret, "armv8(1way,2way)");
        }
#endif

        assert(self_test());
        return ret;
}

/* SHA-256 */

void sha256_init(struct sha256_ctx* ctx)
{
        assert(ctx);
        ctx->bytes = 0;
        Initialize(ctx->s);
}

void sha256_update(struct sha256_ctx* ctx, const void *_data, size_t len)
{
        const unsigned char* data = (const unsigned char*)_data;
        const unsigned char* end = data + len;
        size_t bufsize = ctx->bytes % 64;
        if (bufsize && bufsize + len >= 64) {
                /* Fill the buffer, and process it. */
                memcpy(ctx->buf.u8 + bufsize, data, 64 - bufsize);
                ctx->bytes += 64 - bufsize;
                data += 64 - bufsize;
                transform(ctx->s, ctx->buf.u8, 1);
                bufsize = 0;
        }
        if (end - data >= 64) {
                size_t blocks = (end - data) / 64;
                transform(ctx->s, data, blocks);
                data += 64 * blocks;
                ctx->bytes += 64 * blocks;
        }
        if (end > data) {
                /* Fill the buffer with what remains. */
                memcpy(ctx->buf.u8 + bufsize, data, end - data);
                ctx->bytes += end - data;
        }
}

void sha256_done(struct sha256* hash, struct sha256_ctx* ctx)
{
        static const unsigned char pad[64] = {0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        unsigned char sizedesc[8];
        WriteBE64(sizedesc, ctx->bytes << 3);
        sha256_update(ctx, pad, 1 + ((119 - (ctx->bytes % 64)) % 64));
        sha256_update(ctx, sizedesc, 8);
        WriteBE32(&hash->u8[0], ctx->s[0]);
        WriteBE32(&hash->u8[4], ctx->s[1]);
        WriteBE32(&hash->u8[8], ctx->s[2]);
        WriteBE32(&hash->u8[12], ctx->s[3]);
        WriteBE32(&hash->u8[16], ctx->s[4]);
        WriteBE32(&hash->u8[20], ctx->s[5]);
        WriteBE32(&hash->u8[24], ctx->s[6]);
        WriteBE32(&hash->u8[28], ctx->s[7]);
}

void sha256_double64(struct sha256 out[], const struct sha256 in[], size_t blocks)
{
        if (transform_d64_8way) {
                while (blocks >= 8) {
                        transform_d64_8way(out, in);
                        out += 8;
                        in += 16;
                        blocks -= 8;
                }
        }
        if (transform_d64_4way) {
                while (blocks >= 4) {
                        transform_d64_4way(out, in);
                        out += 4;
                        in += 8;
                        blocks -= 4;
                }
        }
        if (transform_d64_2way) {
                while (blocks >= 2) {
                        transform_d64_2way(out, in);
                        out += 2;
                        in += 4;
                        blocks -= 2;
                }
        }
        while (blocks) {
                transform_d64(out, in);
                ++out;
                in += 2;
                --blocks;
        }
}

void sha256_midstate(struct sha256 out[], const uint32_t midstate[8], const unsigned char in[], size_t blocks)
{
        if (transform_8way) {
                while (blocks >= 8) {
                        transform_8way(out, midstate, in);
                        out += 8;
                        in += 512;
                        blocks -= 8;
                }
        }
        if (transform_4way) {
                while (blocks >= 4) {
                        transform_4way(out, midstate, in);
                        out += 4;
                        in += 256;
                        blocks -= 4;
                }
        }
        if (transform_2way) {
                while (blocks >= 2) {
                        transform_2way(out, midstate, in);
                        out += 2;
                        in += 128;
                        blocks -= 2;
                }
        }
        while (blocks) {
                uint32_t s[8];
                unsigned char* _out = out->u8;
                int i;
                memcpy(s, midstate, 8 * sizeof(uint32_t));
                transform(s, in, 1);
                for (i = 0; i < 8; ++i) {
                        WriteBE32(_out, s[i]);
                        _out += 4;
                }
                out += 1;
                in += 64;
                --blocks;
        }
}

/* End of File
 */
