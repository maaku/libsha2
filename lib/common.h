/* Copyright (c) 2014-2020 The Bitcoin Core developers
 * Copyright (c) 2022 Mark Friedenbach
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SHA2__COMMON_H
#define SHA2__COMMON_H

#if defined(HAVE_CONFIG_H)
#include <config/bitcoin-config.h>
#endif

#include <stdint.h> /* for uint16_t, uint32_t, uint64_t */

uint16_t ReadLE16(const unsigned char* ptr);
uint32_t ReadLE32(const unsigned char* ptr);
uint64_t ReadLE64(const unsigned char* ptr);

void WriteLE16(unsigned char* ptr, uint16_t x);
void WriteLE32(unsigned char* ptr, uint32_t x);
void WriteLE64(unsigned char* ptr, uint64_t x);

uint16_t ReadBE16(const unsigned char* ptr);
uint32_t ReadBE32(const unsigned char* ptr);
uint64_t ReadBE64(const unsigned char* ptr);

void WriteBE16(unsigned char* ptr, uint16_t x);
void WriteBE32(unsigned char* ptr, uint32_t x);
void WriteBE64(unsigned char* ptr, uint64_t x);

/**
 * @brief Retuen the 1-based index of the highest set bit.
 * 
 * @param x any 64-bit number
 * @return uint64_t the highest set bit, or 0 if no bits are set
 * 
 * Return the smallest number n such that (x >> n) == 0 (or 64 if the highest bit in x is set. 
 */
uint64_t CountBits(uint64_t x);

#endif /* SHA2__COMMON_H */

/* End of File
 */
