/* Copyright (c) 2014-2019 The Bitcoin Core developers
 * Copyright (c) 2022 Mark Friedenbach
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SHA2__COMPAT__BYTESWAP_H
#define SHA2__COMPAT__BYTESWAP_H

#if defined(HAVE_CONFIG_H)
#include <libsha2-config.h>
#endif

#include <stdint.h> /* for uint16_t, uint32_t, uint64_t */

#if defined(HAVE_BYTESWAP_H)
#include <byteswap.h>
#endif

#if defined(__APPLE__)

#include <libkern/OSByteOrder.h>
#define bswap_16(x) OSSwapInt16(x)
#define bswap_32(x) OSSwapInt32(x)
#define bswap_64(x) OSSwapInt64(x)

#else
/* Non-MacOS / non-Darwin */

#if HAVE_DECL_BSWAP_16 == 0
uint16_t bswap_16(uint16_t x);
#endif /* HAVE_DECL_BSWAP16 == 0 */

#if HAVE_DECL_BSWAP_32 == 0
uint32_t bswap_32(uint32_t x);
#endif /* HAVE_DECL_BSWAP32 == 0 */

#if HAVE_DECL_BSWAP_64 == 0
uint64_t bswap_64(uint64_t x);
#endif /* HAVE_DECL_BSWAP64 == 0 */

#endif /* defined(__APPLE__) */

#endif /* SHA2__COMPAT__BYTESWAP_H */

/* End of File
 */
