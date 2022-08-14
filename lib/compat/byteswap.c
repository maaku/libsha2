/* Copyright (c) 2014-2019 The Bitcoin Core developers
 * Copyright (c) 2022 Mark Friedenbach
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "byteswap.h"

#if !defined(__APPLE__)
/* Non-MacOS / non-Darwin */

#include <stdint.h> /* for uint16_t, uint32_t, uint64_t */

#if HAVE_DECL_BSWAP_16 == 0
uint16_t bswap_16(uint16_t x)
{
        return (x >> 8) | (x << 8);
}
#endif /* HAVE_DECL_BSWAP16 == 0 */

#if HAVE_DECL_BSWAP_32 == 0
uint32_t bswap_32(uint32_t x)
{
        return (((x & 0xff000000U) >> 24) | ((x & 0x00ff0000U) >>  8) |
                ((x & 0x0000ff00U) <<  8) | ((x & 0x000000ffU) << 24));
}
#endif /* HAVE_DECL_BSWAP32 == 0 */

#if HAVE_DECL_BSWAP_64 == 0
uint64_t bswap_64(uint64_t x)
{
         return (((x & 0xff00000000000000ull) >> 56)
               | ((x & 0x00ff000000000000ull) >> 40)
               | ((x & 0x0000ff0000000000ull) >> 24)
               | ((x & 0x000000ff00000000ull) >>  8)
               | ((x & 0x00000000ff000000ull) <<  8)
               | ((x & 0x0000000000ff0000ull) << 24)
               | ((x & 0x000000000000ff00ull) << 40)
               | ((x & 0x00000000000000ffull) << 56));
}
#endif /* HAVE_DECL_BSWAP64 == 0 */

#endif /* !defined(__APPLE__) */

/* End of File
 */
