/* Copyright (c) 2014-2018 The Bitcoin Core developers
 * Copyright (c) 2022 Mark Friedenbach
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SHA2__COMPAT__ENDIAN_H
#define SHA2__COMPAT__ENDIAN_H

#if defined(HAVE_CONFIG_H)
#include <config/bitcoin-config.h>
#endif

#include "byteswap.h" /* for bswap_16, bswap_32, bswap_64 */

#include <stdint.h> /* for uint16_t, uint32_t, uint64_t */

#if defined(HAVE_ENDIAN_H)
#include <endian.h>
#elif defined(HAVE_SYS_ENDIAN_H)
#include <sys/endian.h>
#endif

#ifndef HAVE_CONFIG_H
/* While not technically a supported configuration, defaulting to defining these
 * DECLs when we were compiled without autotools makes it easier for other build
 * systems to build things like libbitcoinconsensus for strange targets. */
#ifdef htobe16
#define HAVE_DECL_HTOBE16 1
#endif
#ifdef htole16
#define HAVE_DECL_HTOLE16 1
#endif
#ifdef be16toh
#define HAVE_DECL_BE16TOH 1
#endif
#ifdef le16toh
#define HAVE_DECL_LE16TOH 1
#endif

#ifdef htobe32
#define HAVE_DECL_HTOBE32 1
#endif
#ifdef htole32
#define HAVE_DECL_HTOLE32 1
#endif
#ifdef be32toh
#define HAVE_DECL_BE32TOH 1
#endif
#ifdef le32toh
#define HAVE_DECL_LE32TOH 1
#endif

#ifdef htobe64
#define HAVE_DECL_HTOBE64 1
#endif
#ifdef htole64
#define HAVE_DECL_HTOLE64 1
#endif
#ifdef be64toh
#define HAVE_DECL_BE64TOH 1
#endif
#ifdef le64toh
#define HAVE_DECL_LE64TOH 1
#endif

#endif /* HAVE_CONFIG_H */

#if defined(WORDS_BIGENDIAN)

#if HAVE_DECL_HTOBE16 == 0
#define htobe16(host_16bits) host_16bits
#endif /* HAVE_DECL_HTOBE16 */

#if HAVE_DECL_HTOLE16 == 0
#define htole16(host_16bits) bswap_16(host_16bits)
#endif /* HAVE_DECL_HTOLE16 */

#if HAVE_DECL_BE16TOH == 0
#define be16toh(big_endian_16bits) big_endian_16bits
#endif /* HAVE_DECL_BE16TOH */

#if HAVE_DECL_LE16TOH == 0
#define le16toh(little_endian_16bits) bswap_16(little_endian_16bits)
#endif /* HAVE_DECL_LE16TOH */

#if HAVE_DECL_HTOBE32 == 0
#define htobe32(host_32bits) host_32bits
#endif /* HAVE_DECL_HTOBE32 */

#if HAVE_DECL_HTOLE32 == 0
#define htole32(host_32bits) bswap_32(host_32bits)
#endif /* HAVE_DECL_HTOLE32 */

#if HAVE_DECL_BE32TOH == 0
#define be32toh(big_endian_32bits) big_endian_32bits
#endif /* HAVE_DECL_BE32TOH */

#if HAVE_DECL_LE32TOH == 0
#define le32toh(little_endian_32bits) bswap_32(little_endian_32bits)
#endif /* HAVE_DECL_LE32TOH */

#if HAVE_DECL_HTOBE64 == 0
#define htobe64(uint64_t host_64bits) host_64bits
#endif /* HAVE_DECL_HTOBE64 */

#if HAVE_DECL_HTOLE64 == 0
#define htole64(uint64_t host_64bits) bswap_64(host_64bits)
#endif /* HAVE_DECL_HTOLE64 */

#if HAVE_DECL_BE64TOH == 0
#define be64toh(uint64_t big_endian_64bits) big_endian_64bits
#endif /* HAVE_DECL_BE64TOH */

#if HAVE_DECL_LE64TOH == 0
#define le64toh(uint64_t little_endian_64bits) bswap_64(little_endian_64bits)
#endif /* HAVE_DECL_LE64TOH */

#else /* WORDS_BIGENDIAN */

#if HAVE_DECL_HTOBE16 == 0
#define htobe16(host_16bits) bswap_16(host_16bits)
#endif /* HAVE_DECL_HTOBE16 */

#if HAVE_DECL_HTOLE16 == 0
#define htole16(host_16bits) host_16bits
#endif /* HAVE_DECL_HTOLE16 */

#if HAVE_DECL_BE16TOH == 0
#define be16toh(big_endian_16bits) bswap_16(big_endian_16bits)
#endif /* HAVE_DECL_BE16TOH */

#if HAVE_DECL_LE16TOH == 0
#define le16toh(little_endian_16bits) little_endian_16bits
#endif /* HAVE_DECL_LE16TOH */

#if HAVE_DECL_HTOBE32 == 0
#define htobe32(host_32bits) bswap_32(host_32bits)
#endif /* HAVE_DECL_HTOBE32 */

#if HAVE_DECL_HTOLE32 == 0
#define htole32(host_32bits) host_32bits
#endif /* HAVE_DECL_HTOLE32 */

#if HAVE_DECL_BE32TOH == 0
#define be32toh(big_endian_32bits) bswap_32(big_endian_32bits)
#endif /* HAVE_DECL_BE32TOH */

#if HAVE_DECL_LE32TOH == 0
#define le32toh(little_endian_32bits) little_endian_32bits
#endif /* HAVE_DECL_LE32TOH */

#if HAVE_DECL_HTOBE64 == 0
#define htobe64(host_64bits) bswap_64(host_64bits)
#endif /* HAVE_DECL_HTOBE64 */

#if HAVE_DECL_HTOLE64 == 0
#define htole64(host_64bits) host_64bits
#endif /* HAVE_DECL_HTOLE64 */

#if HAVE_DECL_BE64TOH == 0
#define be64toh(big_endian_64bits) bswap_64(big_endian_64bits)
#endif /* HAVE_DECL_BE64TOH */

#if HAVE_DECL_LE64TOH == 0
#define le64toh(little_endian_64bits) little_endian_64bits
#endif /* HAVE_DECL_LE64TOH */

#endif /* WORDS_BIGENDIAN */

#endif /* SHA2__COMPAT__ENDIAN_H */

/* End of File
 */
