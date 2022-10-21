/* Copyright (c) 2014-2018 The Bitcoin Core developers
 * Copyright (c) 2022 Mark Friedenbach
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SHA2__SHA256_H
#define SHA2__SHA256_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h> /* for uint32_t */
#include <stdlib.h> /* for size_t */

/**
 * @brief Autodetect the best available SHA256 implementation.
 *
 * @return const char* an ASCII-encoded string describing the selected
 * algorithm(s)
 *
 * This library supports various vector compute and special-purpose
 * cryptographic accelerators for efficiently calculating SHA256 hashes.  This
 * function queries the host machine capabilities and selects which backend
 * implementation to use.  A pointer to an internal buffer is returned to the
 * caller, containing a string naming the algorithms selected.
 *
 * This API is automatically called the first time a hash update is performed,
 * however it is not thread safe.  To query capabilities in a thread-safe way,
 * and to report the selected algorithm to the user, please call this function
 * during program initialization.
 */
const char* sha256_auto_detect(void);

/**
 * @brief A structure representing a completed SHA256 hash digest value.
 *
 * @u.u8: an unsigned char array
 * @u.u32: a 32-bit integer array
 *
 * A completed SHA256 hash digest is stored in big-endian byte order, regardless
 * of the endianness of the host machine.  Therefore accessing the union via
 * explicitly deserializing u8 is required for cross-platform compatibility.
 */
struct sha256 {
        unsigned char u8[32];
};

/**
 * @brief A structure for storing the running context of a sha256 hash.
 *
 * @s: the intermediate state in host-native byte order
 * @buf: a buffer of up to 63 bytes of unhashed data
 * @len: the total number of bytes hashed, including any buffered data
 *
 * The SHA256 state update function operates in blocks of 64 bytes at a time,
 * but for convenience data of any size can be delivered to the hasher.  The
 * context object contains the midstate data representing the current hash
 * state, and a buffer of length 64 (of which up to 63 bytes might be used), and
 * a count of the number of bytes hashed so far (including any buffered bytes).
 * The count of the total number of bytes hashes is used both to record the
 * number of buffered bytes and to calculate the final padding block.
 */
struct sha256_ctx {
        uint32_t s[8];
        union {
                uint32_t u32[16];
                unsigned char u8[64];
        } buf;
        size_t bytes;
};

/**
 * @brief Initializes a SHA256 context.
 *
 * @param ctx the context to initialize
 *
 * This must be called before any of the other sha256 functions which take a
 * context parameter.  Alternatively you may use the SHA256_INIT initialization
 * constant instead.
 *
 * If the context was already initialized, this forgets any data that was hashed
 * before.
 *
 * SHA256 context state does not include any dynamically allocated data, so
 * destruction of the context once finished is not required.  Contexts may also
 * be safely reset through re-initialization with sha256_init().
 *
 * Example:
 * static void hash_data(const char* data, size_t len, struct sha256* hash)
 * {
 *         struct sha256_ctx ctx;
 *         sha256_init(&ctx);
 *         sha256_update(&ctx, data, len);
 *         sha256_done(&ctx, hash);
 * }
 */
void sha256_init(struct sha256_ctx* ctx);

/**
 * @brief Initialization constant for a SHA256 context.
 *
 * This can be used to statically initialize a SHA256 context, equivalent to
 * calling sha256_init().
 *
 * Example:
 * static void hash_data(const char* data, size_t len, struct sha256* hash)
 * {
 *         struct sha256_ctx ctx = SHA256_INIT;
 *         sha256_update(&ctx, data, len);
 *         sha256_done(&ctx, hash);
 * }
 */
#define SHA256_INIT                                                   \
        { { 0x6a09e667ul, 0xbb67ae85ul, 0x3c6ef372ul, 0xa54ff53aul,   \
            0x510e527ful, 0x9b05688cul, 0x1f83d9abul, 0x5be0cd19ul }, \
          { { 0 } }, 0 }

/**
 * @brief Add some data from memory to the hash.
 *
 * @param ctx the sha256_ctx to use
 * @param data a pointer to data in memory
 * @param len the number of bytes pointed to by \p data
 *
 * Adds data to the hash context, performing hash compressions if a full block
 * of 64 bytes is formed, or storing the bytes in the buffer otherwise.
 *
 * You can call this function multiple times with the same context to hash more
 * data, before calling sha256_done().
 */
void sha256_update(struct sha256_ctx* ctx, const void *data, size_t len);

/**
 * @brief Finalize a SHA256 and return the resulting hash.
 *
 * @param hash the hash to return
 * @param ctx the sha256_ctx to finalize
 *
 * The finalization step involves hashing the remaining bytes in the buffer,
 * plus some padding bytes to finish out the last block and include the number
 * of bytes hashed.
 *
 * Note that the context is used up by this call and must be re-initialized
 * before being used again.  Calling sha256_update() or sha256_done() on an
 * already finalized context results in garbage data.
 */
void sha256_done(struct sha256* hash, struct sha256_ctx* ctx);

/**
 * @brief Perform a Merkle-tree compression step using double-SHA256
 *
 * @param out an array of 2*blocks sha256 hash values
 * @param in an array of 1*blocks sha256 hash values
 * @param blocks the number of double-SHA256 hash operations to perform
 *
 * Each pair of hashes from the input are hashed together, and then the
 * resulting value is hashed one more time before being written to the output.
 * Since two 32-byte hash values take two SHA256 blocks to compress, with
 * padding, and the intermediate value with padding can fit in a single block,
 * this performs 3 compression rounds per inner-node in the tree.  However
 * because of various optimizations that are available due to the fixed format
 * of the input, this only requires about 2.32x more computation compared with a
 * single round of SHA256. In addition, various multi-lane optimizations are
 * available allowing for up to 8 hashes to be performed simultaneously on some
 * architectures.
 *
 * The origin of this primitive is in how Bitcoin and related projects construct
 * Merkle trees by performing double-SHA256 hashes of child leaf values to
 * produce the parent or inner-node value.  Merkle hashing can take a
 * significant amount of time for various cryptocurrency applications, so this
 * library contains vector-optimized implementations for many architectures.
 */
void sha256_double64(struct sha256 out[], const struct sha256 in[], size_t blocks);

/**
 * @brief Performs multiple SHA256 compression rounds in parallel using the same
 * initial state vector but differing data blocks
 *
 * @param out an array of 1*blocks sha256 hash values
 * @param midstate the initial state (e.g. sha256_ctx.s)
 * @param in an array of 64*blocks SHA256 compression round inputs
 * @param blocks the number of parallel SHA256 hash operations to perform
 *
 * This rather specialized API is used to compute multiple SHA256 compression
 * rounds, in parallel, using the same initial state vector but different input
 * blocks.  The typical use case is for grinding hash preimage partial solutions
 * (e.g. proof-of-work), which can be done very efficiently by precomputing a
 * midstate vector and then attempting multiple final compression rounds in
 * parallel.
 *
 * For maximum performance blocks should be a multiple of 8, as that is the
 * highest degree of parallelism on any presently supported architecture.
 *
 * Note that the midstate is delivered as host-ordered unsigned integers, the
 * same as sha256_ctx.s, but the output is a standard SHA256 network-ordered
 * hash.
 *
 * Example:
 * void sha256_write_and_finalize8(struct sha256_ctx* ctx, const unsigned char nonce1[4], const unsigned char nonce2[4], const unsigned char final[4], struct sha256 hashes[8])
 * {
 *         unsigned char blocks[8*64] = { 0 };
 *         int i;
 *         for (i = 0; i < 8; ++i) {
 *                 memcpy(blocks + i*64 + 0, nonce1, 4);
 *                 memcpy(blocks + i*64 + 4, nonce2, 4);
 *                 memcpy(blocks + i*64 + 8, final, 4);
 *                 blocks[i*64 + 12] = 0x80; // padding byte
 *                 WriteBE64(blocks + i*64 + 56, (ctx->bytes + 12) << 3);
 *                 nonce2 += 4;
 *         }
 *         sha256_midstate(hashes, ctx->s, blocks, 8);
 * }
 */
void sha256_midstate(struct sha256 out[], const uint32_t midstate[8], const unsigned char in[], size_t blocks);

#ifdef __cplusplus
}
#endif

#endif /* SHA2__SHA256_H */

/* End of File
 */
