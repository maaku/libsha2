/* Copyright (c) 2022 Mark Friedenbach
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <gtest/gtest.h>

#include <sha2/sha256.h>

TEST(gtest, assert_eq)
{
        ASSERT_EQ(0, 0);
        ASSERT_NE(0, 1);
}

TEST(sha2, self_test)
{
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

        sha256_auto_detect();

        /* Test sha256_init() */
        {
                struct sha256_ctx ctx = SHA256_INIT;
                ASSERT_EQ(memcmp(ctx.s, init, 32), 0);
        }
        {
                struct sha256_ctx ctx;
                sha256_init(&ctx);
                ASSERT_EQ(memcmp(ctx.s, init, 32), 0);
        }

        /* Test transform() for 0 through 8 transformations. */
        for (int i = 0; i <= 8; ++i) {
                struct sha256_ctx ctx = SHA256_INIT;
                sha256_update(&ctx, data + 1, i * 64);
                ASSERT_EQ(memcmp(ctx.s, result[i], 32), 0);
                ASSERT_EQ(ctx.bytes, i * 64);
        }

        /* Test transform_d64 */
        {
                struct sha256 out[1];
                sha256_double64(out, data_d64, 1);
                ASSERT_EQ(memcmp(out, result_d64, 32), 0);
        }

        /* Test transform_d64_2way, if available. */
        {
                struct sha256 out[2];
                sha256_double64(out, data_d64, 2);
                ASSERT_EQ(memcmp(out, result_d64, 64), 0);
        }

        /* Test transform_d64_4way, if available. */
        {
                struct sha256 out[4];
                sha256_double64(out, data_d64, 4);
                ASSERT_EQ(memcmp(out, result_d64, 128), 0);
        }

        /* Test transform_d64_8way, if available. */
        {
                struct sha256 out[8];
                sha256_double64(out, data_d64, 8);
                ASSERT_EQ(memcmp(out, result_d64, 256), 0);
        }
}

int main(int argc, char **argv)
{
        ::testing::InitGoogleTest(&argc, argv);

        return RUN_ALL_TESTS();
}

/* End of File
 */
