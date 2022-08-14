/* Copyright (c) 2022 Mark Friedenbach
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <sha2/sha256.h>

#include <stdio.h>

int main(int argc, char *argv[])
{
        int exit_code = 0;
        struct sha256_ctx ctx = SHA256_INIT;
        struct sha256 hash = { 0 };
        printf("Using SHA256 algorithm: %s\n", sha256_auto_detect());
        sha256_update(&ctx, "a", 1);
        sha256_done(&hash, &ctx);
        return exit_code;
}

/* End of File
 */
