        /* 0 */
        "pclmulqdq $0, %%xmm5, %%xmm5\n\t"

        "pxor %%xmm4, %%xmm4\n\t"
        "palignr %[M0BS], %%xmm5, %%xmm4\n\t"
        "psrlq %[M0BSL], %%xmm4\n\t"
        "pand %[K0MASK], %%xmm5\n\t"
        "pxor %%xmm4, %%xmm5\n\t"
        "psllq %[F0D], %%xmm4\n\t"
        "pxor %%xmm4, %%xmm5\n\t"

        "vpsrlq %[M0], %%xmm5, %%xmm4\n\t"
        "pand %[K0MASK], %%xmm5\n\t"
        "pxor %%xmm4, %%xmm5\n\t"
        "psllq %[F0D], %%xmm4\n\t"
        "pxor %%xmm4, %%xmm5\n\t"

