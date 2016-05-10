        /* 0 */
        "pclmulqdq $0, %%xmm6, %%xmm6\n\t"

        "pxor %%xmm4, %%xmm4\n\t"
        "palignr %[M0BS], %%xmm6, %%xmm4\n\t"
        "psrlq %[M0BSL], %%xmm4\n\t"
        "pand %[K0MASK], %%xmm6\n\t"
        "pxor %%xmm4, %%xmm6\n\t"
        "psllq %[F0D], %%xmm4\n\t"
        "pxor %%xmm4, %%xmm6\n\t"

        "vpsrlq %[M0], %%xmm6, %%xmm4\n\t"
        "pand %[K0MASK], %%xmm6\n\t"
        "pxor %%xmm4, %%xmm6\n\t"
        "psllq %[F0D], %%xmm4\n\t"
        "pxor %%xmm4, %%xmm6\n\t"

        "pclmulqdq $0, %%xmm5, %%xmm6\n\t"

        "pxor %%xmm4, %%xmm4\n\t"
        "palignr %[M0BS], %%xmm6, %%xmm4\n\t"
        "psrlq %[M0BSL], %%xmm4\n\t"
        "pand %[K0MASK], %%xmm6\n\t"
        "pxor %%xmm4, %%xmm6\n\t"
        "psllq %[F0D], %%xmm4\n\t"
        "pxor %%xmm4, %%xmm6\n\t"

        "vpsrlq %[M0], %%xmm6, %%xmm4\n\t"
        "pand %[K0MASK], %%xmm6\n\t"
        "pxor %%xmm4, %%xmm6\n\t"
        "psllq %[F0D], %%xmm4\n\t"
        "pxor %%xmm4, %%xmm6\n\t"

        "pclmulqdq $0, %%xmm6, %%xmm6\n\t"

        "pxor %%xmm4, %%xmm4\n\t"
        "palignr %[M0BS], %%xmm6, %%xmm4\n\t"
        "psrlq %[M0BSL], %%xmm4\n\t"
        "pand %[K0MASK], %%xmm6\n\t"
        "pxor %%xmm4, %%xmm6\n\t"
        "psllq %[F0D], %%xmm4\n\t"
        "pxor %%xmm4, %%xmm6\n\t"

        "vpsrlq %[M0], %%xmm6, %%xmm4\n\t"
        "pand %[K0MASK], %%xmm6\n\t"
        "pxor %%xmm4, %%xmm6\n\t"
        "psllq %[F0D], %%xmm4\n\t"
        "pxor %%xmm4, %%xmm6\n\t"

