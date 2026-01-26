#! /bin/bash 

# Compressing
tar -cf - 4_POST_PROC | zstd -19 --ultra -T0 > 4_POST_PROC_v3.tar.zst

# Uncompressing
#tar --zstd -xf 4_POST_PROC_v3.tar.zst
