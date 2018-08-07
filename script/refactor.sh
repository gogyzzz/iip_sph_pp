#!/bin/sh

for arg in "$@"
do
  echo "$arg"

sed -i 's/alloc_MAT/alloc_mat/g' $arg
sed -i 's/alloc_CMAT/alloc_cmat/g' $arg

sed -i 's/mem_MAT/mpalloc_mat/g' $arg
sed -i 's/mem_CMAT/mpalloc_cmat/g' $arg

sed -i 's/mem_submat/mpsubmat/g' $arg
sed -i 's/memcsubmat/mpcsubmat/g' $arg

sed -i 's/getbydim/get_by_dim/g' $arg
sed -i 's/cgetbydim/cget_by_dim/g' $arg

sed -i 's/id_MAT/ident_mat/g' $arg

sed -i 's/free_MAT/free_mat/g' $arg
sed -i 's/free_CMAT/free_cmat/g' $arg

sed -i 's/print_MAT/print_mat/g' $arg
sed -i 's/print_CMAT/print_cmat/g' $arg

sed -i 's/free_mem_MAT/mpfree_mat/g' $arg
sed -i 's/free_mem_CMAT/mpfree_cmat/g' $arg

sed -i 's/getFILETIMEoffset/get_filetime_offset/g' $arg
sed -i 's/post_append/append_post/g' $arg

sed -i 's/iip_malloc/mpalloc/g' $arg
sed -i 's/iip_free/mpfree/g' $arg

sed -i 's/read_WAV/read_wav/g' $arg
sed -i 's/write_WAV/write_wav/g' $arg
sed -i 's/WAV2MAT/wav2mat/g' $arg
sed -i 's/WAV_BUF2MAT/wav_buf2MAT/g' $arg
sed -i 's/MAT2WAV/mat2wav/g' $arg
sed -i 's/free_WAV/free_wav/g' $arg
sed -i 's/print_WAV/print_wav/g' $arg

done
