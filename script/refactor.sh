#!/bin/sh

for arg in "$@"
do
  echo "$arg"

# sed -i 's/ <func>(/ <func>_mat(/g' $arg


#############################
#### Refactoring History ####
#############################

# sed -i 's/ axpy(/ axpy_mat(/g' $arg
# sed -i 's/ caxpy(/ axpy_cmat(/g' $arg 
# sed -i 's/ asum(/ asum_mat(/g' $arg
# sed -i 's/ casum(/ asum_cmat(/g' $arg
# sed -i 's/ dot(/ dot_mat(/g' $arg
#  sed -i 's/ cdot(/ dot_cmat(/g' $arg
#  sed -i 's/ udot(/ cdot_cmat(/g' $arg
#  sed -i 's/ swap(/ swap_mat(/g' $arg
#  sed -i 's/ cswap(/ swap_cmat(/g' $arg

# sed -i 's/ amax(/ temp_amax_mat(/g' $arg
# sed -i 's/ camax(/ temp_amax_cmat(/g' $arg
# sed -i 's/ amin(/ temp_amin_mat(/g' $arg
# sed -i 's/ camin(/ temp_amin_cmat(/g' $arg

# sed -i 's/ rotg(/ rotg_mat(/g' $arg
# sed -i 's/ crotg(/ rotg_cmat(/g' $arg
# sed -i 's/ nrm2(/ nrm2_mat(/g' $arg
# sed -i 's/ cnrm2(/ nrm2_cmat(/g' $arg
# sed -i 's/ rot(/ rot_mat(/g' $arg
# sed -i 's/ crot(/ rot_cmat(/g' $arg
# sed -i 's/ scal(/ scal_mat(/g' $arg
# sed -i 's/ cscal(/ scal_cmat(/g' $arg
# sed -i 's/ uscal(/ cscal_cmat(/g' $arg
# sed -i 's/ add(/ add_mat(/g' $arg
# sed -i 's/ cadd(/ add_cmat(/g' $arg
# sed -i 's/ uadd(/ cadd_cmat(/g' $arg

# sed -i 's/ gemv(/ gemv_mat(/g' $arg
# sed -i 's/ cgemv(/ gemv_cmat(/g' $arg

# sed -i 's/ gemm(/ gemm_mat(/g' $arg
# sed -i 's/ cgemm(/ gemm_cmat(/g' $arg

####################################

# sed -i 's/alloc_MAT/alloc_mat/g' $arg
# sed -i 's/alloc_CMAT/alloc_cmat/g' $arg

# sed -i 's/mem_MAT/mpalloc_mat/g' $arg
# sed -i 's/mem_CMAT/mpalloc_cmat/g' $arg

# sed -i 's/mem_submat/mpsubmat/g' $arg
# sed -i 's/memcsubmat/mpcsubmat/g' $arg

# sed -i 's/getbydim/get_by_dim/g' $arg
# sed -i 's/cgetbydim/cget_by_dim/g' $arg

# sed -i 's/id_MAT/ident_mat/g' $arg

# sed -i 's/free_MAT/free_mat/g' $arg
# sed -i 's/free_CMAT/free_cmat/g' $arg

# sed -i 's/print_MAT/print_mat/g' $arg
# sed -i 's/print_CMAT/print_cmat/g' $arg

# sed -i 's/free_mem_MAT/mpfree_mat/g' $arg
# sed -i 's/free_mem_CMAT/mpfree_cmat/g' $arg

# sed -i 's/getFILETIMEoffset/get_filetime_offset/g' $arg
# sed -i 's/post_append/append_post/g' $arg

# sed -i 's/iip_malloc/mpalloc/g' $arg
# sed -i 's/iip_free/mpfree/g' $arg

# sed -i 's/read_WAV/read_wav/g' $arg
# sed -i 's/write_WAV/write_wav/g' $arg
# sed -i 's/WAV2MAT/wav2mat/g' $arg
# sed -i 's/WAV_BUF2MAT/wav_buf2mat/g' $arg
# sed -i 's/MAT2WAV/mat2wav/g' $arg
# sed -i 's/free_WAV/free_wav/g' $arg
# sed -i 's/print_WAV/print_wav/g' $arg

done
