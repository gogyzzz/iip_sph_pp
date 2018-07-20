#include "iip_type.h"


/*****************************
 **** MEMORY MANAGER *********
 *****************************/

void init() {
	int i;

#if USE_CUDA
	cudaDeviceProp prop;
	cublasCreate(&handle);
	cudaGetDeviceProperties(&prop,0);
	max_thread = prop.maxThreadsDim[0];
	max_block = prop.maxGridSize[1];
#endif	
	
	for (i = 0; i < MAX_MEM_PAGE; i++) {
		memory_pool[i] = NULL;
		memory_log[i] = NULL;
		log_cnt[i] = 0;
	}

	memory_pool[0] = (void*)malloc(MEM_PAGE_BASIC_SIZE);
	memory_log[0] = (USED_MEM*)malloc(sizeof(USED_MEM) * (MEM_PAGE_BASIC_SIZE / LOG_ALLOC_UNIT));
	memset(memory_log[0], 0, sizeof(USED_MEM) * (MEM_PAGE_BASIC_SIZE / LOG_ALLOC_UNIT));
	pool_cnt = 1;


}

void finit() {
	int i;

#if USE_CUDA
	cublasDestory(handle);
#endif

	for (i = 0; i < MAX_MEM_PAGE; i++) {
		if (memory_pool[i] != NULL) {
			free(memory_pool[i]);
		}
		if (memory_log[i] != NULL) {
			free(memory_log[i]);
		}
		log_cnt[i] = 0;
	}
}

void* iip_malloc(unsigned long int size) {
	signed long int alloc_idx = 0;
	int i;
	unsigned int page_size = 0;

	for (i = 0; i < pool_cnt; i++) {
		alloc_idx = page_alloc_isable(i, size);
		
		#if DEBUG
		printf("alloc_idx(%d) = page_alloc_isable(%lu,%lu)\n",alloc_idx,i,size);
		#endif
		// Allocable, return allocated address.
		if (alloc_idx != -1) {
			memory_log[i][log_cnt[i]].size = size;
			memory_log[i][log_cnt[i]].frag_idx = alloc_idx;
			log_cnt[i]++;
#if DEBUG
			printf("memory_log[%d][%lu=log_cnt[%d]]\n",i,log_cnt[i],i);
			printf("return %lu[%lu] + %ld\n",(unsigned long int)memory_pool[i],i,alloc_idx);
#endif
			return (unsigned long int)memory_pool[i] + alloc_idx;
		}
	}

	do {
		// No more space to alloc, call malloc()..

		page_size = MEM_PAGE_BASIC_SIZE;
		for (i = 0; i < pool_cnt; i++) {
			page_size *= 2;
		}
		memory_pool[pool_cnt] = (void*)malloc(page_size);
		memory_log[pool_cnt] = (USED_MEM*)malloc(sizeof(USED_MEM) * (page_size / LOG_ALLOC_UNIT));
		memset(memory_log[pool_cnt], 0, sizeof(USED_MEM) * (page_size / LOG_ALLOC_UNIT));
		pool_cnt++;
	} while (page_alloc_isable(pool_cnt-1, size) == -1);

	// Allocate
	if (memory_pool[pool_cnt-1] != NULL) {
		memory_log[pool_cnt-1][0].size = size;
		memory_log[pool_cnt-1][0].frag_idx = 0;
		log_cnt[pool_cnt-1]++;
#if DEBUG
			printf("return %lu[%lu]\n",(unsigned long int)memory_pool[pool_cnt-1],pool_cnt-1);
#endif
		return memory_pool[pool_cnt-1];
	}

	printf("Failed to Allocate!\n");
	printf(" Allocation Size : %ld\n", size);

	exit(0);

	return NULL;
}

void iip_free(void *ptr) {
	int i,j,k;
	for (i = 0; i < pool_cnt; i++) {
		for (j = 0; j < log_cnt[i]; j++) {
			// Find!
			if (ptr == (unsigned long int)memory_pool[i] + memory_log[i][j].frag_idx) {
				memory_log[i][j].size = 0;
				memory_log[i][j].frag_idx = 0;
				if (log_cnt[i] != j + 1){
					for (k = j; k < log_cnt[i]-1; k++) {
						memory_log[i][k].size = memory_log[i][k + 1].size;
						memory_log[i][k].frag_idx = memory_log[i][k + 1].frag_idx;
					}
				}
				log_cnt[i]--;

				return;
			}
		}
	}
	printf("Memory redundancy release has been detected.\n");
	exit(0);
	return;
}

signed long int page_alloc_isable(int page_idx, unsigned long int require_size) {
	unsigned int max_remain = 0;
	unsigned int max_remain_start_addr = 0;
	unsigned int page_size = 0;
	unsigned int addr = 0;
	int i;
	page_size = MEM_PAGE_BASIC_SIZE;
	for (i = 0; i < page_idx; i++) {
		page_size *= 2;
	}

	// 지정된 페이지에 이미 할당된 메모리 조각들이 존재.
	
	if (log_cnt[page_idx] > 0) {
		max_remain = page_size;
		max_remain -= memory_log[page_idx][log_cnt[page_idx]-1].frag_idx;
		max_remain -= memory_log[page_idx][log_cnt[page_idx]-1].size;
		addr = memory_log[page_idx][log_cnt[page_idx]-1].frag_idx + memory_log[page_idx][log_cnt[page_idx]-1].size;
	}
	// 지정된 페이지에 할당된 메모리 조각이 없음. 페이지 전체 사용가능.
	else {
		max_remain = page_size;
		addr = 0;
	}

	printf("max_remain %u\n",max_remain);
	if (max_remain >= require_size)
	{
			return addr;
	}else
		return -1;
}
