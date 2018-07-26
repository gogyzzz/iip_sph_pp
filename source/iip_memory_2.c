#include "iip_type.h"

#define MAX_MEM_PAGE 20
#define MEM_PAGE_BASIC_SIZE 256
#define MAX_MEM_BLOCK 16

typedef struct mem_node {
  void *p;
  UINT usedl struct mem_node *next;
} mem_node;

typedef struct memory_list {
  UINT alloced;
  UINT used;
  uint64_t block_size;
  mem_node *front;
} memory_list;

static memory_list *mem_list

    void
    init() {
  ITER i;

  //크게 만들어서 코스트가 크기 않으므로 MAX_MEM_PAGE 개 할당
  // MEM_PAGE_BASIC_SIZE * 2^MAX_MEM_PAGE 크기까지 커버
  mem_list = (memory_list *)malloc(sizeof(memory_list) * MAX_MEM_PAGE);
  printf("list     | alloc | used | block_size\n");
  for (i = 0; i < MAX_MEM_PAGE; i++) {
    mem_list[i].alloced = 0;
    mem_list[i].used = 0;
    mem_list[i].front = NULL;
    mem_list[i].block_size = MEM_PAGE_BASIC_SIZE << i;
    printf("list[%2d] |  %4u | %4u | %9lu\n", i, mem_list[i].alloced,
           mem_list[i].used, mem_list[i].block_size);
  }
}

void *iip_malloc(UINT size) {
  UINT list_idx, temp_num;
  mem_node *temp_node;
  mem_node *cur_node;
  // get proper index
  temp_num = size - 1;

  if (temp_num < MEM_PAGE_BASIC_SIZE)
    list_idx = 0;
  else {
    /* 4096 : -12
     * 256  : - 8
     *
     *
     */
    list_idx = -8 + 1;
    while (temp_num != 1) {
      temp_num = temp_num >> 1;
      list_idx++;
    }
  }
  // printf("malloc %u at list[%u] alloced : %u | used :
  // %u\n",size,list_idx,mem_list[list_idx].alloced,mem_list[list_idx].used);
  //처음임
  if (mem_list[list_idx].front == NULL) {
    printf("create front for list[%d] : %u\n", list_idx, size);
    temp_node = (mem_node *)malloc(sizeof(mem_node));
    temp_node->p = (void *)malloc(mem_list[list_idx].block_size);
    temp_node->used = 1;
    temp_node->next = NULL;

    mem_list[list_idx].used++;
    mem_list[list_idx].alloced++;
    mem_list[list_idx].front = temp_node;
    show_list();
    return mem_list[list_idx].front->p;
  } else {
    //자리 없음
    if (mem_list[list_idx].alloced == mem_list[list_idx].used) {
      printf("create node for list[%d] : %u\n", list_idx, size);
      temp_node = (mem_node *)malloc(sizeof(mem_node));
      temp_node->p = (void *)malloc(mem_list[list_idx].block_size);
      temp_node->used = 1;
      temp_node->next = NULL;
      //	printf("point %u\n",(UINT)temp_node->p);
      cur_node = mem_list[list_idx].front;
      while ((cur_node->next) != NULL) cur_node = cur_node->next;
      cur_node->next = temp_node;
      mem_list[list_idx].alloced++;
      mem_list[list_idx].used++;
      show_list();
      return temp_node->p;
    }
    //자리 있음
    else if (mem_list[list_idx].alloced > mem_list[list_idx].used) {
      printf("get empty node of list[%d] : %u\n", list_idx, size);
      cur_node = mem_list[list_idx].front;
      while (cur_node->used == 1) cur_node = cur_node->next;
      mem_list[list_idx].used++;
      cur_node->used = 1;
      show_list();
      return cur_node->p;
    } else {
      printf("ERROR : memory list invalid\n");
      exit(-1);
    }
  }
}

void iip_free(void *p) {
  mem_node *cur_node;
  mem_node *bef_node;
  mem_node *temp_node;
  ITER i, j;
  UINT found = 0;

  //	printf("freeing %u\n",(UINT)p);
  for (i = 0; i < MAX_MEM_PAGE; i++) {
    // printf("searching list[%d] aloc : %u used :
    // %u\n",i,mem_list[i].alloced,mem_list[i].used);
    cur_node = mem_list[i].front;
    if (cur_node == NULL)
      continue;
    else {
      //첫 노드 해제해야함
      if (cur_node->p == p) {
        //	printf("first %u %u\n",(UINT)cur_node->p,p);
        temp_node = cur_node;
        temp_node->used = 0;
        if (cur_node->next != NULL) {
          mem_list[i].front = cur_node->next;
          // cur_node를 맨 뒤로
          while (cur_node != NULL) {
            bef_node = cur_node;
            cur_node = cur_node->next;
          }
          // bef_node 는 마지막 노드
          bef_node->next = temp_node;
          temp_node->next = NULL;
        }
        //첫노드가 마지막노드
        else {
        }
        mem_list[i].used--;
        found = 1;
        break;
      } else {
        bef_node = cur_node;
        cur_node = cur_node->next;
        while ((cur_node != NULL) && (cur_node->p != p)) {
          bef_node = cur_node;
          cur_node = cur_node->next;
        }
        if (cur_node == NULL) continue;
        //	printf("found %u\n",(UINT)cur_node->p);
        // cur_node 를 해제함
        temp_node = cur_node;
        temp_node->used = 0;

        //마지막 노드가 아님
        if (temp_node->next != NULL) {
          bef_node->next = cur_node->next;

          //마지막 노드로
          while (cur_node != NULL) {
            bef_node = cur_node;
            cur_node = cur_node->next;
          }
          bef_node->next = temp_node;
          temp_node->next = NULL;
        }
        //이미 마지막 노드라면
        else {
        }
        mem_list[i].used--;
        found = 1;
        break;
      }
    }
  }

  printf("freed list[%d] \n", i);
  show_list();

  if (found == 0) printf("WARNING : free-operaton wasn't done completly\n");
}

void finit() {
  ITER i;
  mem_node *cur;
  mem_node *temp;

  for (i = 0; i < MAX_MEM_PAGE; i++) {
    cur = mem_list[i].front;
    while (cur != NULL) {
      temp = cur;
      cur = cur->next;
      free(temp);
    }
  }
  free(mem_list);
}

void show_list() {
  ITER i;
  UINT cnt;
  mem_node *cur;
  mem_node *temp;

  for (i = 0; i < MAX_MEM_PAGE; i++) {
    cur = mem_list[i].front;
    cnt = 0;
    if (mem_list[i].alloced == 0) continue;
    printf("== list[%d] == alloced : %u | used : %u\n", i, mem_list[i].alloced,
           mem_list[i].used);
    while (cur != NULL) {
      printf("%u : %u(%u) \n", cnt, (UINT)(cur->p), cur->used);
      temp = cur;
      cur = cur->next;
    }
    printf("\n");
  }
}
