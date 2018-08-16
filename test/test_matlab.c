/*
 * ===========================================================
 *           Copyright (c) 2018, __IIPLAB__
 *                All rights reserved.
 *
 * This Source Code Form is subject to the terms of
 * the Mozilla Public License, v. 2.0.
 * If a copy of the MPL was not distributed with this file,
 *  You can obtain one at http://mozilla.org/MPL/2.0/.
 * ===========================================================
 */

#include "mother.h"

int aks_preheat();
void preheat();

int main() {
  int selection=0;

  while(1){
    printf("============= TESTING =============\n");
    printf(" IIPLAB - project iip sph pp       \n");
    printf("-----------------------------------\n");
    printf(" 1. Run Test                       \n");
    printf(" 2. View Testing List              \n");
    printf(" 3. Exit                           \n");
    printf("===================================\n");

    scanf_s("%d", &selection);
    switch(selection){
      case 1:
        test_verification(ask_preheat(), 1, 1);
        break;
      // case 2:
      //   test_performance(ask_preheat(), 1);
      //   break;
      case 2:
        test_viewlist();
        break;
      case 3:
        return 0;
        break;
      default:
        printf("Wrong selection.\n");
        break;
    }
    printf("\n\n");
  }

  return 0;
}


int ask_preheat(){
  char temp = 0;
  printf(" Preheat CPU? This may take a few seconds.(y/n) : ");
  while(1){
    temp = getchar();
    if (temp == 'Y' || temp == 'y'){
      return 1;
    }
    else if (temp == 'N' || temp == 'n')
      return 0;
  }
}