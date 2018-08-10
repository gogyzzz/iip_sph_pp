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

int main() {
  int selection=0;

  while(1){
    printf("============= TESTING =============\n");
    printf(" IIPLAB - project iip sph pp       \n");
    printf("-----------------------------------\n");
    printf(" 1. Run Test - Verification        \n");
    printf(" 2. Run Test - Performance         \n");
    printf(" 3. View Testing List              \n");
    printf(" 4. Exit                           \n");
    printf("===================================\n");

    scanf_s("%d", &selection);
    switch(selection){
      case 1:
        test_verification();
        break;
      case 2:
        test_performance();
        break;
      case 3:
        test_viewlist();
        break;
      case 4:
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
