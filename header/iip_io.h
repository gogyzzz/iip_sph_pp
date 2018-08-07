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
#ifndef IIP_IO_H
#define IIP_IO_H

#include "iip_type.h"

/*** READ MATLAB .bin FILE ****/
void read_mat(const char* filename, MAT* mat);
void read_cmat(const char* filename, CMAT* mat);

/*** WRTIE MATLAB .bin FILE ****/
void write_mat(const char* filename, MAT* mat);
void write_cmat(const char* filename, CMAT* mat);
#endif
