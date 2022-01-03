/* memwrapper.c
*
* Copyright (C) <2016>  Giuseppe Marco Randazzo
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "memwrapper.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#endif


void *xmalloc(size_t size)
{
  void *ptr = NULL;
  ptr = malloc(size);
  if(ptr == NULL){
    fprintf(stderr, "Memory Exhausted!\n");
    abort();
  }
  return ptr;
}

void *xrealloc(void *ptr, size_t size)
{
  register void *value = realloc (ptr, size);
  if (value == 0){
    fprintf(stderr, "Memory Exhausted!\n");
    abort();
  }
  return value;
}

void xfree(void *ptr)
{
  free(ptr);
}
