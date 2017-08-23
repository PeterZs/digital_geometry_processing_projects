//
//    File: misc.h
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#ifndef _MISC_H
#define _MISC_H

#include <stdio.h>
#include <math.h>
#include <wctype.h>

using namespace std;

#ifndef TRUE
   #define TRUE 1
#endif

#ifndef FALSE
   #define FALSE 0
#endif

#ifndef M_PI
   #define M_PI 3.14159265359
#endif

#define MIN(TYPE, A, B) ((TYPE) A < (TYPE) B ? A : B)
#define MAX(TYPE, A, B) ((TYPE) A > (TYPE) B ? A : B)

#define GRAD2RAD(G) (G/180.0*M_PI)
#define RAD2GRAD(R) (R/M_PI*180.0)

#define POSITIVE(A) (A > 0 ? A : 0)

#define FILE_ERROR(F, T)  (fprintf(stderr, "\n\nFile descriptor %d: %s\n\n", fileno(F), T))

#define SKIP_COMMENTS(F) \
{ \
  int c; \
\
  while ((c=fgetc(F)) == '#') \
  { \
    do \
    { \
      c=fgetc(F); \
    } \
    while (c != EOF && c != '\n'); \
  } \
\
  ungetc(c, F); \
}

#endif
