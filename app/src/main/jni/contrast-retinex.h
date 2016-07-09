#ifndef CONTRASTRETINEX
#define CONTRASTRETINEX

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*  Provide type definitions for commonly used types.
 *  These are useful because a "gint8" can be adjusted
 *  to be 1 byte (8 bits) on all platforms. Similarly and
 *  more importantly, "gint32" can be adjusted to be
 *  4 bytes (32 bits) on all platforms.
 */
typedef char   gchar;
typedef short  gshort;
typedef long   glong;
typedef int    gint;
typedef gint   gboolean;
typedef unsigned char   guchar;
typedef unsigned short  gushort;
typedef unsigned long   gulong;
typedef unsigned int    guint;
typedef float   gfloat;
typedef double  gdouble;
/*
 * MSRCR = MultiScale Retinex with Color Restoration
 */
void MSRCR (
    guchar *src,
    gint width,
    gint height,
    gint bytes
);
#ifdef __cplusplus
}
#endif

#endif // CONTRASTRETINEX

