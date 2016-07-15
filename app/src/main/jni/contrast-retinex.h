#ifndef CONTRASTRETINEX
#define CONTRASTRETINEX

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
//
#include <opencv2/core/core_c.h>
#include <opencv2/imgproc/imgproc_c.h>


typedef char gchar;
typedef short gshort;
typedef long glong;
typedef int gint;
typedef gint gboolean;
typedef unsigned char guchar;
typedef unsigned short gushort;
typedef unsigned long gulong;
typedef unsigned int guint;
typedef float gfloat;
typedef double gdouble;


//==========================================================
// src = RGB  image, gint bytes = 3
// src = RGBA image, gint bytes = 4
//----------------------------------------------------------
// MSRCP = MultiScale Retinex with Chromaticity Preservation
//
void MSRCP( guchar *src, gint width, gint height, gint bytes, gfloat cvar );
//----------------------------------------------------------
// GLAY_BLUR
//
void GLAY_BLUR( guchar *src, gint width, gint height, gint bytes, gfloat segma );
//==========================================================
#ifdef __cplusplus
}
#endif

#endif // CONTRASTRETINEX

