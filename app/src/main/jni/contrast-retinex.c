/* GIMP - The GNU Image Manipulation Program
 * Copyright (C) 1995 Spencer Kimball and Peter Mattis
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
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
#include "contrast-retinex.h"

#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

#define MAX_RETINEX_SCALES    8
#define MIN_GAUSSIAN_SCALE   16
#define MAX_GAUSSIAN_SCALE  250
#define SCALE_WIDTH         150
#define ENTRY_WIDTH           4

typedef struct{
    gint scale;
    gint nscales;
    gint scales_mode;
    gfloat cvar;
} RetinexParams;

#define RETINEX_UNIFORM 0
#define RETINEX_LOW     1
#define RETINEX_HIGH    2

static gfloat RetinexScales[ MAX_RETINEX_SCALES ];

typedef struct{
    gint N;
    gfloat sigma;
    gdouble B;
    gdouble b[4];
} gauss3_coefs;
//
static void retinex_scales_distribution ( gfloat *scales, gint nscales, gint mode, gint s );
static void compute_mean_var ( gfloat *src, gfloat *mean, gfloat *var, gint size, gint bytes );
static void gausssmooth ( gfloat *in, gfloat *out, gint size, gint rowtride, gauss3_coefs *c );
//
static void compute_coefs3 (gauss3_coefs *c, gfloat sigma);
static RetinexParams rvals ={
    200,             /* Scale */
    4,               /* Scales */
    RETINEX_UNIFORM, /* Echelles reparties uniformement */
    3             /* A voir */
};
// #define RETINEX_UNIFORM 0
//  #define RETINEX_LOW     1
//  #define RETINEX_HIGH    2

/*
 * This function is the heart of the algo.
 * (a)  Filterings at several scales and sumarize the results.
 * (b)  Calculation of the final values.
 */
void
MSRCR ( guchar *src, gint width, gint height, gint bytes )
{
    gint scale,row,col;
    gint i,j;
    gint size;
    gint channel;
    guchar *psrc = NULL;    /* backup pointer for src buffer */
    gfloat *dst  = NULL;    /* float buffer for algorithm */
    gfloat *pdst = NULL;    /* backup pointer for float buffer */
    gfloat *in, *out;
    gint channelsize;   /* Float memory cache for one channel */
    gfloat weight;
    gauss3_coefs coef;
    gfloat mean, var;
    gfloat mini, range, maxi;
    gfloat alpha;
    gfloat gain;
    gfloat offset;
    gdouble max_preview = 0.0;
    /*
            Allocate all the memory needed for algorithm
        */
    size = width * height * bytes;
    dst = malloc ( size * sizeof (gfloat) );
    if (dst == NULL){
        //  g_warning ("Failed to allocate memory");
        return;
    }
    memset ( dst, 0, size * sizeof (gfloat) );
    channelsize  = ( width * height );
    in  = (gfloat *) malloc (channelsize * sizeof (gfloat));
    if (in == NULL){
        free (dst);
        // g_warning ("Failed to allocate memory");
        return; /* do some clever stuff */
    }
    out  = (gfloat *) malloc (channelsize * sizeof (gfloat));
    if (out == NULL){
        free (in);
        free (dst);
        //  g_warning ("Failed to allocate memory");
        return; /* do some clever stuff */
    }
    /*
            Calculate the scales of filtering according to the
            number of filter and their distribution.
        */
    retinex_scales_distribution (
        RetinexScales, rvals.nscales, rvals.scales_mode, rvals.scale );
    /*
            Filtering according to the various scales.
            Summerize the results of the various filters according to a
            specific weight(here equivalent for all).
        */
    weight = 1./ (gfloat) rvals.nscales;
    /*
            The recursive filtering algorithm needs different coefficients according
            to the selected scale (~ = standard deviation of Gaussian).
        */
    for (channel = 0; channel < 3; channel++){
        gint pos;
        for (i = 0, pos = channel; i < channelsize ; i++, pos += bytes){
            /* 0-255 => 1-256 */
            in[i] = (gfloat)(src[pos] + 1.0);
        }
        for (scale = 0; scale < rvals.nscales; scale++){
            compute_coefs3 (&coef, RetinexScales[scale]);
            /*
                        *  Filtering (smoothing) Gaussian recursive.
                        *
                        *  Filter rows first
                        */
            for (row=0 ;row < height; row++){
              pos =  row * width;
              gausssmooth (in + pos, out + pos, width, 1, &coef);
            }
            memcpy(in,  out, channelsize * sizeof(gfloat));
            memset(out, 0  , channelsize * sizeof(gfloat));
            /*
                        *  Filtering (smoothing) Gaussian recursive.
                        *
                        *  Second columns
                        */
            for (col=0; col < width; col++){
                pos = col;
                gausssmooth(in + pos, out + pos, height, width, &coef);
            }
            /*
                            Summarize the filtered values.
                            In fact one calculates a ratio between the original values and the filtered values.
                         */
            for (i = 0, pos = channel; i < channelsize; i++, pos += bytes){
                dst[pos] += weight * (log (src[pos] + 1.) - log (out[i]));
            }
        }
    }
    free(in);
    free(out);
    /*
            Final calculation with original value and cumulated filter values.
            The parameters gain, alpha and offset are constants.
        */
    /* Ci(x,y)=log[a Ii(x,y)]-log[ Ei=1-s Ii(x,y)] */
    alpha  = 128.;
    gain   = 1.;
    offset = 0.;
    for (i = 0; i < size; i += bytes){
        gfloat logl;
        psrc = src+i;
        pdst = dst+i;
        logl = log((gfloat)psrc[0] + (gfloat)psrc[1] + (gfloat)psrc[2] + 3.);
        pdst[0] = gain * ((log(alpha * (psrc[0]+1.)) - logl) * pdst[0]) + offset;
        pdst[1] = gain * ((log(alpha * (psrc[1]+1.)) - logl) * pdst[1]) + offset;
        pdst[2] = gain * ((log(alpha * (psrc[2]+1.)) - logl) * pdst[2]) + offset;
    }
    /*
            Adapt the dynamics of the colors according to the statistics of the first and second order.
            The use of the variance makes it possible to control the degree of saturation of the colors.
        */
    pdst = dst;
    compute_mean_var (pdst, &mean, &var, size, bytes);
    mini = mean - rvals.cvar*var;
    maxi = mean + rvals.cvar*var;
    range = maxi - mini;
    if (!range) range = 1.0;
    for (i = 0; i < size; i+= bytes){
        psrc = src + i;
        pdst = dst + i;
        for (j = 0 ; j < 3 ; j++){
            gfloat c = 255 * ( pdst[j] - mini ) / range;
            psrc[j] = (guchar) CLAMP (c, 0, 255);
        }
    }
    free (dst);
}
/*
 * calculate scale values for desired distribution.
 */
static void
retinex_scales_distribution(gfloat* scales, gint nscales, gint mode, gint s)
{
    if (nscales == 1){
    /* For one filter we choose the median scale */
        scales[0] = (gint) s / 2;
    }
    else if (nscales == 2){
     /* For two filters whe choose the median and maximum scale */
      scales[0] = (gint) s / 2;
      scales[1] = (gint) s;
    }
    else{
        gfloat size_step = (gfloat) s / (gfloat) nscales;
        gint i;
        switch(mode){

            case RETINEX_UNIFORM:
                for(i = 0; i < nscales; ++i)
                scales[i] = 2. + (gfloat) i * size_step;
            break;

            case RETINEX_LOW:
                size_step = (gfloat) log(s - 2.0) / (gfloat) nscales;
                for (i = 0; i < nscales; ++i)
                scales[i] = 2. + pow (10, (i * size_step) / log (10));
            break;

            case RETINEX_HIGH:
                size_step = (gfloat) log(s - 2.0) / (gfloat) nscales;
                for (i = 0; i < nscales; ++i)
                scales[i] = s - pow (10, (i * size_step) / log (10));
            break;

            default:
            break;
        }
    }
}

/*
 * Calculate the coefficients for the recursive filter algorithm
 * Fast Computation of gaussian blurring.
 */
static void
compute_coefs3 (gauss3_coefs *c, gfloat sigma)
{
    /*
        * Papers:  "Recursive Implementation of the gaussian filter.",
        *          Ian T. Young , Lucas J. Van Vliet, Signal Processing 44, Elsevier 1995.
        * formula: 11b       computation of q
        *          8c        computation of b0..b1
        *          10        alpha is normalization constant B
        */
    gfloat q, q2, q3;

    if (sigma >= 2.5){
      q = 0.98711 * sigma - 0.96330;
    }
    else if ((sigma >= 0.5) && (sigma < 2.5)){
      q = 3.97156 - 4.14554 * (gfloat) sqrt ((double) 1 - 0.26891 * sigma);
    }
    else{
      q = 0.1147705018520355224609375;
    }
    q2 = q * q;
    q3 = q * q2;
    c->b[0] = (1.57825+(2.44413*q)+(1.4281 *q2)+(0.422205*q3));
    c->b[1] = (        (2.44413*q)+(2.85619*q2)+(1.26661 *q3));
    c->b[2] = (                   -((1.4281*q2)+(1.26661 *q3)));
    c->b[3] = (                                 (0.422205*q3));
    c->B = 1.0-((c->b[1]+c->b[2]+c->b[3])/c->b[0]);
    c->sigma = sigma;
    c->N = 3;
}
/*
  Definit comment sont repartis les
  differents filtres en fonction de
  l'echelle (~= ecart type de la gaussienne)
 */
static void
gausssmooth (gfloat *in, gfloat *out, gint size, gint rowstride, gauss3_coefs *c){
  /*
   * Papers:  "Recursive Implementation of the gaussian filter.",
   *          Ian T. Young , Lucas J. Van Vliet, Signal Processing 44, Elsevier 1995.
   * formula: 9a        forward filter
   *          9b        backward filter
   *          fig7      algorithm
   */
    gint i,n, bufsize;
    gfloat *w1,*w2;
    /* forward pass */
    bufsize = size+3;
    size -= 1;
    w1 = (gfloat *) malloc (bufsize * sizeof (gfloat));
    w2 = (gfloat *) malloc (bufsize * sizeof (gfloat));
    w1[0] = in[0];
    w1[1] = in[0];
    w1[2] = in[0];
    for ( i = 0 , n=3; i <= size ; i++, n++){
        w1[n] = (gfloat)(c->B*in[i*rowstride] +
                       ((c->b[1]*w1[n-1] +
                         c->b[2]*w1[n-2] +
                         c->b[3]*w1[n-3] ) / c->b[0]));
    }

    /* backward pass */
    w2[size+1]= w1[size+3];
    w2[size+2]= w1[size+3];
    w2[size+3]= w1[size+3];
    for (i = size, n = i; i >= 0; i--, n--){
        w2[n]= out[i * rowstride] = (gfloat)(c->B*w1[n+3] +
                                           ((c->b[1]*w2[n+1] +
                                             c->b[2]*w2[n+2] +
                                             c->b[3]*w2[n+3] ) / c->b[0]));
    }
    free (w1);
    free (w2);
}
/*
 * Calculate the average and variance in one go.
 */
static void
compute_mean_var (gfloat *src, gfloat *mean, gfloat *var, gint size, gint bytes){
    gfloat vsquared;
    gint i,j;
    gfloat *psrc;
    vsquared = 0;
    *mean = 0;
    for (i = 0; i < size; i+= bytes){
        psrc = src+i;
        for (j = 0 ; j < 3 ; j++){
            *mean += psrc[j];
            vsquared += psrc[j] * psrc[j];
        }
    }
    *mean /= (gfloat) size; /* mean */
    vsquared /= (gfloat) size; /* mean (x^2) */
    *var = ( vsquared - (*mean * *mean) );
    *var = sqrt(*var); /* var */
}
