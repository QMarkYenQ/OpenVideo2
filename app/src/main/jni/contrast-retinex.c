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
typedef struct{ gint scale; gint nscales; gint scales_mode; gfloat cvar;} RetinexParams;
enum{ RETINEX_UNIFORM, RETINEX_LOW, RETINEX_HIGH }; // scales_mode
//
typedef struct{ gfloat *src_float; gfloat *src_gray; gfloat *src_smooth; gfloat *src_scale; gfloat *src_temp; } image_pointer;
static gboolean allocate_image_memory( gint channelsize, image_pointer *img_g );
static void deallocate_image_memory( image_pointer *img_g );


//
static void float_channel( guchar *src, gfloat *src_float, gfloat *src_gray, gint channelsize , gint bytes );
static void log_channel( gfloat * src_gray, gfloat * src_scale,  gint channelsize );
static void log_channel_lookup( gfloat * src_gray, gfloat * src_scale, gfloat* lookup_table_log, gint channelsize );
//
static void retinex_scales_distribution ( gfloat *scales, gint nscales, gint mode, gint s );
//gaussian_smoothing
typedef struct{ gint N; gfloat sigma; gdouble B; gdouble b[4];} gauss3_coefs;
static void gaussian_smoothing_2D( gfloat *src, gfloat *dst, gfloat sigma, gint width, gint height, gfloat *temp );
static void gaussian_smoothing_1D ( gfloat *src, gfloat *dst, gint size, gint rowtride, gauss3_coefs *c );
static void compute_coefs3 ( gauss3_coefs *c, gfloat sigma);
//
static void accumulation_minus_log(  gfloat * src_scale, gfloat *src_smooth , gfloat weight, gint channelsize  );
static void accumulation_minus_log_lookup(  gfloat * src_scale, gfloat * src_smooth, gfloat weight,gfloat* lookup_table_log, gint channelsize );
// simplest_color_balance
static void simplest_color_balance( gfloat *src_gray, gfloat mean, gfloat var, gfloat cvar, gint channelsize );
static void compute_mean_var(gfloat *src, gfloat *mean, gfloat *var, gint channelsize );
//
static void compute_color_from_grayscale( gfloat *src, gfloat *src_gray, gfloat *src1_gray, gint channelsize );
static void uchar_color( gfloat * src_float, guchar * dst,  gint channelsize,  gint bytes );
//==========================================================
//
// MSRCP = MultiScale Retinex with Chromaticity Preservation
void MSRCP( guchar *src, gint width, gint height, gint bytes, gfloat cvar  )
{
    gint j;
    gint size, channelsize;
    /* Float memory cache for one channel */
    gfloat weight;
    gfloat mean, var;
    gboolean successAllocate;
    image_pointer img_g = { NULL, NULL, NULL, NULL, NULL };
    RetinexParams rvals = { 256, 3, RETINEX_UNIFORM, cvar };
    gfloat RetinexScales[ 25 ];

    gfloat lookup_table_log[ 256 ];
    clock_t start, end1;
    double diff;
    IplImage *src_ipl, *dst_ipl;

    /*
        lookup_table
     */
    lookup_table_log[0] = 0;
    for ( j = 1; j < 256; j++ ) lookup_table_log[j] = (gfloat) log(j);



    /*
        Allocate Memory
     */
    size = width * height * bytes;
    channelsize  =  width * height ;
    //
    successAllocate = allocate_image_memory(  channelsize, &img_g );
    if(!successAllocate) return;
    //
    src_ipl = cvCreateImageHeader	(	cvSize(width ,height  ), IPL_DEPTH_32F, 1 );
    src_ipl->imageData = (char *)img_g.src_gray;
    dst_ipl = cvCreateImageHeader	(  cvSize(width ,  height),  IPL_DEPTH_32F, 1 );
    dst_ipl->imageData = (char *)img_g.src_smooth;



    /*
        Image Process
     */
    start = clock(); // Start Record the time
    //
    float_channel( src, img_g.src_float, img_g.src_gray, channelsize, bytes );
    end1 = clock();
    diff = end1 - start; // ms
    printf("%f  ms - float_channel\n" , diff);

    //log_channel( img_g.src_gray, img_g.src_scale, channelsize );
    log_channel_lookup(img_g.src_gray, img_g.src_scale,lookup_table_log, channelsize );

    end1 = clock();
    diff = end1 - start; // ms
    printf("%f  ms - log_channel\n" , diff);
    //
    retinex_scales_distribution( RetinexScales, rvals.nscales, rvals.scales_mode, rvals.scale );
    //= fixed
    rvals.nscales = 3;
    RetinexScales[0] = 15.;
    RetinexScales[1] = 80.;
    RetinexScales[2] = 250.;
    // 15 80 250
    //for (j = 0; j < rvals.nscales; j++) printf("%f\n", (float)RetinexScales[j]);
    //
    weight = 1.F/ (gfloat) rvals.nscales;

    for ( j = 0; j < rvals.nscales; j++ ){
       // gaussian_smoothing_2D( img_g.src_gray, img_g.src_smooth, RetinexScales[j], width, height, img_g.src_temp );
        cvSmooth( src_ipl, dst_ipl, CV_BLUR ,(int)RetinexScales[j],(int)RetinexScales[j],0,0);

        end1 = clock();
        diff = end1 - start; // ms
        printf("%f  ms - gaussian %d\n" , diff,j);

       //  accumulation_minus_log( img_g.src_scale, img_g.src_smooth, weight, channelsize );
        accumulation_minus_log_lookup( img_g.src_scale, img_g.src_smooth, weight,lookup_table_log, channelsize );

        end1 = clock();
        diff = end1 - start; // ms
        printf("%f  ms - accumulation %d\n" , diff,j);
    }
    //
    compute_mean_var( img_g.src_scale, &mean, &var, channelsize );
    //
    end1 = clock();
    diff = end1 - start; // ms
    printf("%f  ms - compute_mean_var\n" , diff);
    //
    simplest_color_balance( img_g.src_scale, mean, var, rvals.cvar, channelsize );
    //
    end1 = clock();
    diff = end1 - start; // ms
    printf("%f  ms - simplest_color_balance\n" , diff);
    //
    compute_color_from_grayscale( img_g.src_float, img_g.src_gray, img_g.src_scale, channelsize );
    //
    end1 = clock();
    diff = end1 - start; // ms
    printf("%f  ms - compute_color_from_grayscale\n" , diff);
    //
    uchar_color( img_g.src_float, src, channelsize, bytes );
    //
    end1 = clock();
    diff = end1 - start; // ms
    printf("%f  ms - uchar_color\n" , diff);
    //
    deallocate_image_memory( &img_g );
}

// GLAY_BLUR
void GLAY_BLUR( guchar *src, gint width, gint height, gint bytes, gfloat segma )
{
    gint  i,j,pos;
    gint size, channelsize;
    gboolean successAllocate;
    image_pointer img_g = { NULL, NULL, NULL, NULL, NULL };
    /*
        Allocate Memory
     */
    size = width * height * bytes;
    channelsize  =  width * height ;
    //
    successAllocate = allocate_image_memory( channelsize, &img_g );
    if(!successAllocate) return;
    /*
        Image Process
     */
    float_channel( src, img_g.src_float, img_g.src_gray, channelsize, bytes );
    gaussian_smoothing_2D( img_g.src_gray, img_g.src_smooth, segma, width, height, img_g.src_temp );
    for ( i = 0, pos = 0; i < channelsize ; i++, pos += 3 ){
        for ( j = 0 ; j < 3 ; j++){
            img_g.src_float[pos +j] = img_g.src_smooth[i];
        }
    }
    uchar_color( img_g.src_float, src, channelsize, bytes );
    deallocate_image_memory( &img_g );
}
//
//
//==========================================================
//
static gboolean allocate_image_memory( gint channelsize, image_pointer *img_g )
{
    img_g->src_float  = malloc ( 3 * channelsize * sizeof (gfloat) );
    if (img_g->src_float == NULL){
        //  g_warning ("Failed to allocate memory");
        return 0;
    }
    //
    img_g->src_gray = malloc ( channelsize * sizeof (gfloat) );
    if (img_g->src_gray == NULL){
        free (img_g->src_float);
        //  g_warning ("Failed to allocate memory");
         return 0;
    }
    //
    img_g->src_scale = malloc ( channelsize * sizeof (gfloat) );
    if (img_g->src_scale == NULL){
        free (img_g->src_gray);
        free (img_g->src_float);
        //  g_warning ("Failed to allocate memory");
         return 0;
    }
    //
    img_g->src_temp  = (gfloat *) malloc (channelsize * sizeof (gfloat));
    if (img_g->src_temp == NULL){
        free (img_g->src_scale);
        free (img_g->src_gray);
        free (img_g->src_float);
        // g_warning ("Failed to allocate memory");
         return 0; /* do some clever stuff */
    }
    //
    img_g->src_smooth  = (gfloat *) malloc (channelsize * sizeof (gfloat));
    if (img_g->src_smooth == NULL){
        free (img_g->src_temp);
        free (img_g->src_scale);
        free (img_g->src_gray);
        free (img_g->src_float);
        //  g_warning ("Failed to allocate memory");
         return 0;/* do some clever stuff */
    }
    return 1;

    //IplImage* img = cvCreateImage( cvSize( height, width ), IPL_DEPTH_32F, 3 );
    //uchar* ptr = (uchar*) img->imageData;
    // img->imageData
    //

    // IplImage*img = cvCreateImageHeader	(	cvSize( height, width ), IPL_DEPTH_32F, 3 )
    //  img->imageSize = sizeof(arr);
    //  img->imageData = (char *)arr;

    //cvCreateData
    //IPL_DEPTH_8U - 8位无符号整数
    //IPL_DEPTH_8S - 8位符号整数
    //IPL_DEPTH_16U - 16位无符号整数
    //IPL_DEPTH_16S - 16位符号整数
    //IPL_DEPTH_32S - 32位符号整数
    //IPL_DEPTH_32F - 单精度浮点数
    //IPL_DEPTH_64F - 双精度浮点数
    //cvSetZero(img);
    //cvReleaseImage( &img );
}
static void deallocate_image_memory(image_pointer *img_g)
{
    free(img_g->src_float);
    free(img_g->src_scale);
    free(img_g->src_temp);
    free(img_g->src_gray);
    free(img_g->src_smooth);
}
//
static void retinex_scales_distribution( gfloat* scales, gint nscales, gint mode, gint s )
{
    if (nscales == 1){
    /* For one filter we choose the median scale */
        scales[0] = (gfloat)( s / 2 );
    }
    else if (nscales == 2){
     /* For two filters whe choose the median and maximum scale */
      scales[0] = (gfloat) ( s / 2 );
      scales[1] = (gfloat) s;
    }
    else{
        gfloat size_step = (gfloat) s / (gfloat) nscales;
        gint i;
        switch(mode){

            case RETINEX_UNIFORM:
                for(i = 0; i < nscales; ++i)
                scales[i] = 2.F + (gfloat) i * size_step;
            break;

            case RETINEX_LOW:
                size_step = (gfloat) log(s - 2.0) / (gfloat) nscales;
                for (i = 0; i < nscales; ++i)
                scales[i] = 2.F + (gfloat) pow (10, (i * size_step) / log (10));
            break;

            case RETINEX_HIGH:
                size_step = (gfloat) log(s - 2.0) / (gfloat) nscales;
                for (i = 0; i < nscales; ++i)
                scales[i] = s - (gfloat) pow (10, (i * size_step) / log (10));
            break;

            default:
            break;
        }
    }
}
static void float_channel( guchar *src, gfloat *src_float, gfloat *src_gray, gint channelsize , gint bytes )
{
    gint i,j, pos,pos2;

    for ( i = 0, pos = 0, pos2 = 0; i < channelsize ; i++, pos += bytes, pos2 +=3 ){
        for ( j = 0 ; j < 3 ; j++){
            src_float[pos2 +j] = (gfloat)src[pos+j];
        }
        src_gray[i] = 0.3333F *(src_float[pos2]+src_float[pos2+1]+src_float[pos2+2]);
    }
}
static void log_channel( gfloat * src_gray, gfloat * src_scale,  gint channelsize )
{
    gint i;
    for (i = 0 ; i < channelsize; i++ ){
        src_scale[i] = (gfloat) log ( src_gray[i] + 1. );
    }
}

static void log_channel_lookup( gfloat * src_gray, gfloat * src_scale, gfloat* lookup_table_log, gint channelsize )
{
    gint i;
    gfloat c;
    for ( i = 0 ; i < channelsize; i++ ){
        c = src_gray[i] + 1.;
        src_scale[i] = lookup_table_log [ (guchar) (((c) > (256.)) ? (256.) : (((c) < (0.)) ? (0.) : (c))) ];
    }
}
static void gaussian_smoothing_2D( gfloat *src, gfloat *dst, gfloat sigma, gint width, gint height, gfloat *temp )
{
    gauss3_coefs coef;
    gint row,col, pos;
    compute_coefs3 ( &coef, sigma );
    /* Filtering (smoothing) Gaussian recursive.
     * Filter rows first
     */
    for (row=0 ;row < height; row++){
        pos =  row * width;
        gaussian_smoothing_1D ( src + pos, temp + pos, width, 1, &coef);
    }
    /* Filtering (smoothing) Gaussian recursive.
     * Second columns
     */
    for (col=0; col < width; col++){
        pos = col;
        gaussian_smoothing_1D( temp + pos, dst + pos, height, width, &coef);
    }
}
static void gaussian_smoothing_1D( gfloat *src, gfloat *dst, gint size, gint rowstride, gauss3_coefs *c )
{
    /*
     * Papers:  "Recursive Implementation of the gaussian filter.",
     * Ian T. Young , Lucas J. Van Vliet, Signal Processing 44, Elsevier 1995.
     * formula:
     *  9a forward filter
     *  9b backward filter
     *  fig7 algorithm
     */
    gint i,n, bufsize;
    gfloat *w1,*w2;
    /* forward pass */
    bufsize = size+3;
    size -= 1;
    w1 = (gfloat *) malloc (bufsize * sizeof (gfloat));
    w2 = (gfloat *) malloc (bufsize * sizeof (gfloat));
    w1[0] = src[0];
    w1[1] = src[0];
    w1[2] = src[0];
    for ( i = 0 , n = 3; i <= size ; i++, n++){
        w1[n] =
            ( gfloat )( c->B*src[i*rowstride] +
            (
                ( c->b[1] * w1[n-1] + c->b[2] * w1[n-2] + c->b[3] * w1[n-3] ) / c->b[0])
            );
    }

    /* backward pass */
    w2[size+1]= w1[size+3];
    w2[size+2]= w1[size+3];
    w2[size+3]= w1[size+3];
    for (i = size, n = i; i >= 0; i--, n--){
        w2[n] =
        dst[i * rowstride] =
            (gfloat)(c->B*w1[n+3] +
            (
                ( c->b[1]*w2[n+1] +c->b[2]*w2[n+2] +c->b[3]*w2[n+3] ) / c->b[0])
            );
    }
    free (w1);
    free (w2);
}
static void compute_coefs3( gauss3_coefs *c, gfloat sigma )
{
    /*
     * Papers:  "Recursive Implementation of the gaussian filter.",
     * Ian T. Young , Lucas J. Van Vliet, Signal Processing 44, Elsevier 1995.
     * formula:
     * 11b   computation of q
     * 8c    computation of b0..b1
     * 10    alpha is normalization constant B
     */
    gfloat q, q2, q3;

    if (sigma >= 2.5){
      q = 0.98711F * sigma - 0.96330F;
    }
    else if ((sigma >= 0.5) && (sigma < 2.5)){
      q = 3.97156F - 4.14554F * (gfloat) sqrt ((double) 1 - 0.26891 * sigma);
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
static void accumulation_minus_log(  gfloat * src_scale, gfloat * src_smooth , gfloat weight, gint channelsize )
{
    gint i;
    for (i = 0 ; i < channelsize; i++ ){
        src_scale[i] -= weight * (gfloat)log ( src_smooth[i]+ 1. );
    }
}
static void accumulation_minus_log_lookup(  gfloat * src_scale, gfloat * src_smooth , gfloat weight,gfloat* lookup_table_log, gint channelsize )
{
    gint i;
    gfloat c;

    for (i = 0 ; i < channelsize; i++ ){
        c = src_smooth[i] + 1.;


        src_scale[i] -= weight *  lookup_table_log [ (guchar) (((c) > (256.)) ? (256.) : (((c) < (0.)) ? (0.) : (c))) ];


    }


}
static void simplest_color_balance( gfloat *src_gray, gfloat mean, gfloat var, gfloat cvar, gint channelsize )
{
    gfloat mini,maxi,range;
    gint i;
    mini = mean - cvar * var;
    maxi = mean + cvar * var;
    range = maxi - mini;
    if (!range) range = 1.0;
    for (i = 0; i < channelsize; i++ ){
        src_gray[i] = 255 * ( src_gray[i] - mini ) / range;
    }
}
static void compute_mean_var(gfloat *src, gfloat *mean, gfloat *var, gint channelsize)
{
    gfloat vsquared;
    gint i;
    vsquared = 0;
    *mean = 0;
    for (i = 0; i < channelsize; i++){
        *mean += src[i];
        vsquared += src[i] * src[i];
    }
    *mean /= (gfloat) channelsize; /* mean */
    vsquared /= (gfloat) channelsize; /* mean (x^2) */
    *var = ( vsquared - (*mean * *mean) );
    *var = (gfloat)  sqrt(*var); /* var */
}
static void compute_color_from_grayscale( gfloat *src, gfloat *src_gray, gfloat *src1_gray, gint channelsize )
{
    gfloat *psrc;
    gfloat *psrc_gray;
    gfloat *psrc1_gray;
    gfloat  factor, max, c[3];
    gint pos;
    gint i, j;
    for ( i = 0, pos = 0; i < channelsize; i++, pos += 3){
        psrc = src + pos;
        psrc_gray = src_gray + i;
        psrc1_gray = src1_gray + i;
        if( psrc_gray[0] <= 1.) psrc_gray[0] = 1.;
        factor = psrc1_gray[0] / psrc_gray[0];
        if( factor > 3.) factor = 3.;
        for (j = 0 ; j < 3 ; j++) c[j] = factor * psrc[j];
        if( c[0] > 255. || c[1] > 255. || c[2] > 255.){
            max = psrc[0];
            if( psrc[1] > max) max=psrc[1];
            if( psrc[2] > max) max=psrc[2];
            factor = 255.F /max;
            for (j = 0 ; j < 3 ; j++) psrc[j] = factor * psrc[j];
        }
        else{
            for (j = 0 ; j < 3 ; j++) psrc[j] = c[j];
        }
    }
}
static void uchar_color( gfloat * src_float, guchar * src,  gint channelsize , gint bytes )
{
    gint i, j, pos,pos2;
    gfloat c;
    for ( i = 0, pos = 0, pos2=0; i < channelsize; i++, pos += bytes, pos2+=3){
        for (j = 0 ; j < 3 ; j++){
            c = src_float[j+pos2];
            src[j + pos] = (guchar) (((c) > (255.)) ? (255.) : (((c) < (0.)) ? (0.) : (c))) ;
        }
    }
}
