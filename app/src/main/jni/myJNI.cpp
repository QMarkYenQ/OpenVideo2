//
// Created by Mark on 2016/5/6.
//
//#include <opencv2/opencv.hpp>
#include <android/log.h>
//#include "com_example_mark_openvideo_myNDK.h"
//#include "EDLineDetector.h"
#include"HazeRemove.h"
#include "dehaze.h"
#include"contrast-retinex.h"


//#define min(a,b) (((a)<(b))?(a):(b))


#define LOG "Demo_JNI"
#define LOGD(...)  __android_log_print(ANDROID_LOG_DEBUG,LOG,__VA_ARGS__) // 定义LOGD类型
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <opencv2/opencv.hpp>
//#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "linefinder.h"
#include "edgedetector.h"
#include "com_example_mark_openvideo_myNDK.h"

using namespace cv;
//
// Create LineFinder instance
LineFinder ld;
//
Mat trans_refine;		//refine of trans
//
//
JNIEXPORT jint JNICALL Java_com_example_mark_openvideo_myNDK_iniProcess
    (JNIEnv *, jobject, jlong){
    // Set probabilistic Hough parameters  ---- we can fiddle with this later in order to improve accuracy (2nd semester) ...Alex
    ld.setLineLengthAndGap(40,5);//40,5
    ld.setMinVote(10); //10
    return 0;
    }


JNIEXPORT jint JNICALL Java_com_example_mark_openvideo_myNDK_imgProcess
    (JNIEnv *, jobject, jlong, jlong img){
    //-- Get a frame from the capture
    Mat* myrgba = (Mat*) img; // WT
    //-- Frame Captured from capture object
    Mat frame;
	cv::resize( *myrgba, frame, Size( 640,480 ), 0, 0, INTER_NEAREST);  // WT
    LOGD( "frame, rows =  %d,  cols = %d, channels = %d, depth = %d",
        frame.rows,frame.cols,frame.channels(),frame.depth());
    //-- If fram is not empty
    if( !frame.empty() ){
        //
		cv::Mat image;
		cv::resize(frame, image, Size(), 0.5, 0.5, INTER_NEAREST);
		//
		cv::Mat currentFrame;
		cv::resize(frame, currentFrame, Size(), 0.5, 0.5, INTER_NEAREST);
		//
		int height    = image.size().height;
  		int width     = image.size().width;
  		cv::Rect myROI(width*0.167, height/2, width*.833, height/2);
		image = image(myROI);
		//
		if (!image.data) return 0;
		// Compute Sobel
		EdgeDetector ed;
		ed.computeSobel(image);
		// Apply Canny algorithm
		cv::Mat contours;
		cv::Canny(image,contours, 125, 350);
		//cv::Canny(image, contours, 50, 200);
		// Detect lines
		std::vector<cv::Vec4i> li = ld.findLines(contours);
		// eliminate inconsistent lines   ----maybe use the slopes found within the while loops to eliminate certain lines?? ...Alex
		li = ld.removeLinesOfInconsistentOrientations( ed.getOrientation(), 0.4,0.15, image );
		//
		ld.drawDetectedLines(currentFrame);
        LOGD( "resizeFrame, rows =  %d,  cols = %d, channels = %d, depth = %d",
            currentFrame.rows,currentFrame.cols,currentFrame.channels(),currentFrame.depth());
        //
        LOGD( "inFrame, rows =  %d,  cols = %d, channels = %d, depth = %d",
            myrgba->rows, myrgba->cols, myrgba->channels(), myrgba->depth() );
        //
        //
		cv::resize( currentFrame, *myrgba, myrgba->size(), 0, 0, INTER_NEAREST);
		//
        LOGD( "outFrame, rows =  %d,  cols = %d, channels = %d, depth = %d",
            myrgba->rows, myrgba->cols, myrgba->channels(), myrgba->depth() );
	}else{
        LOGD(" --(!) No captured frame -- Break!"); //break;
	}
	return 0;
}
JNIEXPORT jint JNICALL Java_com_example_mark_openvideo_myNDK_HazeRemove
    (JNIEnv *, jobject, jlong, jlong img){
	const int radius = 3;
    //-- Get a frame from the capture
    Mat* myrgba = (Mat*) img; // WT
    //-- Frame Captured from capture object
    Mat frame;

    cvtColor( *myrgba, frame, CV_RGBA2RGB );
	cv::resize( frame, frame, Size( 320,240 ), 0, 0, INTER_NEAREST);  // WT
        LOGD( "resizeFrame, rows =  %d,  cols = %d, channels = %d, depth = %d",
            frame.rows,frame.cols,frame.channels(),frame.depth());
    //-- If fram is not empty
    if( !frame.empty() ){
    		HazeRemove hazeRemove(frame, radius, 0.95);
        	hazeRemove.getDarkChannelPrior();
		hazeRemove.getTransmissionMap();
		hazeRemove.getEstimatedTransmissionMap();
		Mat defog = hazeRemove.getHazeRemoveImg(hazeRemove.estimatedTransmissionMap);

		LOGD( "defog, rows =  %d,  cols = %d, channels = %d, depth = %d",
                    defog.rows,defog.cols,defog.channels(),defog.depth());
	//	imshow("defog", hazeRemove.darkChannelImg);
    	//	imshow("defog", hazeRemove.darkChannelImg);

    defog.convertTo( defog ,CV_8U,255.0,0.0);

    cvtColor( defog, defog, CV_RGB2RGBA );
		LOGD( "defog, rows =  %d,  cols = %d, channels = %d, depth = %d",
                    defog.rows,defog.cols,defog.channels(),defog.depth());
    	cv::resize( defog, *myrgba, myrgba->size(), 0, 0, INTER_NEAREST);

    }


  	return 0;
  }


JNIEXPORT jint JNICALL Java_com_example_mark_openvideo_myNDK_DehazeSimple
    (JNIEnv *, jobject, jlong, jlong myImg){
    //-- Get a frame from the capture
     Mat* myrgba = (Mat*) myImg; // WT
    //-- Frame Captured from capture object
     Mat frame;
    cvtColor( *myrgba, frame, CV_RGBA2RGB );
	cv::resize( frame, frame, Size( 320 ,240 ), 0, 0, INTER_LINEAR );  // WT
       LOGD( "resizeFrame, rows =  %d,  cols = %d, channels = %d, depth = %d",
            frame.rows,frame.cols,frame.channels(),frame.depth());

	Mat dark_channel;
	Mat trans;
	Mat img;
	Mat free_img;
	char filename[100];

	//clock_t start , finish ;
	//double duration1,duration2,duration3,duration4,duration5,duration6,duration7;
img = ReadImage( frame );

	//img = frame;

	//Calculate DarkChannelPrior
    dark_channel=DarkChannelPrior(img);
        LOGD( "dark_channel, rows =  %d,  cols = %d, channels = %d, depth = %d",
            dark_channel.rows,dark_channel.cols,dark_channel.channels(),dark_channel.depth());
	//Calculate Airlight
	Vec<float,3> a=Airlight(img,dark_channel);
	//Reading Refine Trans
	//Calculate Transmission Matrix
    trans = TransmissionMat( dark_channel );
	//trans_refine=ReadTransImage();
	//Haze Free
    free_img = hazefree(img,trans,a,0.2);
    LOGD( "hazefree, rows =  %d,  cols = %d, channels = %d, depth = %d",
        free_img.rows,free_img.cols,free_img.channels(),free_img.depth());
	free_img.convertTo(free_img,CV_8U,255.0,0.0);
    cvtColor( free_img, free_img, CV_RGB2RGBA );
    cv::resize( free_img, *myrgba, myrgba->size(), 0, 0, INTER_LINEAR);
    //
	LOGD( "free_img, rows =  %d,  cols = %d, channels = %d, depth = %d",
        free_img.rows,free_img.cols,free_img.channels(),free_img.depth());

  return 0;
}



JNIEXPORT jint JNICALL Java_com_example_mark_openvideo_myNDK_MSR_1original
    (JNIEnv *, jobject, jlong, jlong myImg){
    //-- Get a frame from the capture
    Mat* myrgba = (Mat*) myImg;
    //-- Frame Captured from capture object
    Mat frame;
    if( myrgba->rows >640 || myrgba->cols> 480){
        cv::resize( *myrgba, frame, Size( 1280 ,960 ), 0, 0, INTER_NEAREST );  // WT
        cv::pyrDown( frame, frame,  Size( 640, 480 ));
    }else{
        cv::resize( *myrgba, frame, Size( 640 ,480 ), 0, 0, INTER_LINEAR );  // WT
    }
    LOGD( "resizeFrame, rows =  %d,  cols = %d, channels = %d, depth = %d",
            frame.rows,frame.cols,frame.channels(),frame.depth());
    guchar *src = frame.ptr<uchar>(0);
    gint width  = frame.cols;
    gint height = frame.rows;
    gint bytes = 4;
    MSRCP( src, width, height, bytes, 2.3 );
    cv::resize( frame, *myrgba, myrgba->size(), 0, 0, INTER_LINEAR);
	LOGD( "free_img, rows =  %d,  cols = %d, channels = %d, depth = %d",
         myrgba->rows, myrgba->cols, myrgba->channels(), myrgba->depth());
    return 0;
        //cvtColor( *myrgba, frame, CV_RGBA2RGB );
    	//cv::resize( *myrgba, frame, Size( 320 ,240 ), 0, 0, INTER_LINEAR );  // WT
        //INTER_NEAREST , INTER_LINEAR and INTER_CUBIC
        //int fheight    = frame.size().height;
        //int fwidth     = frame.size().width;
        //cv::Rect myROI(fwidth*0, fheight*0, fwidth/2, fheight);
        //Mat frame2 = frame(myROI);
        //frame = frame2;
    	//free_img.convertTo(free_img,CV_8U,255.0,0.0);
        //cvtColor( frame, frame, CV_RGB2RGBA );
  }