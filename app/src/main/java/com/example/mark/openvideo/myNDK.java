package com.example.mark.openvideo;
import org.opencv.core.Mat;

/**
 * Created by Mark on 2016/5/7.
 */
public class myNDK
{
    static{
        System.loadLibrary("myJNI");
    }
    public native int imgProcess( long thiz, long inputImage );
    public native int iniProcess( long thiz );

    public native int HazeRemove( long thiz, long inputImage );

    public native int DehazeSimple( long thiz, long inputImage );

    public native int MSR_original( long thiz, long inputImage );


}
