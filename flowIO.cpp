// flow_io.cpp (modified from Middlebury tool)
//
// read and write our simple .flo flow file format

// ".flo" file format used for optical flow evaluation
//
// Stores 2-band float image for horizontal (u) and vertical (v) flow components.
// Floats are stored in little-endian order.
// A flow value is considered "unknown" if either |u| or |v| is greater than 1e9.
//
//  bytes  contents
//
//  0-3     tag: "PIEH" in ASCII, which in little endian happens to be the float 202021.25
//          (just a sanity check that floats are represented correctly)
//  4-7     width as an integer
//  8-11    height as an integer
//  12-end  data (width*height*2*4 bytes total)
//          the float values for u and v, interleaved, in row ordeatr, i.e.,
//          u[row0,col0], v[row0,col0], u[row0,col1], v[row0,col1], ...
//


// first four bytes, should be the same in little endian
#define TAG_FLOAT 202021.25  // check for this when READING the file
#define TAG_STRING "PIEH"    // use this when WRITING the file

//#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "imageLib.h"
#include "flowIO.h"
//#include <iostream>

bool unknown_flow(float u, float v) {
    return (fabs(u) >  UNKNOWN_FLOW_THRESH) 
	|| (fabs(v) >  UNKNOWN_FLOW_THRESH)
	|| std::isnan(u) || std::isnan(v);
}

bool unknown_flow(float *f) {
    return unknown_flow(f[0],f[1]);
}

// read a flow file into 2-band image
int ReadFlowFile(cv::Mat2f& img, const char* filename)
{
    //printf("bp1.1\n");
    if (filename == NULL) {
        printf("ReadFlowFile: empty filename\n");
        return 0;
    }

    const char *dot = strrchr(filename, '.');
    if (strcmp(dot, ".flo") != 0) {
        printf("ReadFlowFile (%s): extension .flo expected\n", filename);
        return 0;
    }

    FILE *stream = fopen(filename, "rb");
    if (stream == 0) {
        printf("ReadFlowFile: could not open %s\n", filename);
        return 0;
    }
    
    int width, height;
    float tag;

    if ((int)fread(&tag,    sizeof(float), 1, stream) != 1 ||
        (int)fread(&width,  sizeof(int),   1, stream) != 1 ||
        (int)fread(&height, sizeof(int),   1, stream) != 1) {    // read file head
        printf("ReadFlowFile: problem reading file %s\n", filename);
        return 0;
    }
    // printf("width in file is %d\n",width);
    // printf("height in file is %d\n",height);
    //printf("channels is %d\n",tag);


    if (tag != TAG_FLOAT) { // simple test for correct endian-ness
        printf("ReadFlowFile(%s): wrong tag (possibly due to big-endian machine?)\n", filename);
        return 0;
    }

    // another sanity check to see that integers were read correctly (99999 should do the trick...)
    if (width < 1 || width > 99999) {
        printf("ReadFlowFile(%s): illegal width %d\n", filename, width);
        return 0;
    }

    if (height < 1 || height > 99999) {
        printf("ReadFlowFile(%s): illegal height %d\n", filename, height);
        return 0;
    }

    int nBands = 2;

    cv::Mat2f img2(height,width,nBands);

    //printf("reading %d x %d x 2 = %d floats\n", width, height, width*height*2);
    
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            float tmp[nBands];
            if ((int)fread(tmp, sizeof(float), nBands, stream) != nBands){
                printf("ReadFlowFile(%s): file is too short\n", filename);
                return 0;
            }
	    
            img2.at<cv::Vec2f>(j,i)[0] = tmp[0];
            img2.at<cv::Vec2f>(j,i)[1] = tmp[1];
        }
    }

    if (fgetc(stream) != EOF) {
        printf("ReadFlowFile(%s): file is too long\n", filename);
        return 0;
    }
    
    fclose(stream);
    img = img2.clone();
    return 1;

}

// write a 2-band image into flow file
int WriteFlowFile(cv::Mat2f& img, const char* filename)
{
    if (filename == NULL) {
        printf("WriteFlowFile: empty filename\n");
        return 0;
    }

    const char *dot = strrchr(filename, '.');
    if (dot == NULL) {
        printf("WriteFlowFile: extension required in filename '%s'\n", filename);
        return 0;
    }

    if (strcmp(dot, ".flo") != 0) {
        printf("WriteFlowFile: filename '%s' should have extension '.flo'\n", filename);
        return 0;
    }

    //CShape sh = img.Shape();
    //int width = sh.width, height = sh.height, nBands = sh.nBands;
    cv::Size sze = img.size();
    int width = sze.width, height = sze.height, nBands = img.channels();
    //cout<<"width is"<<width<<endl;
    // printf("width is %d\n",width);
    // printf("height is %d\n",height);
    // printf("channels is %d\n",nBands);

    if (nBands != 2) {
        printf("WriteFlowFile(%s): image must have 2 bands\n", filename);
        return 0;
    }

    FILE *stream = fopen(filename, "wb");
    if (stream == 0) {
        printf("WriteFlowFile: could not open %s\n", filename);
        return 0;
    }

        // write the header
    fprintf(stream, TAG_STRING);
    if ((int)fwrite(&width,  sizeof(int),   1, stream) != 1 ||
        (int)fwrite(&height, sizeof(int),   1, stream) != 1) {
        printf("WriteFlowFile(%s): problem writing header\n", filename);
        return 0;
    }
    
    
        // write the rows
        for (int j=0; j<height; j++) {
            for (int i=0; i<width; i++) {
                float tmp[nBands];
                tmp[0] = img.at<cv::Vec2f>(j,i)[0];
                tmp[1] = img.at<cv::Vec2f>(j,i)[1];
                if ((int)fwrite(tmp, sizeof(float), nBands, stream) != nBands){
                    printf("WriteFlowFile(%s): problem writing data\n", filename);
                    return 0;
                }
            }
        }
    
        fclose(stream);
	return 1;
}


