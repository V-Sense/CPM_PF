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
#include "flow.h"
//#include <iostream>

bool unknown_flow(float u, float v) {
    return (fabs(u) >  UNKNOWN_FLOW_THRESH) 
	|| (fabs(v) >  UNKNOWN_FLOW_THRESH)
	|| std::isnan(u) || std::isnan(v);
}

bool unknown_flow(float *f) {
    return unknown_flow(f[0],f[1]);
}

/* ------------------------ OPTICAL FLOW INPUT / OUTPUT ----------------------- */
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


void setcols(int* colorwheel, int r, int g, int b, int k)
{
	colorwheel[k*3+0] = r;
	colorwheel[k*3+1] = g;
	colorwheel[k*3+2] = b;
}

int createcolorwheel(int* colorwheel)
{
	// relative lengths of color transitions:
	// these are chosen based on perceptual similarity
	// (e.g. one can distinguish more shades between red and yellow 
	//  than between yellow and green)
	int RY = 15;
	int YG = 6;
	int GC = 4;
	int CB = 11;
	int BM = 13;
	int MR = 6;
	int ncols = RY + YG + GC + CB + BM + MR;
	//printf("ncols = %d\n", ncols);
	if (ncols > MAXWHEELCOLS){
		printf("Too Many Columns in ColorWheel!\n");
		//exit(1);
	}
	int i;
	int k = 0;
	for (i = 0; i < RY; i++) setcols(colorwheel, 255,	   255*i/RY,	 0,	       k++);
	for (i = 0; i < YG; i++) setcols(colorwheel, 255-255*i/YG, 255,		 0,	       k++);
	for (i = 0; i < GC; i++) setcols(colorwheel, 0,		   255,		 255*i/GC,     k++);
	for (i = 0; i < CB; i++) setcols(colorwheel, 0,		   255-255*i/CB, 255,	       k++);
	for (i = 0; i < BM; i++) setcols(colorwheel, 255*i/BM,	   0,		 255,	       k++);
	for (i = 0; i < MR; i++) setcols(colorwheel, 255,	   0,		 255-255*i/MR, k++);

	return ncols;
}

void computeColor(double fx, double fy, unsigned char *pix, int* colorwheel, int ncols)
{
	double rad = sqrt(fx * fx + fy * fy);
	double a = atan2(-fy, -fx) / M_PI;
	double fk = (a + 1.0) / 2.0 * (ncols-1);
	int k0 = (int)fk;
	int k1 = (k0 + 1) % ncols;
	double f = fk - k0;
	//f = 0; // uncomment to see original color wheel
	for (int b = 0; b < 3; b++) {
		double col0 = colorwheel[k0*3+b] / 255.0;
		double col1 = colorwheel[k1*3+b] / 255.0;
		double col = (1 - f) * col0 + f * col1;
		if (rad <= 1)
			col = 1 - rad * (1 - col); // increase saturation with radius
		else
			col *= .75; // out of range
		pix[2 - b] = (int)(255.0 * col);
	}
	pix[3] = 0xff; // alpha channel, only for alignment
}

// convert flow to color
double MotionToColor(cv::Mat2f& flow, unsigned char* fillPix, float range /*= -1*/)
{
    // image size
    int w = flow.cols;
    int h = flow.rows;

	// determine motion range:
	double maxrad;

	if (range > 0) {
		maxrad = range;
	}else{	// obtain the motion range according to the max flow
		double maxu = -999, maxv = -999;
		double minu = 999, minv = 999;
		maxrad = -1;
		for (int i = 0; i < h; i++){
			for (int j = 0; j < w; j++){
				double u = flow(i,j)[0];
                double v = flow(i,j)[1];
				if (unknown_flow(u, v))
					continue;
				maxu = std::max(maxu, u);
				maxv = std::max(maxv, v);
				minu = std::min(minu, u);
				minv = std::min(minv, v);
				double rad = sqrt(u * u + v * v);
				maxrad = std::max(maxrad, rad);
			}
		}
		if (maxrad == 0) // if flow == 0 everywhere
			maxrad = 1;
	}

	//printf("max motion: %.2f  motion range: u = [%.2f,%.2f];  v = [%.2f,%.2f]\n",
	//	maxrad, minu, maxu, minv, maxv);

	int colorwheel[MAXWHEELCOLS*3];
	int ncols = createcolorwheel(colorwheel);

	for(int i=0; i<h; i++){
		for(int j=0; j<w; j++){
			int idx = i*w+j;
			double u = flow(i,j)[0];
			double v = flow(i,j)[1];
			if (unknown_flow(u, v)){
				memset(fillPix+idx*4, 0, 4);
				fillPix[idx*4 + 3] = 0xff; // alpha channel, only for alignment
			}else{
				double dx = std::min(std::max(u / maxrad, -1.0), 1.0);
				double dy = std::min(std::max(v / maxrad, -1.0), 1.0);
				computeColor(dx, dy, (unsigned char*)(fillPix + idx * 4), colorwheel, ncols);
			}
		}
	}

	return maxrad;
}

// write a 2-band image into a color coded png flow file
void WriteFlowAsImage(cv::Mat2f& flow, const char* imgName, float range /*= -1*/)
{
    int w = flow.cols;
    int h = flow.rows;

	cv::Mat img(h, w, CV_8UC4);
	float maxFlow = MotionToColor(flow, img.data, range);

#if 0
	// get corner color
	int x = 10, y = 20;
	unsigned char color[3];
	unsigned char* pSrc = img.data + y*img.step + x * 4;
	color[0] = 255 - pSrc[0];
	color[1] = 255 - pSrc[1];
	color[2] = 255 - pSrc[2];
	char info[256];
	sprintf(info, "max: %.1f", maxFlow);
    //cv::putText(img, info, cvPoint(x, y), CV_FONT_HERSHEY_SIMPLEX, 0.5, cvScalar(color[0], color[1], color[2]));
    cv::putText(img, info, cv::Point(x, y), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(color[0], color[1], color[2]));
#endif

	cv::imwrite(imgName, img);
}

/* ------------------------ OPERATIONS ON OPTICAL FLOWS  ----------------------- */


// Computes the normalized confidence map C between forward flow F and backward flow B
// First, the disntance D at every position (X,Y) is computed :
//     D(X,Y) = ||F(X,Y) - B(X + Fx(X,Y),y + Fy(X,Y))||
// Then, normalized confidence map C is computed :
//     C = 1 - D / max(D)
cv::Mat1f getFlowConfidence(cv::Mat2f forward_flow, cv::Mat2f backward_flow)
{
    cv::Mat1f distances = cv::Mat1f(forward_flow.rows, forward_flow.cols, -1);

    int h = forward_flow.rows;
    int w = forward_flow.cols;
    float max_distance = -1;

    // Computes the distance between forward and backwards flow
    #pragma omp parallel for
    for(int y = 0; y < h ; ++y)
    {
        for(int x = 0; x < w ; ++x)
        {
            // If there is forward flow for the position F(x,y)
            cv::Vec2f forward = forward_flow.at<cv::Vec2f>(y, x);
            if(forward[0] != kFLOW_UNKNOWN[0] && forward[1] != kFLOW_UNKNOWN[1])
            {
                cv::Vec2i next_position = getAbsoluteFlow(x, y, forward, h, w);
                if(next_position != kPOSITION_INVALID)
                {
                    // If there is backward flow for the refered position B(x + F(x,y).x,y + F(x,y).y)
                    cv::Vec2f backward = backward_flow.at<cv::Vec2f>(next_position[0], next_position[1]);
                    if(backward[0] != kFLOW_UNKNOWN[0] && backward[1] != kFLOW_UNKNOWN[1])
                    {
                        // computes the distance
                        float distance = (float)norm(forward + backward);

                        // Updates the max distance, if required
                        if(distance > max_distance)
                        {
                            #pragma omp critical
                            {
                                if(distance > max_distance)
                                    max_distance = distance;
                            }
                        }

                        // Updates the distance map
                        distances.at<float>(y, x) = distance;
                    }
                }
            }
        }
    }

    // If there is a difference between F and B
    if(max_distance > 0)
    {
        // Computes the normalized confidence map
        #pragma omp parallel for
        for(int y = 0; y < h ; ++y)
        {
            for(int x = 0; x < w ; ++x)
            {
                
                if(distances.at<float>(y, x) < 0)
                {
                    // Unknown flow, C = 0
                    distances.at<float>(y, x) = 0;
                    
                }
                else
                {
                    // C = 1 - normalized distance
                    //printf("max distance is %f\n", max_distance);
                    distances.at<float>(y, x) = (max_distance - distances.at<float>(y, x)) / max_distance;
                    //printf("result distance is %f!\n\n", distances.at<float>(y, x));
                }
                if (distances.at<float>(y, x) > 1 || distances.at<float>(y, x) < 0) printf("Error in confidence flow computation!\n");
            }
        }
        //printf("max distance is %d\n", max_distance);
        return distances;
    }
    else
    {
        // Forward and backwards flow are the same, C = 1
        return cv::Mat1f::ones(forward_flow.rows, forward_flow.cols);
    }
}
