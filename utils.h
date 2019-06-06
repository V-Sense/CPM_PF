/*  Name:
 *      utils.h
 *
 *  Description:
 *      Miscellaneous utilities, mainly conversion between different image formats / data types
 */
#pragma once
#ifndef UTILS_H
#define UTILS_H
#include "CPM_Tip2017Mod/CPM.h"
#include "CPM_Tip2017Mod/OpticFlowIO.h"
#include "PFilter/PermeabilityFilter.h"
extern "C" {
#include "PFilter/variational/variational.h"
}

void Match2Flow(FImage& inMat, FImage& ou, FImage& ov, int w, int h);

void Match2Mat2f(FImage matches, Mat2f &flow);

void WriteMatches(const char *filename, FImage& inMat);

void FImage2image_t(FImage mu, image_t* flow_x);

void Mat2f2image_t_uv(Mat2f flow, image_t* flow_x, image_t* flow_y);

void Mat3f2FImage(Mat3f in, FImage &out);

void Mat3f2color_image_t(Mat3f in, color_image_t *out);

void image_t2Mat(image_t* in, Mat &out);


#endif //!UTILS_H