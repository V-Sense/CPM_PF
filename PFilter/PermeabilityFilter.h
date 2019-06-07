#pragma once

#define _USE_MATH_DEFINES
#include <opencv2/opencv.hpp>
#include <cmath>
#include <assert.h>

#include "globals.h"
#include "flowIO.h"
#include "ImageIOpfm.h"

using namespace cv;
using namespace std;


const Vec2i kPOSITION_INVALID = Vec2i(-1, -1);
const float kMOVEMENT_UNKNOWN = 1e10;
const Vec2f kFLOW_UNKNOWN = Vec2f(kMOVEMENT_UNKNOWN,kMOVEMENT_UNKNOWN);

// Transforms a relative flow R into an absolute flow A, checking for margins
// Ax = X + Rx , Ay = Y + Ry
// if either
//   - Ax is outside [0,w] or
//   - Ay outside [0,y]
// it returns kPOSITION_INVALID
inline Vec2i getAbsoluteFlow(int x, int y, const Vec2f& flow, int h, int w)
{
    Vec2i result(cvRound(y + flow[1]), cvRound(x + flow[0]));
    if(result[0] >= 0 && result[0] < h && result[1] >= 0 && result[1] < w)
        return result;
    else
        return kPOSITION_INVALID;
}

// Computes the normalized confidence map C between forward flow F and backward flow B
// First, the disntance D at every position (X,Y) is computed :
//     D(X,Y) = ||F(X,Y) - B(X + Fx(X,Y),y + Fy(X,Y))||
// Then, normalized confidence map C is computed :
//     C = 1 - D / max(D)
inline Mat1f getFlowConfidence(Mat2f forward_flow, Mat2f backward_flow)
{
    Mat1f distances = Mat1f(forward_flow.rows, forward_flow.cols, -1);

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
            Vec2f forward = forward_flow.at<Vec2f>(y, x);
            if(forward[0] != kFLOW_UNKNOWN[0] && forward[1] != kFLOW_UNKNOWN[1])
            {
                Vec2i next_position = getAbsoluteFlow(x, y, forward, h, w);
                if(next_position != kPOSITION_INVALID)
                {
                    // If there is backward flow for the refered position B(x + F(x,y).x,y + F(x,y).y)
                    Vec2f backward = backward_flow.at<Vec2f>(next_position[0], next_position[1]);
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
        return Mat1f::ones(forward_flow.rows, forward_flow.cols);
    }
}



template <class TSrc>
Mat1f computeSpatialPermeability(Mat_<TSrc> src, float delta_XY, float alpha_XY)
{
    Mat_<TSrc> I = src;
    int h = I.rows;
    int w = I.cols;
    int num_channels = I.channels();

    // horizontal & vertical difference
    Mat_<TSrc> I_shifted = Mat_<TSrc>::zeros(h,w);
    Mat warp_mat = (Mat_<double>(2,3) <<1,0,-1, 0,1,0);
    warpAffine(I, I_shifted, warp_mat, I_shifted.size()); // could also use rect() to translate image

    // Equation 3.2 in Michel's Paper
    Mat_<TSrc> diff_perm = Mat_<TSrc>::zeros(h,w);
    diff_perm = I - I_shifted; // Ip - Ip' with 3 channels
    std::vector<Mat1f> diff_channels(num_channels);
    split(diff_perm, diff_channels);
    Mat1f dJdx = Mat1f::zeros(h, w);
    
    for (int c = 0; c < num_channels; c++) // ||Ip - Ip'|| via spliting 3 channels and do pow() and then sum them up and then sqrt()
    {
        Mat1f temp = Mat1f::zeros(h,w);
        temp = diff_channels[c];
        pow(temp, 2, temp);
        dJdx = dJdx + temp;
    }
    
    sqrt(dJdx, dJdx);
    
    Mat1f result = Mat1f::zeros(h,w); // finish the rest of the equation and return permeability map stored in "result"
    result = abs(dJdx / (sqrt(3) * delta_XY) );
    pow(result, alpha_XY, result);
    pow(1 + result, -1, result);
    
    return result;
}


// For single-channel target image J
template <class TSrc>
Mat1f filterXY(const Mat_<TSrc> I, const Mat1f J, const cpm_pf_params_t &cpm_pf_params)
{
    // Input image
    int h = I.rows;
    int w = I.cols;

    // initializations (move outside later)
    float iterations = cpm_pf_params.iterations_input_int;
    int lambda_XY = cpm_pf_params.lambda_XY_input_float;
    float delta_XY = cpm_pf_params.delta_XY_input_float;
    float alpha_XY = cpm_pf_params.alpha_XY_input_float;

    // spatial filtering
    Mat1f J_XY;
    J.copyTo(J_XY);
    Mat1f Mat_Ones = Mat1f::ones(1,1);

    //compute spatial permeability
    //compute horizontal filtered image
    Mat1f perm_horizontal;
    Mat1f perm_vertical;
    perm_horizontal = computeSpatialPermeability<TSrc>(I, delta_XY, alpha_XY);

    //compute vertical filtered image
    Mat_<TSrc> I_t = I.t();
    perm_vertical = computeSpatialPermeability<TSrc>(I_t, delta_XY, alpha_XY);
    perm_vertical = perm_vertical.t();

    for (int i = 0; i < iterations; ++i) {
        // spatial filtering
        // Equation 3.7~3.9 in Michel's thesis (lambda is 0 for flow map)

        // horizontal
        Mat1f J_XY_upper_h = Mat1f::zeros(h,w); //upper means upper of fractional number
        Mat1f J_XY_lower_h = Mat1f::zeros(h,w);
        for (int y = 0; y < h; y++) {
            Mat1f lp = Mat1f::zeros(1,w);
            Mat1f lp_normal = Mat1f::zeros(1,w);
            Mat1f rp = Mat1f::zeros(1,w);
            Mat1f rp_normal = Mat1f::zeros(1,w);

            // left pass
            for (int x = 1; x <= w-1; x++) {
                    lp(0, x) = perm_horizontal(y, x - 1) * (lp(0, x - 1) + J_XY(y, x - 1));
                    lp_normal(0, x) = perm_horizontal(y, x - 1) * (lp_normal(0, x - 1) + 1.0);
            }
        
            // right pass & combining
            for (int x = w-2; x >= 0; x--) {
                rp(0, x) = perm_horizontal(y, x) * (rp(0, x + 1) + J_XY(y, x + 1));
                rp_normal(0, x) = perm_horizontal(y, x) * (rp_normal(0, x + 1) + 1.0);
                
                //combination in right pass loop on-the-fly & deleted source image I
                if (x == w-2) {
                    J_XY(y, x+1) = (lp(0, x+1) + (1 - lambda_XY) * J_XY(y, x+1) + rp(0, x+1)) / (lp_normal(0, x+1) + 1.0 + rp_normal(0, x+1));
                }

                J_XY(y, x) = (lp(0, x) + (1 - lambda_XY) * J_XY(y, x) + rp(0, x)) / (lp_normal(0, x) + 1.0 + rp_normal(0, x));
            }
        }
            

        //vertical
        Mat1f J_XY_upper_v = Mat1f::zeros(h,w);
        Mat1f J_XY_lower_v = Mat1f::zeros(h,w);
        for (int x = 0; x < w; x++) {
            Mat1f dp = Mat1f::zeros(h,1);
            Mat1f dp_normal = Mat1f::zeros(h,1);
            Mat1f up = Mat1f::zeros(h,1);
            Mat1f up_normal = Mat1f::zeros(h,1);

            // (left pass) down pass
            for (int y = 1; y <= h-1; y++) {
                    dp(y, 0) = perm_vertical(y - 1, x) * (dp(y - 1, 0) + J_XY(y - 1, x));
                    dp_normal(y, 0) = perm_vertical(y - 1, x) * (dp_normal(y - 1, 0) + 1.0);
                }
    
            // (right pass) up pass & combining
            for (int y = h-2; y >= 0; y--) {
                up(y, 0) = perm_vertical(y, x) * (up(y + 1, 0) + J_XY(y + 1, x));
                up_normal(y, 0) = perm_vertical(y, x) * (up_normal(y + 1, 0) + 1.0);

                if(y == h-2) {
                    J_XY(y+1, x) = (dp(y+1, 0) + (1 - lambda_XY) * J_XY(y+1, x) + up(y+1, 0)) / (dp_normal(y+1, 0) + 1.0 + up_normal(y+1, 0));
                }
                J_XY(y, x) = (dp(y, 0) + (1 - lambda_XY) * J_XY(y, x) + up(y, 0)) / (dp_normal(y, 0) + 1.0 + up_normal(y, 0));
            }
        }
    }
    return J_XY;
}

// For multi-channel target image J
template <class TSrc, class TValue>
Mat_<TValue> filterXY(const Mat_<TSrc> I, const Mat_<TValue> J, const cpm_pf_params_t &cpm_pf_params)
{
    // Input image
    int h = I.rows;
    int w = I.cols;

    // initializations (move outside later)
    float iterations = cpm_pf_params.iterations_input_int;
    int lambda_XY = cpm_pf_params.lambda_XY_input_float;
    float delta_XY = cpm_pf_params.delta_XY_input_float;
    float alpha_XY = cpm_pf_params.alpha_XY_input_float;

    // spatial filtering
    int num_chs = J.channels();
    Mat_<TValue> J_XY;
    J.copyTo(J_XY);
    Mat_<TValue> Mat_Ones = Mat_<TValue>::ones(1,1);

    
    //compute spatial permeability
    //compute horizontal filtered image
    Mat1f perm_horizontal;
    Mat1f perm_vertical;
    perm_horizontal = computeSpatialPermeability<TSrc>(I, delta_XY, alpha_XY);

    //compute vertical filtered image
    Mat_<TSrc> I_t = I.t();
    perm_vertical = computeSpatialPermeability<TSrc>(I_t, delta_XY, alpha_XY);
    perm_vertical = perm_vertical.t();

    for (int i = 0; i < iterations; ++i) {
        // spatial filtering
        // Equation 3.7~3.9 in Michel's thesis (lambda is 0 for flow map)

        // horizontal
        Mat_<TValue> J_XY_upper_h = Mat_<TValue>::zeros(h,w); //upper means upper of fractional number
        Mat_<TValue> J_XY_lower_h = Mat_<TValue>::zeros(h,w);
        for (int y = 0; y < h; y++) {
            Mat_<TValue> lp = Mat_<TValue>::zeros(1,w);
            Mat_<TValue> lp_normal = Mat_<TValue>::zeros(1,w);
            Mat_<TValue> rp = Mat_<TValue>::zeros(1,w);
            Mat_<TValue> rp_normal = Mat_<TValue>::zeros(1,w);

            for (int c = 0; c < num_chs; c++) {
                // left pass
                for (int x = 1; x <= w-1; x++) {
                        lp(0, x)[c] = perm_horizontal(y, x - 1) * (lp(0, x - 1)[c] + J_XY(y, x - 1)[c]);
                        lp_normal(0, x)[c] = perm_horizontal(y, x - 1) * (lp_normal(0, x - 1)[c] + 1.0);
                }
            
                // right pass & combining
                for (int x = w-2; x >= 0; x--) {
                    rp(0, x)[c] = perm_horizontal(y, x) * (rp(0, x + 1)[c] + J_XY(y, x + 1)[c]);
                    rp_normal(0, x)[c] = perm_horizontal(y, x) * (rp_normal(0, x + 1)[c] + 1.0);
                    
                    //combination in right pass loop on-the-fly & deleted source image I
                    if (x == w-2) {
                        J_XY(y, x+1)[c] = (lp(0, x+1)[c] + (1 - lambda_XY) * J_XY(y, x+1)[c] + rp(0, x+1)[c]) / (lp_normal(0, x+1)[c] + 1.0 + rp_normal(0, x+1)[c]);
                    }

                    J_XY(y, x)[c] = (lp(0, x)[c] + (1 - lambda_XY) * J_XY(y, x)[c] + rp(0, x)[c]) / (lp_normal(0, x)[c] + 1.0 + rp_normal(0, x)[c]);
                }
            }
        }

        //vertical
        Mat_<TValue> J_XY_upper_v = Mat_<TValue>::zeros(h,w);
        Mat_<TValue> J_XY_lower_v = Mat_<TValue>::zeros(h,w);
        for (int x = 0; x < w; x++) {
            Mat_<TValue> dp = Mat_<TValue>::zeros(h,1);
            Mat_<TValue> dp_normal = Mat_<TValue>::zeros(h,1);
            Mat_<TValue> up = Mat_<TValue>::zeros(h,1);
            Mat_<TValue> up_normal = Mat_<TValue>::zeros(h,1);

            for (int c = 0; c < num_chs; c++) {
                // (left pass) down pass
                for (int y = 1; y <= h-1; y++) {
                        dp(y, 0)[c] = perm_vertical(y - 1, x) * (dp(y - 1, 0)[c] + J_XY(y - 1, x)[c]);
                        dp_normal(y, 0)[c] = perm_vertical(y - 1, x) * (dp_normal(y - 1, 0)[c] + 1.0);
                    }
        
                // (right pass) up pass & combining
                for (int y = h-2; y >= 0; y--) {
                    up(y, 0)[c] = perm_vertical(y, x) * (up(y + 1, 0)[c] + J_XY(y + 1, x)[c]);
                    up_normal(y, 0)[c] = perm_vertical(y, x) * (up_normal(y + 1, 0)[c] + 1.0);

                    if(y == h-2) {
                        J_XY(y+1, x)[c] = (dp(y+1, 0)[c] + (1 - lambda_XY) * J_XY(y+1, x)[c] + up(y+1, 0)[c]) / (dp_normal(y+1, 0)[c] + 1.0 + up_normal(y+1, 0)[c]);
                    }
                    J_XY(y, x)[c] = (dp(y, 0)[c] + (1 - lambda_XY) * J_XY(y, x)[c] + up(y, 0)[c]) / (dp_normal(y, 0)[c] + 1.0 + up_normal(y, 0)[c]);
                }
            }
        }
    }
    return J_XY;
}


template <class TSrc>
Mat1f computeTemporalPermeability(Mat_<TSrc> I, Mat_<TSrc> I_prev, Mat2f flow_XY, Mat2f flow_prev_XYT, float delta_photo, float delta_grad, float alpha_photo, float alpha_grad)
{
    int h = I.rows;
    int w = I.cols;
    int num_channels = I.channels();
    int num_channels_flow = flow_XY.channels();

    Mat1f perm_temporal = Mat1f::zeros(h,w);
    Mat1f perm_gradient = Mat1f::zeros(h,w);
    Mat1f perm_photo = Mat1f::zeros(h,w);

    Mat_<TSrc> I_prev_warped = Mat_<TSrc>::zeros(h,w);

    Mat2f prev_map = Mat2f::zeros(h,w);
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++){
            prev_map(y,x)[0] = x - flow_prev_XYT(y,x)[0];
            prev_map(y,x)[1] = y - flow_prev_XYT(y,x)[1];
        }
    }
    std::vector<Mat1f> prev_maps(num_channels_flow);
    split(prev_map, prev_maps);
    remap(I_prev, I_prev_warped, prev_maps[0], prev_maps[1], cv::INTER_CUBIC);

    // Equation 11
    Mat_<TSrc> diff_I = I - I_prev_warped;
    std::vector<Mat1f> diff_I_channels(num_channels);
    split(diff_I, diff_I_channels);
    Mat1f sum_diff_I = Mat1f::zeros(h, w);
    for (int c = 0; c < num_channels; c++)
    {
        Mat1f temp = Mat1f::zeros(h,w);
        temp = diff_I_channels[c];
        pow(temp, 2, temp);
        sum_diff_I = sum_diff_I + temp;
    }
    sqrt(sum_diff_I, sum_diff_I);
    pow(abs(sum_diff_I / (sqrt(3) * delta_photo)), alpha_photo, perm_photo);
    pow(1 + perm_photo, -1, perm_photo);

    Mat2f flow_prev_XYT_warped  = Mat2f::zeros(h,w);;
    remap(flow_prev_XYT, flow_prev_XYT_warped, prev_maps[0], prev_maps[1], cv::INTER_CUBIC);
    // Equation 12
    Mat2f diff_flow = flow_XY - flow_prev_XYT_warped;
    std::vector<Mat1f> diff_flow_channels(num_channels_flow);
    split(diff_flow, diff_flow_channels);
    Mat1f sum_diff_flow = Mat1f::zeros(h, w);
    for (int c = 0; c < num_channels_flow; c++)
    {
        Mat1f temp = Mat1f::zeros(h,w);
        temp = diff_flow_channels[c];
        pow(temp, 2, temp);
        sum_diff_flow = sum_diff_flow + temp;
    }
    sqrt(sum_diff_flow, sum_diff_flow);
    pow(abs(sum_diff_flow / (sqrt(2) * delta_grad)), alpha_grad, perm_gradient);
    pow(1 + perm_gradient, -1, perm_gradient);

    Mat perm_temporalMat = perm_photo.mul(perm_gradient);
    perm_temporalMat.convertTo(perm_temporal, CV_32FC1);

    return perm_temporal;
}

template <class TSrc, class TValue>
vector<Mat_<TValue> > filterT(Mat_<TSrc> src, Mat_<TSrc> src_prev, Mat_<TValue> J_XY, Mat_<TValue> J_prev_XY, Mat2f flow_XY, Mat2f flow_prev_XYT, Mat_<TValue> l_t_prev, Mat_<TValue> l_t_normal_prev)
{
    //store result variable
    vector<Mat_<TValue> > result;

    // Input image
    Mat_<TSrc> I = src;
    Mat_<TSrc> I_prev = src_prev;
    int h = I.rows;
    int w = I.cols;
    int num_channels = I.channels();
    int num_channels_flow = J_XY.channels();
/*
    // Joint image (optional).
    Mat_<TRef> A = I;
    if (!joint_image.empty())
    {
        // Input and joint images must have equal width and height.
        assert(src.size() == joint_image.size());
        A = joint_image;
    }
*/

    Mat_<TValue> J_XYT = J_XY;

    // Initialization parameters
    float lambda_T = 0;
    float delta_photo = 0.3;
    float delta_grad = 1.0;
    float alpha_photo = 2.0;
    float alpha_grad = 2.0;

    // temporal filtering
    int iterations = 1;


    //compute temporal permeability in this iteration
    Mat1f perm_temporal;
    perm_temporal = computeTemporalPermeability<TSrc>(I, I_prev, flow_XY, flow_prev_XYT, delta_photo, delta_grad, alpha_photo, alpha_grad);


    for (int i = 0; i < iterations; ++i)
    {
        //split flow
        //change flow file format to remap function required format
        Mat2f prev_map = Mat2f::zeros(h,w);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++){
                prev_map(y,x)[0] = x - flow_prev_XYT(y,x)[0];
                prev_map(y,x)[1] = y - flow_prev_XYT(y,x)[1];
            }
        }
        std::vector<Mat1f> flow_XYT_prev_maps(num_channels_flow);
        split(prev_map, flow_XYT_prev_maps);

        // temporal filtering
        // Equation 3.7~3.9 in Michel's thesis (lambda is 0 for flow map)

        // horizontal
        Mat_<TValue> l_t = Mat_<TValue>::zeros(h, w);
        Mat_<TValue> l_t_normal = Mat_<TValue>::zeros(h, w);
        
        // no need to do pixel-wise operation (via J(y,x)) since all operation is based on same location pixels
        // (left pass) forward pass & combining
        Mat_<TValue> temp_l_t_prev = Mat_<TValue>::zeros(h, w);
        Mat_<TValue> temp_l_t_prev_warped = Mat_<TValue>::zeros(h, w);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < num_channels_flow; c++) {
                    temp_l_t_prev(y,x)[c] = l_t_prev(y,x)[c] + J_prev_XY(y,x)[c];
                }
            }
        }
        //temp_l_t_prev = l_t_prev + J_prev_XY;
        remap(temp_l_t_prev, temp_l_t_prev_warped, flow_XYT_prev_maps[0], flow_XYT_prev_maps[1], cv::INTER_CUBIC);
        
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < num_channels_flow; c++) {
                    l_t(y,x)[c] = perm_temporal(y,x) * temp_l_t_prev_warped(y,x)[c];
                }
            }
        }
        
        Mat_<TValue> temp_l_t_normal_prev = Mat_<TValue>::zeros(h, w);
        Mat_<TValue> temp_l_t_normal_prev_warped = Mat_<TValue>::zeros(h, w);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < num_channels_flow; c++) {
                    temp_l_t_normal_prev(y,x)[c] = l_t_normal_prev(y,x)[c] + 1.0;
                }
            }
        }
        
        remap(temp_l_t_normal_prev, temp_l_t_normal_prev_warped, flow_XYT_prev_maps[0], flow_XYT_prev_maps[1], cv::INTER_CUBIC);
        
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < num_channels_flow; c++) {
                    l_t_normal(y,x)[c] = perm_temporal(y,x) * temp_l_t_normal_prev_warped(y,x)[c];
                    J_XYT(y,x)[c] = (l_t(y,x)[c] + (1 - lambda_T) * J_XY(y,x)[c]) / (l_t_normal(y,x)[c]+ 1.0);
                }
            }
        }

        result.push_back(l_t);
        result.push_back(l_t_normal);
    }

    result.push_back(J_XYT);
    return result;
}


