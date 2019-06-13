/**
 * ORIGINAL RELEASE November 2018
 * Author:   Yang Chen
 * Contact:  cheny5@scss.tcd.ie 
 * Institution:  V-SENSE, School of Computer Science, Trinity College Dublin
  * 
 * UPDATED June 2019
 * Author:   Martin Alain
 * Contact:  alainm@scss.tcd.ie 
 * Institution:  V-SENSE, School of Computer Science, Trinity College Dublin
 */

#pragma once
#ifndef PF_H
#define PF_H

#define _USE_MATH_DEFINES
#include <opencv2/opencv.hpp>
#include <cmath>
#include <assert.h>

using namespace cv;
using namespace std;

template <class TI>
class PermeabilityFilter
{
private:
    // Guide images
    Mat_<TI> I_XY; // spatial
    Mat_<TI> I_T0, I_T1; // temporal
    bool is_I_XY_set, is_I_T_set;

    // Permeability maps
    bool is_perm_xy_set, is_perm_t_set;
        
public:
    PermeabilityFilter(); // Initializes default parameters

    // Guide images
    void set_I_XY(const Mat_<TI> I); // spatial
    void set_I_T(const Mat_<TI> I0, const Mat_<TI> I1); // temporal

    // Permeability maps
    Mat1f perm_h, perm_v, perm_t;

    // Spatial parameters
    int iter_XY;
    int lambda_XY;
    float sigma_XY;
    float alpha_XY;

    Mat1f computeSpatialPermeability(const Mat_<TI> I);
    void computeSpatialPermeabilityMaps();
    
    Mat1f filterXY(const Mat1f J); // For single-channel target image J
    template <class TJ>
    Mat_<TJ> filterXY(const Mat_<TJ> J); // For multi-channel target image J


    // Temporal parameters
    int iter_T;
    float lambda_T;
    float sigma_photo;
    float sigma_grad;
    float alpha_photo;
    float alpha_grad;
    
    void computeTemporalPermeability(const Mat2f flow_XY, const Mat2f flow_prev_XYT);
    
    template <class TJ>
    vector<Mat_<TJ> > filterT(  const Mat_<TJ> J_XY, const Mat_<TJ> J_prev_XY,
                                const Mat2f flow_XY, const Mat2f flow_prev_XYT,
                                const Mat_<TJ> l_t_prev, const Mat_<TJ> l_t_normal_prev);
};
#endif //! PF_H



/* ---------------- Set parameters --------------------------- */
template <class TI>
PermeabilityFilter<TI>::PermeabilityFilter()
{
    // Init default parameters
    // Spatial
    iter_XY = 5;
    lambda_XY = 0.0f;
    sigma_XY = 0.017f;
    alpha_XY = 2.0f;

    // Temporal
    iter_T = 1;
    lambda_T = 0.0f;
    sigma_photo = 0.3f;
    sigma_grad = 1.0f;
    alpha_photo = 2.0f;
    alpha_grad = 2.0f;

    // Guide images
    is_I_XY_set = false;
    is_I_T_set = false;

    // Permeability maps
    is_perm_xy_set = false;
    is_perm_t_set = false;
}

// Guide images
// spatial
template <class TI>
void PermeabilityFilter<TI>::set_I_XY(const Mat_<TI> I)
{
    I_XY = I;
    is_I_XY_set = true;
}

// temporal
template <class TI>
void PermeabilityFilter<TI>::set_I_T(const Mat_<TI> I0, const Mat_<TI> I1)
{
    I_T0 = I0;
    I_T1 = I1;
    is_I_T_set = true;
}

/* ---------------- Spatial filtering --------------------------- */
template <class TI>
Mat1f PermeabilityFilter<TI>::computeSpatialPermeability(Mat_<TI> I)
{
    int h = I.rows;
    int w = I.cols;
    int num_channels = I.channels();

    // horizontal & vertical difference
    Mat_<TI> I_shifted = Mat_<TI>::zeros(h,w);
    Mat warp_mat = (Mat_<double>(2,3) <<1,0,-1, 0,1,0);
    warpAffine(I, I_shifted, warp_mat, I_shifted.size()); // could also use rect() to translate image

    // Equation (2) in paper "Towards Edge-Aware Spatio-Temporal Filtering in Real-Time"
    Mat_<TI> diff_perm = Mat_<TI>::zeros(h,w);
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
    result = abs(dJdx / (sqrt(3) * sigma_XY) );
    pow(result, alpha_XY, result);
    // pow(1 + result, -1, result);
    result = 1 / (1 + result);
    
    return result;
}

template <class TI>
void PermeabilityFilter<TI>::computeSpatialPermeabilityMaps()
{
    if(! is_I_XY_set)
    {
        cerr << "Can not compute spatial permeability maps, spatial guide image is not defined." << endl;
        exit (EXIT_FAILURE);
        // return;
    }

    //compute horizontal filtered image
    perm_h = computeSpatialPermeability(I_XY);

    //compute vertical filtered image
    perm_v = computeSpatialPermeability(I_XY.t());
    perm_v = perm_v.t();
    
    is_perm_xy_set = true;
}


// For single-channel target image J
template <class TI>
Mat1f PermeabilityFilter<TI>::filterXY(const Mat1f J)
{
    if(! is_perm_xy_set)
    {
        cerr << "Can not compute spatial permeability filtering, spatial permeability maps are not computed." << endl;
        exit (EXIT_FAILURE);
        // return Mat1f::zeros(J.rows, J.cols);
    }

    // Input image
    int h = J.rows;
    int w = J.cols;

    // spatial filtering
    Mat1f J_XY;
    J.copyTo(J_XY);
    
    for (int i = 0; i < iter_XY; ++i) {
        // spatial filtering
        // Equation (5) and (6) in paper "Towards Edge-Aware Spatio-Temporal Filtering in Real-Time"

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
                    lp(0, x) = perm_h(y, x - 1) * (lp(0, x - 1) + J_XY(y, x - 1));
                    lp_normal(0, x) = perm_h(y, x - 1) * (lp_normal(0, x - 1) + 1.0);
            }
        
            // right pass & combining
            for (int x = w-2; x >= 0; x--) {
                rp(0, x) = perm_h(y, x) * (rp(0, x + 1) + J_XY(y, x + 1));
                rp_normal(0, x) = perm_h(y, x) * (rp_normal(0, x + 1) + 1.0);
                
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
                    dp(y, 0) = perm_v(y - 1, x) * (dp(y - 1, 0) + J_XY(y - 1, x));
                    dp_normal(y, 0) = perm_v(y - 1, x) * (dp_normal(y - 1, 0) + 1.0);
                }
    
            // (right pass) up pass & combining
            for (int y = h-2; y >= 0; y--) {
                up(y, 0) = perm_v(y, x) * (up(y + 1, 0) + J_XY(y + 1, x));
                up_normal(y, 0) = perm_v(y, x) * (up_normal(y + 1, 0) + 1.0);

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
template <class TI>
template <class TJ>
Mat_<TJ> PermeabilityFilter<TI>::filterXY(const Mat_<TJ> J)
{
    if(! is_perm_xy_set)
    {
        cerr << "Can not compute spatial permeability filtering, spatial permeability maps are not computed." << endl;
        exit (EXIT_FAILURE);
        // return Mat1f::zeros(J.rows, J.cols);
    }

    // Input image
    int h = J.rows;
    int w = J.cols;

    // spatial filtering
    int num_chs = J.channels();
    Mat_<TJ> J_XY;
    J.copyTo(J_XY);
    
    
    for (int i = 0; i < iter_XY; ++i) {
        // spatial filtering
        // Equation (5) and (6) in paper "Towards Edge-Aware Spatio-Temporal Filtering in Real-Time"

        // horizontal
        Mat_<TJ> J_XY_upper_h = Mat_<TJ>::zeros(h,w); //upper means upper of fractional number
        Mat_<TJ> J_XY_lower_h = Mat_<TJ>::zeros(h,w);
        for (int y = 0; y < h; y++) {
            Mat_<TJ> lp = Mat_<TJ>::zeros(1,w);
            Mat_<TJ> lp_normal = Mat_<TJ>::zeros(1,w);
            Mat_<TJ> rp = Mat_<TJ>::zeros(1,w);
            Mat_<TJ> rp_normal = Mat_<TJ>::zeros(1,w);

            for (int c = 0; c < num_chs; c++) {
                // left pass
                for (int x = 1; x <= w-1; x++) {
                        lp(0, x)[c] = perm_h(y, x - 1) * (lp(0, x - 1)[c] + J_XY(y, x - 1)[c]);
                        lp_normal(0, x)[c] = perm_h(y, x - 1) * (lp_normal(0, x - 1)[c] + 1.0);
                }
            
                // right pass & combining
                for (int x = w-2; x >= 0; x--) {
                    rp(0, x)[c] = perm_h(y, x) * (rp(0, x + 1)[c] + J_XY(y, x + 1)[c]);
                    rp_normal(0, x)[c] = perm_h(y, x) * (rp_normal(0, x + 1)[c] + 1.0);
                    
                    //combination in right pass loop on-the-fly & deleted source image I
                    if (x == w-2) {
                        J_XY(y, x+1)[c] = (lp(0, x+1)[c] + (1 - lambda_XY) * J_XY(y, x+1)[c] + rp(0, x+1)[c]) / (lp_normal(0, x+1)[c] + 1.0 + rp_normal(0, x+1)[c]);
                    }

                    J_XY(y, x)[c] = (lp(0, x)[c] + (1 - lambda_XY) * J_XY(y, x)[c] + rp(0, x)[c]) / (lp_normal(0, x)[c] + 1.0 + rp_normal(0, x)[c]);
                }
            }
        }

        //vertical
        Mat_<TJ> J_XY_upper_v = Mat_<TJ>::zeros(h,w);
        Mat_<TJ> J_XY_lower_v = Mat_<TJ>::zeros(h,w);
        for (int x = 0; x < w; x++) {
            Mat_<TJ> dp = Mat_<TJ>::zeros(h,1);
            Mat_<TJ> dp_normal = Mat_<TJ>::zeros(h,1);
            Mat_<TJ> up = Mat_<TJ>::zeros(h,1);
            Mat_<TJ> up_normal = Mat_<TJ>::zeros(h,1);

            for (int c = 0; c < num_chs; c++) {
                // (left pass) down pass
                for (int y = 1; y <= h-1; y++) {
                        dp(y, 0)[c] = perm_v(y - 1, x) * (dp(y - 1, 0)[c] + J_XY(y - 1, x)[c]);
                        dp_normal(y, 0)[c] = perm_v(y - 1, x) * (dp_normal(y - 1, 0)[c] + 1.0);
                    }
        
                // (right pass) up pass & combining
                for (int y = h-2; y >= 0; y--) {
                    up(y, 0)[c] = perm_v(y, x) * (up(y + 1, 0)[c] + J_XY(y + 1, x)[c]);
                    up_normal(y, 0)[c] = perm_v(y, x) * (up_normal(y + 1, 0)[c] + 1.0);

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



/* ---------------- Temporal filtering --------------------------- */
template <class TI>
void PermeabilityFilter<TI>::computeTemporalPermeability(const Mat2f flow_XY, const Mat2f flow_prev_XYT)
{
    if(! is_I_T_set)
    {
        cerr << "Can not compute temporal permeability maps, spatial guide images are not defined." << endl;
        exit (EXIT_FAILURE);
        // return;
    }

    int h = I_T1.rows;
    int w = I_T1.cols;
    int num_channels = I_T1.channels();
    int num_channels_flow = flow_XY.channels();

    Mat1f perm_temporal = Mat1f::zeros(h,w);
    Mat1f perm_gradient = Mat1f::zeros(h,w);
    Mat1f perm_photo = Mat1f::zeros(h,w);

    Mat_<TI> I_T0_warped = Mat_<TI>::zeros(h,w);

    Mat prev_map_x = Mat::zeros(h,w, CV_32FC1);
    Mat prev_map_y = Mat::zeros(h,w, CV_32FC1);
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++){
            prev_map_x.at<float>(y,x) = x - flow_prev_XYT(y,x)[0];
            prev_map_y.at<float>(y,x) = y - flow_prev_XYT(y,x)[1];
        }
    }
    
    Mat map_x = Mat::zeros(h,w, CV_16SC2);
    Mat map_y = Mat::zeros(h,w, CV_16UC1);
    convertMaps(prev_map_x, prev_map_y, map_x, map_y, CV_16SC2);
    
    remap(I_T0, I_T0_warped, map_x, map_y, cv::INTER_CUBIC);

    // Mat2f prev_map = Mat2f::zeros(h,w);
    // for (int y = 0; y < h; y++) {
    //     for (int x = 0; x < w; x++){
    //         prev_map(y,x)[0] = x - flow_prev_XYT(y,x)[0];
    //         prev_map(y,x)[1] = y - flow_prev_XYT(y,x)[1];
    //     }
    // }
    // std::vector<Mat1f> prev_maps(num_channels_flow);
    // split(prev_map, prev_maps);
    // remap(I_prev, I_T0_warped, prev_maps[0], prev_maps[1], cv::INTER_CUBIC);

    // Mat I_org, I_warp;
    // I.convertTo(I_org, CV_8UC3, 255);
    // I_T0_warped.convertTo(I_warp, CV_8UC3, 255);
    
    // int rnum = rand();
    // std::ostringstream oss;
    // oss << rnum;

    // imwrite("original_image" + oss.str() + ".png", I_org);
    // imwrite("warp_image" + oss.str() + ".png", I_warp);

    // Equation (11) in paper "Towards Edge-Aware Spatio-Temporal Filtering in Real-Time"
    Mat_<TI> diff_I = I_T1 - I_T0_warped;
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
    pow(abs(sum_diff_I / (sqrt(3) * sigma_photo)), alpha_photo, perm_photo);
    // pow(1 + perm_photo, -1, perm_photo);
    perm_photo = 1 / (1 + perm_photo);

    Mat2f flow_prev_XYT_warped  = Mat2f::zeros(h,w);;
    // remap(flow_prev_XYT, flow_prev_XYT_warped, prev_maps[0], prev_maps[1], cv::INTER_CUBIC);
    remap(flow_prev_XYT, flow_prev_XYT_warped, map_x, map_y, cv::INTER_CUBIC);
    // Equation (12) in paper "Towards Edge-Aware Spatio-Temporal Filtering in Real-Time"
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
    pow(abs(sum_diff_flow / (sqrt(2) * sigma_grad)), alpha_grad, perm_gradient);
    // pow(1 + perm_gradient, -1, perm_gradient);
    perm_gradient = 1 / (1 + perm_gradient);

    Mat perm_temporalMat = perm_photo.mul(perm_gradient);

    perm_temporalMat.convertTo(perm_t, CV_32FC1);
    is_perm_t_set = true;
}


template <class TI>
template <class TJ>
vector<Mat_<TJ> > PermeabilityFilter<TI>::filterT(  const Mat_<TJ> J_XY, const Mat_<TJ> J_prev_XY,
                                                    const Mat2f flow_XY, const Mat2f flow_prev_XYT,
                                                    const Mat_<TJ> l_t_prev, const Mat_<TJ> l_t_normal_prev)
{
    if(! is_perm_t_set)
    {
        cerr << "Can not compute temporal permeability filtering, temporal permeability map is not computed." << endl;
        exit (EXIT_FAILURE);
        // return Mat1f::zeros(J.rows, J.cols);
    }

    //store result variable
    vector<Mat_<TJ> > result;

    // Input image
    int h = J_XY.rows;
    int w = J_XY.cols;
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

    Mat_<TJ> J_XYT;
    J_XY.copyTo(J_XYT);


    for (int i = 0; i < iter_T; ++i)
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
        Mat_<TJ> l_t = Mat_<TJ>::zeros(h, w);
        Mat_<TJ> l_t_normal = Mat_<TJ>::zeros(h, w);
        
        // no need to do pixel-wise operation (via J(y,x)) since all operation is based on same location pixels
        // forward pass & combining
        Mat_<TJ> temp_l_t_prev = Mat_<TJ>::zeros(h, w);
        Mat_<TJ> temp_l_t_prev_warped = Mat_<TJ>::zeros(h, w);
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
                    l_t(y,x)[c] = perm_t(y,x) * temp_l_t_prev_warped(y,x)[c];
                }
            }
        }
        
        Mat_<TJ> temp_l_t_normal_prev = Mat_<TJ>::zeros(h, w);
        Mat_<TJ> temp_l_t_normal_prev_warped = Mat_<TJ>::zeros(h, w);
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
                    l_t_normal(y,x)[c] = perm_t(y,x) * temp_l_t_normal_prev_warped(y,x)[c];
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


