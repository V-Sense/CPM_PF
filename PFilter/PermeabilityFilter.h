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

const float sqrt_2 = sqrt(2);
const float sqrt_3 = sqrt(3);

template <class TI>
class PermeabilityFilter
{
private:
    
    // Spatial parameters
    Mat_<TI> I_XY; // Guide image
    
    bool is_I_XY_set;
    bool is_perm_xy_set;

    // Temporal parameters
    Mat _l_t0_num, _l_t0_den; // Accumulated left pass buffer
    Mat2f flow_t0, flow_t1;
    Mat1f disp_t0, disp_t1;
    Mat map_t0_to_t1_x, map_t0_to_t1_y;
    Mat_<TI> I_t0, I_t1; // Guide images
    
    bool is_l_init;
    bool is_flow_set;
    bool is_disp_set;
    bool is_I_T_set;
    bool is_perm_t_set;
        
public:
    PermeabilityFilter(); // Initializes default parameters

    
    // Spatial parameters
    int iter_XY;
    int lambda_XY;
    float sigma_XY;
    float alpha_XY;

    Mat1f perm_v, perm_h; // Permeability maps
    void set_I_XY(const Mat_<TI> I); // Guide images

    void computeSpatialPermeability(string spatialDir);
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

    Mat1f perm_t; // Permeability map

    template <class TJ>
    void init_T(int h, int w);
    void set_I_T(const Mat_<TI> I0, const Mat_<TI> I1); // Guide images
    void set_flow_T(const Mat2f f0, const Mat2f f1);
    void set_flow_off();
    void set_disp_T(const Mat1f d0, const Mat1f d1, const string parallax);
    void set_disp_off();
    
    void computeTemporalPermeability();
    
    Mat1f filterT(const Mat1f J_t0_XYT, const Mat1f J_t1_XY, Mat1f &l_t0_num, Mat1f &l_t0_den); // For single-channel target images J
    template <class TJ>
    Mat_<TJ> filterT(const Mat_<TJ> J_t0_XYT, const Mat_<TJ> J_t1_XY, Mat_<TJ> &l_t0_num, Mat_<TJ> &l_t0_den); // For multi-channel target images J
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

    is_I_XY_set = false;
    is_perm_xy_set = false;

    // Temporal
    iter_T = 1;
    lambda_T = 0.0f;
    sigma_photo = 0.3f;
    sigma_grad = 1.0f;
    alpha_photo = 2.0f;
    alpha_grad = 2.0f;

    is_I_T_set = false;
    is_l_init = false;
    is_flow_set = false;
    is_disp_set = false;
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
    I_t0 = I0;
    I_t1 = I1;
    is_I_T_set = true;
}


// Temporal filter init
template <class TI>
template <class TJ>
void PermeabilityFilter<TI>::init_T(int h, int w)
{
    _l_t0_num = Mat_<TJ>::zeros(h, w);
    _l_t0_den = Mat_<TJ>::zeros(h, w);
    is_l_init = true;
}


template <class TI>
void PermeabilityFilter<TI>::set_flow_T(const Mat2f f0, const Mat2f f1)
{
    if(is_disp_set)
    {
        cerr << "Disparity map is already set in permeability filter, cannot set both flow and disparity." << endl;
        exit(EXIT_FAILURE);
    }

    flow_t0 = f0;
    flow_t1 = f1;
    int h = f0.rows;
    int w = f0.cols;

    Mat map_x_32FC1 = Mat::zeros(h,w, CV_32FC1);
    Mat map_y_32FC1 = Mat::zeros(h,w, CV_32FC1);
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++){
            map_x_32FC1.at<float>(y,x) = x - flow_t0(y,x)[0];
            map_y_32FC1.at<float>(y,x) = y - flow_t0(y,x)[1];
        }
    }
    
    map_t0_to_t1_x = Mat::zeros(h,w, CV_16SC2);
    map_t0_to_t1_y = Mat::zeros(h,w, CV_16UC1);
    convertMaps(map_x_32FC1, map_y_32FC1, map_t0_to_t1_x, map_t0_to_t1_y, CV_16SC2);

    is_flow_set = true;
}

template <class TI>
void PermeabilityFilter<TI>::set_flow_off(){
    is_flow_set = false;
}


template <class TI>
void PermeabilityFilter<TI>::set_disp_T(const Mat1f d0, const Mat1f d1, const string parallax)
{
    if(is_flow_set)
    {
        cerr << "Optical flow is already set in permeability filter, cannot set both flow and disparity." << endl;
        exit(EXIT_FAILURE);
    }

    if(parallax != "ver" || parallax != "vertical" || parallax != "hor" || parallax != "horizontal") 
	{
		cerr << "Wrong parallax direction when setting disparity maps in permeabilty filter, should be ver or vertical or hor or horizontal" << endl;
		exit(EXIT_FAILURE);
	}

    disp_t0 = d0;
    disp_t1 = d1;
    int h = d0.rows;
    int w = d0.cols;

    Mat map_x_32FC1 = Mat::zeros(h,w, CV_32FC1);
    Mat map_y_32FC1 = Mat::zeros(h,w, CV_32FC1);

    if(parallax == "hor" || parallax == "horizontal") 
    {
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++){
                map_x_32FC1.at<float>(y,x) = x - disp_t0(y,x);
                map_y_32FC1.at<float>(y,x) = y;
            }
        }
    }
    else if(parallax == "ver" || parallax == "horizontal")
    {
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++){
                map_x_32FC1.at<float>(y,x) = x;
                map_y_32FC1.at<float>(y,x) = y - disp_t0(y,x);
            }
        }
    }
    
    map_t0_to_t1_x = Mat::zeros(h,w, CV_16SC2);
    map_t0_to_t1_y = Mat::zeros(h,w, CV_16UC1);
    convertMaps(map_x_32FC1, map_y_32FC1, map_t0_to_t1_x, map_t0_to_t1_y, CV_16SC2);

    is_disp_set = true;
}

template <class TI>
void PermeabilityFilter<TI>::set_disp_off(){
    is_disp_set = false;
}


/* ---------------- Spatial filtering --------------------------- */
template <class TI>
void PermeabilityFilter<TI>::computeSpatialPermeability(string spatialDir)
{
    Mat_<TI> I;
    if (spatialDir == "ver" || spatialDir == "vertical")
        I = I_XY.t();
    else if (spatialDir == "hor" || spatialDir == "horizontal")
        I = I_XY;
    else
    {
        cerr << "Wrong spatial direction for spatial permeability computation, should be hor, or horizontal, or ver, or vertical." << endl;
        exit(EXIT_FAILURE);
    }
    

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
    
    if (spatialDir == "ver" || spatialDir == "vertical")
        perm_v = result.t();
    else if (spatialDir == "hor" || spatialDir == "horizontal")
        perm_h = result;
    
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

    computeSpatialPermeability("hor");
    computeSpatialPermeability("ver");

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
void PermeabilityFilter<TI>::computeTemporalPermeability()
{
    if(! is_I_T_set)
    {
        cerr << "Can not compute temporal permeability maps, spatial guide images are not defined." << endl;
        exit (EXIT_FAILURE);
    }
    if(! is_flow_set && ! is_disp_set)
    {
        cerr << "Can not compute temporal permeability maps, optical flows / disparity maps are not defined." << endl;
        exit (EXIT_FAILURE);
    }
    if( is_flow_set && is_disp_set)
    {
        cerr << "Can not compute temporal permeability maps, optical flows and disparity maps are both set." << endl;
        exit (EXIT_FAILURE);
    }

    int h = I_t1.rows;
    int w = I_t1.cols;
    int num_channels = I_t1.channels();
    int num_channels_flow = flow_t1.channels(); // 2

    Mat1f perm_temporal = Mat1f::zeros(h,w);
    Mat1f perm_gradient = Mat1f::zeros(h,w);
    Mat1f perm_photo = Mat1f::zeros(h,w);

    // Equation (11) in paper "Towards Edge-Aware Spatio-Temporal Filtering in Real-Time"
    Mat_<TI> I_t0_warp_2_t1 = Mat_<TI>::zeros(h,w);
    remap(I_t0, I_t0_warp_2_t1, map_t0_to_t1_x, map_t0_to_t1_y, cv::INTER_CUBIC);
    
    Mat_<TI> diff_I2 = I_t1 - I_t0_warp_2_t1;
    diff_I2 = diff_I2.mul(diff_I2); // Square

    std::vector<Mat1f> diff_I2_channels(num_channels);
    split(diff_I2, diff_I2_channels);
    Mat1f sum_diff_I2 = Mat1f::zeros(h, w);
    for (int c = 0; c < num_channels; c++) sum_diff_I2 = sum_diff_I2 + diff_I2_channels[c];
    
    sqrt(sum_diff_I2, sum_diff_I2);
    
    pow(abs(sum_diff_I2 / (sqrt_3 * sigma_photo)), alpha_photo, perm_photo);
    perm_photo = 1 / (1 + perm_photo);


    // Equation (12) in paper "Towards Edge-Aware Spatio-Temporal Filtering in Real-Time"
    Mat1f sum_diff_flow2;
    if( is_flow_set )
    {
        Mat2f flow_XY_t0_warp_2_t1  = Mat2f::zeros(h,w);;
        remap(flow_t0, flow_XY_t0_warp_2_t1, map_t0_to_t1_x, map_t0_to_t1_y, cv::INTER_CUBIC);
        
        Mat2f diff_flow2 = flow_t1 - flow_XY_t0_warp_2_t1;
        diff_flow2 = diff_flow2.mul(diff_flow2);
        
        std::vector<Mat1f> diff_flow2_channels(num_channels_flow);
        split(diff_flow2, diff_flow2_channels);
    
        sum_diff_flow2 = Mat1f::zeros(h, w);
        for (int c = 0; c < num_channels_flow; c++) sum_diff_flow2 = sum_diff_flow2 + diff_flow2_channels[c];
    }
    else if( is_disp_set )
    {
        Mat1f disp_XY_t0_warp_2_t1  = Mat1f::zeros(h,w);;
        remap(disp_t0, disp_XY_t0_warp_2_t1, map_t0_to_t1_x, map_t0_to_t1_y, cv::INTER_CUBIC);
        
        Mat1f diff_disp2 = disp_t1 - disp_XY_t0_warp_2_t1;
        sum_diff_flow2 = diff_disp2.mul(diff_disp2);
    }

    sqrt(sum_diff_flow2, sum_diff_flow2);
    
    pow(abs(sum_diff_flow2 / (sqrt_2 * sigma_grad)), alpha_grad, perm_gradient);
    perm_gradient = 1 / (1 + perm_gradient);

    perm_t = perm_photo.mul(perm_gradient);

    perm_t.convertTo(perm_t, CV_32FC1);
    is_perm_t_set = true;

    // FOR DEBUG PURPOSE
    // Mat2f prev_map = Mat2f::zeros(h,w);
    // for (int y = 0; y < h; y++) {
    //     for (int x = 0; x < w; x++){
    //         prev_map(y,x)[0] = x - flow_t0(y,x)[0];
    //         prev_map(y,x)[1] = y - flow_t0(y,x)[1];
    //     }
    // }
    // std::vector<Mat1f> prev_maps(num_channels_flow);
    // split(prev_map, prev_maps);
    // remap(I_prev, I_t0_warp_2_t1, prev_maps[0], prev_maps[1], cv::INTER_CUBIC);

    // Mat I_org, I_warp;
    // I.convertTo(I_org, CV_8UC3, 255);
    // I_t0_warp_2_t1.convertTo(I_warp, CV_8UC3, 255);
    
    // int rnum = rand();
    // std::ostringstream oss;
    // oss << rnum;

    // imwrite("original_image" + oss.str() + ".png", I_org);
    // imwrite("warp_image" + oss.str() + ".png", I_warp);

}



// For single-channel target images J
template <class TI>
Mat1f PermeabilityFilter<TI>::filterT(const Mat1f J_t0_XYT, const Mat1f J_t1_XY, Mat1f &l_t0_num, Mat1f &l_t0_den)
{
    if(! is_perm_t_set)
    {
        cerr << "Can not compute temporal permeability filtering, temporal permeability map is not computed." << endl;
        exit (EXIT_FAILURE);
    }
    if(! is_flow_set && ! is_disp_set)
    {
        cerr << "Can not compute temporal permeability filtering, optical flows / disparity maps are not defined." << endl;
        exit (EXIT_FAILURE);
    }
    if( is_flow_set && is_disp_set)
    {
        cerr << "Can not compute temporal permeability filtering, optical flows and disparity maps are both set." << endl;
        exit (EXIT_FAILURE);
    }

    // Input image
    int h = J_t1_XY.rows;
    int w = J_t1_XY.cols;
    int num_channels_flow = J_t1_XY.channels();

    Mat1f J_t1_XYT;
    J_t1_XY.copyTo(J_t1_XYT);

    Mat1f l_t1_num = Mat1f::zeros(h, w);
    Mat1f l_t1_den = Mat1f::zeros(h, w);

    for (int i = 0; i < iter_T; ++i)
    {
        // temporal filtering
        // no need to do pixel-wise operation (via J(y,x)) since all operation is based on same location pixels
        // forward pass & combining
        Mat1f temp_l_t0_num = l_t0_num + J_t0_XYT;
       
        Mat1f temp_l_t0_num_warp_2_t1 = Mat1f::zeros(h, w);
        //temp_l_t0_num = l_t0_num + J_t0_XYT;
        remap(temp_l_t0_num, temp_l_t0_num_warp_2_t1, map_t0_to_t1_x, map_t0_to_t1_y, cv::INTER_CUBIC);
        
        l_t1_num = perm_t.mul(temp_l_t0_num_warp_2_t1);

        
        Mat1f temp_l_t0_den = l_t0_den + 1.0;
        
        Mat1f temp_l_t0_den_warp_2_t1 = Mat1f::zeros(h, w);
        remap(temp_l_t0_den, temp_l_t0_den_warp_2_t1, map_t0_to_t1_x, map_t0_to_t1_y, cv::INTER_CUBIC);

        l_t1_den = perm_t.mul(temp_l_t0_den_warp_2_t1);

        J_t1_XYT = (l_t1_num + (1 - lambda_T) * J_t1_XY) / (l_t1_den + 1.0);
        
        l_t0_num = l_t1_num;
        l_t0_den = l_t1_den;
    }

    return J_t1_XYT;
}



// For multi-channel target images J
template <class TI>
template <class TJ>
Mat_<TJ> PermeabilityFilter<TI>::filterT(const Mat_<TJ> J_t0_XYT, const Mat_<TJ> J_t1_XY, Mat_<TJ> &l_t0_num, Mat_<TJ> &l_t0_den)
{
    if(! is_perm_t_set)
    {
        cerr << "Can not compute temporal permeability filtering, temporal permeability map is not computed." << endl;
        exit (EXIT_FAILURE);
    }
    if(! is_flow_set && ! is_disp_set)
    {
        cerr << "Can not compute temporal permeability filtering, optical flows / disparity maps are not defined." << endl;
        exit (EXIT_FAILURE);
    }
    if( is_flow_set && is_disp_set)
    {
        cerr << "Can not compute temporal permeability filtering, optical flows and disparity maps are both set." << endl;
        exit (EXIT_FAILURE);
    }
    // if(! is_l_init)
    // {
    //     cerr << "Can not compute temporal permeability filtering, left pass accumulated buffer are not initialized." << endl;
    //     exit (EXIT_FAILURE);
    // }

    // Input image
    int h = J_t1_XY.rows;
    int w = J_t1_XY.cols;
    int num_channels_flow = J_t1_XY.channels();

    Mat_<TJ> J_t1_XYT;
    J_t1_XY.copyTo(J_t1_XYT);

    Mat_<TJ> l_t1_num = Mat_<TJ>::zeros(h, w);
    Mat_<TJ> l_t1_den = Mat_<TJ>::zeros(h, w);

    for (int i = 0; i < iter_T; ++i)
    {
        // temporal filtering
        // no need to do pixel-wise operation (via J(y,x)) since all operation is based on same location pixels
        // forward pass & combining
        Mat_<TJ> temp_l_t0_num = Mat_<TJ>::zeros(h, w);
        Mat_<TJ> temp_l_t0_num_warp_2_t1 = Mat_<TJ>::zeros(h, w);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < num_channels_flow; c++) {
                    temp_l_t0_num(y,x)[c] = l_t0_num(y,x)[c] + J_t0_XYT(y,x)[c];
                }
            }
        }
        //temp_l_t0_num = l_t0_num + J_t0_XYT;
        remap(temp_l_t0_num, temp_l_t0_num_warp_2_t1, map_t0_to_t1_x, map_t0_to_t1_y, cv::INTER_CUBIC);
        
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < num_channels_flow; c++) {
                    l_t1_num(y,x)[c] = perm_t(y,x) * temp_l_t0_num_warp_2_t1(y,x)[c];
                }
            }
        }
        
        Mat_<TJ> temp_l_t0_den = Mat_<TJ>::zeros(h, w);
        Mat_<TJ> temp_l_t0_den_warp_2_t1 = Mat_<TJ>::zeros(h, w);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < num_channels_flow; c++) {
                    temp_l_t0_den(y,x)[c] = l_t0_den(y,x)[c] + 1.0;
                }
            }
        }
        
        remap(temp_l_t0_den, temp_l_t0_den_warp_2_t1, map_t0_to_t1_x, map_t0_to_t1_y, cv::INTER_CUBIC);
        
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < num_channels_flow; c++) {
                    l_t1_den(y,x)[c] = perm_t(y,x) * temp_l_t0_den_warp_2_t1(y,x)[c];
                    J_t1_XYT(y,x)[c] = (l_t1_num(y,x)[c] + (1 - lambda_T) * J_t1_XY(y,x)[c]) / (l_t1_den(y,x)[c]+ 1.0);
                }
            }
        }

        l_t0_num = l_t1_num;
        l_t0_den = l_t1_den;
    }

    return J_t1_XYT;
}


