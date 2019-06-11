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

#ifndef CPMPF_PARAMS
#define CPMPF_PARAMS

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "Variational_refinement/variational.h"

class cpmpf_parameters
{
public:
	/* Public Variables */
    // Output folders
    std::string output_CPM_dir;
    std::string output_PF_dir;
    std::string output_VR_dir;
    bool write_intermediate_results;

    // Addition naming options for input images
    int img_idx_width;
    int img_skip;
    std::string img_suf;

    // CPM parameters
    int CPM_max_displacement_input_int;
    int CPM_check_threshold_input_int;
    int CPM_cost_threshold_input_int;
    int CPM_stereo_flag;
    int CPM_step;

    // Permeability filter 
    // spatial parameters
    int PF_iter_XY;
    float PF_lambda_XY;
    float PF_sigma_XY;
    float PF_alpha_XY;

    // temporal parameters
    int PF_iter_T;
    float PF_lambda_T;
    float PF_sigma_photo;
    float PF_sigma_grad;
    float PF_alpha_photo;
    float PF_alpha_grad;

    // Variational refinement parameters
    float VR_alpha;
    float VR_gamma;
    float VR_delta;
    float VR_sigma;
    int VR_niter_outer;
    int VR_niter_inner;  
    int VR_niter_solver;
    float VR_sor_omega;

    /* Public Methods */
	// Constructor
    cpmpf_parameters();

    // Operator overloading
	friend std::ostream& operator<< (std::ostream& os, cpmpf_parameters& cpmpf_param); // For debug purpose mostly

    void set_dataset(std::string dataset_name);

    // Conversion to variational parameters
    void to_variational_params(variational_params_t *v_params);
};

#endif CPMPF_PARAMS
