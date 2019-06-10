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
    // CPM parameters
    int CPM_max_displacement_input_int;
    int CPM_check_threshold_input_int;
    int CPM_cost_threshold_input_int;
    int CPM_stereo_flag;
    int CPM_step;

    // Permeability filter parameters
    int PF_iterations_input_int;
    float PF_lambda_XY_input_float;
    float PF_delta_XY_input_float;
    float PF_alpha_XY_input_float;

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
