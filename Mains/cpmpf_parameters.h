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

#ifndef CPMPF_PARAMS_H
#define CPMPF_PARAMS_H

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "CPM/CPM.h"
#include "PFilter/PermeabilityFilter.h"
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
    bool write_color_png;

    // Addition naming options for input images
    int img_idx_width;
    int img_skip;
    std::string img_suf;

    // CPM parameters
    int CPM_max_displacement;
    float CPM_check_threshold;
    int CPM_cost_threshold;
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
    cpmpf_parameters(std::string dataset_name);
    
    // Operator overloading
	friend std::ostream& operator<< (std::ostream& os, cpmpf_parameters& cpmpf_param); // For debug purpose mostly

    void set_dataset(std::string dataset_name);

    // Conversion to specific steps parameters
    void to_CPM_params(CPM &cpm);
    template <class TI>
    void to_PF_params(PermeabilityFilter<TI> &PF);
    void to_variational_params(variational_params_t *v_params);
};

// Conversion to specific steps parameters
// Define in .h because of template
template <class TI>
void cpmpf_parameters::to_PF_params(PermeabilityFilter<TI> &PF){
    // Spatial parameters
    PF.iter_XY = PF_iter_XY;
    PF.lambda_XY = PF_lambda_XY;
    PF.sigma_XY = PF_sigma_XY;
    PF.alpha_XY = PF_alpha_XY;

    // Temporal parameters
    PF.iter_T = PF_iter_T;
    PF.lambda_T = PF_lambda_T;
    PF.sigma_photo = PF_sigma_photo;
    PF.sigma_grad = PF_sigma_grad;
    PF.alpha_photo = PF_alpha_photo;
    PF.alpha_grad = PF_alpha_grad;
}

#endif //! CPMPF_PARAMS_H
