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

#include "cpmpf_parameters.h"

// default constructor initialized with Sintel default parameters
cpmpf_parameters::cpmpf_parameters() {
    // Set default parameters
    set_dataset(std::string("Sintel"));

    img_idx_width = 4;
    img_skip = 1;
    img_suf = "";

    // Set output directories to local directory
    output_CPM_dir = std::string("./");
    output_PF_dir  = std::string("./");
    output_VR_dir  = std::string("./");
    write_intermediate_results = false;
}

cpmpf_parameters::cpmpf_parameters(std::string dataset_name) {
    // Set default parameters
    set_dataset(dataset_name);

    img_idx_width = 4;
    img_skip = 1;
    img_suf = "";

    // Set output directories to local directory
    output_CPM_dir = std::string("./");
    output_PF_dir  = std::string("./");
    output_VR_dir  = std::string("./");
    write_intermediate_results = false;
}

// Operators overloading
// For debug purpose mostly
std::ostream& operator<< (std::ostream& os, cpmpf_parameters &cpmpf_param)
{
	os << std::endl << "[CPMPF parameters]" << std::endl;
	os << "CPM_max_displacement: "    << cpmpf_param.CPM_max_displacement << std::endl;
	os << "CPM_check_threshold: "     << cpmpf_param.CPM_check_threshold << std::endl;
	os << "CPM_stereo_flag: "         << cpmpf_param.CPM_stereo_flag << std::endl;
	os << "CPM_step: "                << cpmpf_param.CPM_step << std::endl;
	
	os << "PF_iter_XY: "    << cpmpf_param.PF_iter_XY << std::endl;
	os << "PF_lambda_XY: "  << cpmpf_param.PF_lambda_XY << std::endl;
	os << "PF_sigma_XY: "   << cpmpf_param.PF_sigma_XY << std::endl;
	os << "PF_alpha_XY: "   << cpmpf_param.PF_alpha_XY << std::endl;
    
	os << "VR_alpha: "          << cpmpf_param.VR_alpha << std::endl;
	os << "VR_gamma: "          << cpmpf_param.VR_gamma << std::endl;
	os << "VR_delta: "          << cpmpf_param.VR_delta << std::endl;
	os << "VR_sigma: "          << cpmpf_param.VR_sigma << std::endl;
    os << "VR_niter_outer: "    << cpmpf_param.VR_niter_outer << std::endl;
	os << "VR_niter_inner: "    << cpmpf_param.VR_niter_inner << std::endl;
	os << "VR_niter_solver: "   << cpmpf_param.VR_niter_solver << std::endl;
	os << "VR_sor_omega: "      << cpmpf_param.VR_sor_omega << std::endl;

	return os;
}

// Initialization parameters for some common datasets
void cpmpf_parameters::set_dataset(std::string dataset_name) {
    
    // First init with default parameters used for Sintel dataset
    CPM_max_displacement = 400;
    CPM_check_threshold = 1;
    CPM_cost_threshold = 1880;
    CPM_stereo_flag = 0;
    CPM_step = 3;

    PF_iter_XY = 5;
    PF_lambda_XY = 0;
    PF_sigma_XY = 0.017;
    PF_alpha_XY = 2;
    PF_iter_T = 1;
    PF_lambda_T = 0.0;
    PF_sigma_photo = 0.3;
    PF_sigma_grad = 1.0;
    PF_alpha_photo = 2.0;
    PF_alpha_grad = 2.0;

    VR_alpha = 1.0f;
    VR_gamma = 0.71f;
    VR_delta = 0.0f;
    VR_sigma = 1.00f;
    VR_niter_outer = 5;
    VR_niter_inner = 1;
    VR_niter_solver = 30;
    VR_sor_omega  = 1.9f;

    if(dataset_name == "Sintel")
    {} // Nothing to do as default parameters are defined from this dataset
    else if(dataset_name == "HCI"){
        CPM_max_displacement = 4;
        CPM_stereo_flag = 1;
    }
    else if( dataset_name == "Stanford" ) {
        CPM_max_displacement = 40;
        CPM_stereo_flag = 1;
    }
    else if( dataset_name == "TCH" ) {
        CPM_max_displacement = 400;
        CPM_stereo_flag = 1;
    }
    else {
        std::cout << "Dataset name unknown for CPMPF parameter initialization (Sintel dataset parameter used)."  << std::endl;
    }
}

// Conversion to variational parameters
void cpmpf_parameters::to_CPM_params(CPM &cpm){
    cpm.SetStereoFlag(CPM_stereo_flag);
	cpm.SetStep(CPM_step);
	cpm.SetMaxDisplacement(CPM_max_displacement);
	cpm.SetCheckThreshold(CPM_check_threshold);
}

void cpmpf_parameters::to_variational_params(variational_params_t *v_params){
    if(!v_params){
        fprintf(stderr,"Error optical_flow_params_default: argument is null\n");
        exit(1);
    }
    v_params->alpha = VR_alpha;
    v_params->gamma = VR_gamma;
    v_params->delta = VR_delta;
    v_params->sigma = VR_sigma;
    v_params->niter_outer = VR_niter_outer;
    v_params->niter_inner = VR_niter_inner;  
    v_params->niter_solver = VR_niter_solver;
    v_params->sor_omega = VR_sor_omega;
}
