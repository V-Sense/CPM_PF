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

#include "CPM/CPM.h"
#include "PFilter/PermeabilityFilter.h"
extern "C" {
#include "Variational_refinement/variational.h"
}

#include "flow.h"
#include "ImageIOpfm.h"
#include "utils.h"
#include "cpmpf_parameters.h"



void Usage()
{
    std::cout<< endl
        << "Usage:" << endl
        << "  ./CPMPF <input_image_folder> <img_pre> <img_ext> <start_idx> <nb_imgs> <ang_dir> [options]" << endl
        << "Options:" << endl
        << "    -h, -help                                  print this message" << endl
        << "  Additional image naming options:" << endl
        << "    -img_idx_width                             length of the image index number" << endl
        << "    -img_skip                                  index number skip" << endl
        << "    -img_suf                                   suffix to add before image format extension" << endl
        << "  Output result folders:" << endl
        << "    -o, -output_VR                             set the final output folder (after variational refinement), default is <input_image_folder>" << endl
        << "    -save_intermediate                         use this flag to save results from CPM and PF steps, use the following flages to set the output folders" << endl
        << "    -output_CPM                                set the output folder for the Coarse-to-fine Patchmatch step, default is <input_image_folder>" << endl
        << "    -output_PF                                 set the output folder for the Permeability Filter steps (both spatial and temporal, default is <input_image_folder>" << endl
        << "  CPM parameters:" << endl
        << "    -CPM_max                                   outlier handling maxdisplacement threshold" << endl
        << "    -CPM_fbth                                  forward and backward consistency threshold" << endl
        << "    -CPM_cth                                   matching cost check threshold" <<endl
        << "    -CPM_stereo                                stereo flag" <<endl
        << "    -CPM_nstep                                 number of step giving the final result resolution" <<endl
        << "  PF:" << endl
        << "    Spatial parameters:" << endl
        << "    -PF_iter_XY                                number of iterations" << endl
        << "    -PF_lambda_XY                              lagrangian factor to balance fidelity to the input data" << endl
        << "    -PF_sigma_XY                               transition point of the edge-stopping function" << endl
        << "    -PF_alpha_XY                               falloff rate of the edge-stopping function" << endl
        << "    Temporal parameters:" << endl
        << "    -PF_iter_T                                 number of iterations" << endl
        << "    -PF_lambda_T                               lagrangian factor to balance fidelity to the input data" << endl
        << "    -PF_sigma_photo                            transition point of the edge-stopping function based on color consistency" << endl
        << "    -PF_alpha_photo                            falloff rate of the edge-stopping function based on color consistency" << endl
        << "    -PF_sigma_grad                             transition point of the edge-stopping function based on flow-gradient magnitude" << endl
        << "    -PF_alpha_grad                             falloff rate of the edge-stopping function based on flow-gradient magnitude" << endl
        << "  VR parameters:" << endl
        << "    -VR_alpha                                  smoothness weight" << endl
        << "    -VR_gamma                                  gradient constancy assumption weight" << endl
        << "    -VR_delta                                  color constancy assumption weight" << endl
        << "    -VR_sigma                                  presmoothing of the images" << endl
        << "    -VR_niter_outer                            number of outer fixed point iterations" << endl
        << "    -VR_niter_inner                            number of inner fixed point iterations" << endl
        << "    -VR_niter_solver                           number of solver iterations " << endl
        << "    -VR_sor_omega                              omega parameter of sor method" << endl
        << "  Predefined parameters:" << endl
        << "    -HCI                                       parameters for the HCI synthetic light field dataset" << endl
        << "    -Stanford                                  parameters for the Stanford gantry light field dataset" << endl
        << "    -TCH                                       parameters for the Technicolor camera array light field dataset" << endl
        << endl;
}

void parse_cmd(int argc, char** argv, int current_arg, cpmpf_parameters &cpm_pf_params) {

    #define isarg(key)  !strcmp(a,key)
    
    while(current_arg < argc){
        const char* a = argv[current_arg++];
        if( isarg("-h") || isarg("-help") )
            Usage();
        
        // Coarse-to-fine Patchmatch parameters
        else if( isarg("-CPM_max") )
            cpm_pf_params.CPM_max_displacement = atoi(argv[current_arg++]);
        else if( isarg("-CPM_fbth") )
            cpm_pf_params.CPM_check_threshold = atoi(argv[current_arg++]);
        else if( isarg("-CPM_cth") )
            cpm_pf_params.CPM_cost_threshold = atoi(argv[current_arg++]);
        else if( isarg("-CPM_stereo") )
            cpm_pf_params.CPM_stereo_flag = atoi(argv[current_arg++]);
        else if( isarg("-CPM_nstep") )
            cpm_pf_params.CPM_step = atoi(argv[current_arg++]);
        
        // Permeability Filter 
        // spatial parameters
        else if( isarg("-PF_iter_XY") )
            cpm_pf_params.PF_iter_XY = atof(argv[current_arg++]);
        else if( isarg("-PF_lambda_XY") )
            cpm_pf_params.PF_lambda_XY = atof(argv[current_arg++]);
        else if( isarg("-PF_sigma_XY") )
            cpm_pf_params.PF_sigma_XY = atof(argv[current_arg++]);
        else if( isarg("-PF_alpha_XY") )
            cpm_pf_params.PF_alpha_XY = atof(argv[current_arg++]);
        // temporal parameters
        else if( isarg("-PF_iter_T") )
            cpm_pf_params.PF_iter_T = atof(argv[current_arg++]);
        else if( isarg("-PF_lambda_T") )
            cpm_pf_params.PF_lambda_T = atof(argv[current_arg++]);
        else if( isarg("-PF_sigma_photo") )
            cpm_pf_params.PF_sigma_photo = atof(argv[current_arg++]);
        else if( isarg("-PF_sigma_grad") )
            cpm_pf_params.PF_sigma_grad = atof(argv[current_arg++]);
        else if( isarg("-PF_alpha_photo") )
            cpm_pf_params.PF_alpha_photo = atof(argv[current_arg++]);
        else if( isarg("-PF_alpha_grad") )
            cpm_pf_params.PF_alpha_grad = atof(argv[current_arg++]);

        // Variational refinement parameters
        else if( isarg("-VR_alpha") )
            cpm_pf_params.VR_alpha = atof(argv[current_arg++]);
        else if( isarg("-VR_gamma") )
            cpm_pf_params.VR_gamma = atof(argv[current_arg++]);
        else if( isarg("-VR_delta") )
            cpm_pf_params.VR_delta = atof(argv[current_arg++]);
        else if( isarg("-VR_sigma") )
            cpm_pf_params.VR_sigma = atof(argv[current_arg++]);
        else if( isarg("-VR_niter_outer") )
            cpm_pf_params.VR_niter_outer = atof(argv[current_arg++]);
        else if( isarg("-VR_niter_inner") )
            cpm_pf_params.VR_niter_inner = atof(argv[current_arg++]);
        else if( isarg("-VR_niter_solver") )
            cpm_pf_params.VR_niter_solver = atof(argv[current_arg++]);
        else if( isarg("-VR_sor_omega") )
            cpm_pf_params.VR_sor_omega = atof(argv[current_arg++]);

        
        // Predefined parameters for common test datasets
        // Video dataset for optical flow estimation
        else if( isarg("-Sintel") )
            cpm_pf_params.set_dataset(string("Sintel"));
        
        // Various light field datasets for disparity estimation
        // Each dataset has a different disparity range with a different order of magnitude
        else if( isarg("-HCI") )
            cpm_pf_params.set_dataset(string("HCI"));
        else if( isarg("-Stanford") )
            cpm_pf_params.set_dataset(string("Stanford"));
        else if( isarg("-TCH") ) 
            cpm_pf_params.set_dataset(string("TCH"));

        // Output folders
        // Final results
        else if( isarg("-o") || isarg("-output_VR") )
            cpm_pf_params.output_VR_dir = string(argv[current_arg++]);
        // Intermediate results
        else if( isarg("-save_intermediate") )
            cpm_pf_params.write_intermediate_results = true;
        else if( isarg("-output_CPM") )
            cpm_pf_params.output_CPM_dir = string(argv[current_arg++]);
        else if( isarg("-output_PF") )
            cpm_pf_params.output_PF_dir = string(argv[current_arg++]);

        // Additional options for input image naming
        else if( isarg("-img_idx_width") ) // Number of zeros in the index
            cpm_pf_params.img_idx_width = atoi(argv[current_arg++]);
        else if( isarg("-img_skip") )
            cpm_pf_params.img_skip = atoi(argv[current_arg++]);
        else if( isarg("-img_suf") )
            cpm_pf_params.img_suf = string(argv[current_arg++]);            


        else {
            fprintf(stderr, "unknown argument %s\n", a);
            Usage();
            exit(1);
        }
    }
}

int main(int argc, char** argv)
{
    if (argc < 7){
        if (argc > 1) fprintf(stderr, "Error: not enough arguments\n");
        Usage();
        exit(1);
	}

    CTimer total_time;

    /* ---------------- PARSE COMMAND LINE --------------------------- */
    // Mandatory parameters
    int current_arg = 1;
    string input_images_folder = string(argv[current_arg++]);
    string img_pre = string(argv[current_arg++]);
    if(img_pre == "none") img_pre = "";
    string img_ext = string(argv[current_arg++]);
    int start_idx = atoi(argv[current_arg++]);
    int nb_imgs = atoi(argv[current_arg++]);
    string ang_dir = string(argv[current_arg++]);

    cpmpf_parameters cpm_pf_params("HCI"); // Initiates default params
    
    // Results folder, default is the input folder
    cpm_pf_params.output_CPM_dir = input_images_folder;
    cpm_pf_params.output_PF_dir  = input_images_folder;
    cpm_pf_params.output_VR_dir  = input_images_folder;
    
    // Optional parameters
    parse_cmd(argc, argv, current_arg, cpm_pf_params);

    vector<string> input_images_name_vec(nb_imgs);

    
    /* ---------------- READ INPUT RBG IMAGES --------------------------- */
    CTimer CPM_input_time;
    std::cout << "Reading input RGB images... " << endl;
    vector<Mat3f> input_RGB_images_vec(nb_imgs);
    int width = -1 , height = -1, nch = -1;

    for (size_t i = 0; i < nb_imgs; i++) {
        std::stringstream ss_idx;
        ss_idx << std::setw(cpm_pf_params.img_idx_width) << std::setfill('0') << start_idx + i * cpm_pf_params.img_skip;
		std::string s_idx = ss_idx.str();
        String img_name = img_pre + s_idx + cpm_pf_params.img_suf + img_ext;
        input_images_name_vec[i] = img_name;
        
        std::cout << input_images_folder + "/" + img_name << endl;
        Mat tmp_img = imread(input_images_folder + "/" + img_name);
        if ( tmp_img.empty() ) {
            std::cout << input_images_name_vec[i] << " is invalid!" << endl;
            continue;
        }

        if( i == 0 ) {
            width  = tmp_img.cols;
            height = tmp_img.rows;
            nch = tmp_img.channels();
        }
        else
        {
            if(width != tmp_img.cols || height != tmp_img.rows || nch != tmp_img.channels()) {
                std::cout << "All input images should have the same size; size of image " << input_images_name_vec[i] << " is different from previous image" << endl;
                return -1;
            }
        }
        
        Mat3f tmp_img3f;
        tmp_img.convertTo(tmp_img3f, CV_32F, 1/255.);

        if(ang_dir == "ver") // Rotate image 90 degress and process them as horizontal parallax (allows to use stereo_flag=1 for CPM)
        {
            cv::rotate(tmp_img3f, tmp_img3f, cv::ROTATE_90_COUNTERCLOCKWISE);
        }

        input_RGB_images_vec[i] = tmp_img3f;
    }
    CPM_input_time.toc(" done in: ");


    /* ---------------- RUN COARSE-TO-FINE PATCHMATCH --------------------------- */
    CTimer CPM_time;
    std::cout << "Running CPM... " << flush;
    
    CPM cpm;
    cpm_pf_params.to_CPM_params(cpm);

    vector<Mat1f> cpm_disp_fwd(nb_imgs-1), cpm_disp_bwd(nb_imgs-1);
    
    #pragma omp parallel for 
    for (size_t i = 0; i < nb_imgs - 1; ++i) {
        FImage img1(width, height, nch);
        FImage img2(width, height, nch);
        
        Mat3f2FImage(input_RGB_images_vec[i],   img1);
        Mat3f2FImage(input_RGB_images_vec[i+1], img2);

        // Forward flow
        FImage matches;
        cpm.Matching(img1, img2, matches);

        Mat1f disp_fwd(height, width, kMOVEMENT_UNKNOWN);
        Match2Disp(matches, disp_fwd, "hor");
        cpm_disp_fwd[i] = disp_fwd;

        // Backward flow
        matches.clear();
        cpm.Matching(img2, img1, matches);

        Mat1f disp_bwd(height, width, kMOVEMENT_UNKNOWN);
        Match2Disp(matches, disp_bwd, "hor");
        cpm_disp_bwd[i] = disp_bwd;
    }
    CPM_time.toc(" done in: ");


    // Write CPM results on disk
    if(cpm_pf_params.write_intermediate_results) {
        CTimer CPM_write_time;
        std::cout << "Writing CPM results... " << flush;
        #pragma omp parallel for 
        for (size_t i = 0; i < nb_imgs - 1; ++i) {
            // Remove image file extension
            string img_name1 = input_images_name_vec[i];
            string img_name2 = input_images_name_vec[i+1];
            str_replace(img_name1, img_ext, "");
            str_replace(img_name2, img_ext, "");

            // Forward matching
            Mat1f cpm_disp = cpm_disp_fwd[i];
            Mat1f mask_flow_unknown = cpm_disp != kMOVEMENT_UNKNOWN;
            cpm_disp = cpm_disp.mul(mask_flow_unknown); // Set unkown flow to 0 before writing to pfm file

            if(ang_dir == "ver") // Rotate back to original orientation
                cv::rotate(-cpm_disp, cpm_disp, cv::ROTATE_90_CLOCKWISE);
            
            string disp_file = cpm_pf_params.output_CPM_dir +  "/CPM__" + img_name1 + "__TO__" + img_name2 + ".pfm";
            WriteFilePFM(-cpm_disp, disp_file.c_str(), 1/255.0);
            
            // Backward matching
            Mat1f cpm_disp = cpm_disp_bwd[i];
            Mat1f mask_flow_unknown = cpm_disp != kMOVEMENT_UNKNOWN;
            cpm_disp = cpm_disp.mul(mask_flow_unknown); // Set unkown flow to 0 before writing to pfm file

            if(ang_dir == "ver") // Rotate back to original orientation
                cv::rotate(-cpm_disp, cpm_disp, cv::ROTATE_90_CLOCKWISE);

            disp_file = cpm_pf_params.output_CPM_dir +  "/CPM__" + img_name2 + "__TO__" + img_name1 + ".pfm";
            WriteFilePFM(-cpm_disp, disp_file.c_str(), 1/255.0);
        }
        CPM_write_time.toc(" done in: ");
    }


    /* ---------------- RUN PERMEABILITY FILTER --------------------------- */
    PermeabilityFilter<Vec3f> PF;
    cpm_pf_params.to_PF_params<Vec3f>(PF);

    // spatial filter
    CTimer sPF_time;
    std::cout << "Running spatial permeability filter... " << flush;
    vector<Mat1f> pf_spatial_flow_vec(nb_imgs);

    for (size_t i = 0; i < nb_imgs; ++i) {
        PF.set_I_XY(input_RGB_images_vec[i]); // Set guide image
        Mat1f disp_forward, disp_backward, disp_confidence;
        if(i == (nb_imgs-1)) // For the last image, associate the backward flow with a minus sign
        {
            disp_forward  = cpm_disp_fwd[i-1];
            disp_backward = cpm_disp_bwd[i-1];

            // compute flow confidence map
            disp_confidence = getHorDispConfidence(disp_backward, disp_forward);
        }
        else
        {
            disp_forward  = cpm_disp_fwd[i];
            disp_backward = cpm_disp_bwd[i];

            // compute flow confidence map
            disp_confidence = getHorDispConfidence(disp_forward, disp_backward);
        }
        
        
        // Apply spatial permeability filter on confidence
        PF.computeSpatialPermeabilityMaps();
        Mat1f disp_confidence_filtered = PF.filterXY(disp_confidence);

        // multiply initial confidence and sparse flow
        Mat1f confidenced_disp;
        if(i == (nb_imgs-1)) // For the last image, associate the backward flow with a minus sign
        {
            confidenced_disp = disp_backward.mul(-disp_confidence);    
        }
        else
        {
            confidenced_disp = disp_forward.mul(disp_confidence);    
        }
            
        
        // filter confidenced sparse flow
        Mat1f confidenced_disp_XY = PF.filterXY(confidenced_disp);

        // compute normalized spatial filtered flow FXY by division
        Mat1f normalized_confidenced_disp_filtered = confidenced_disp_XY.mul(1 / disp_confidence_filtered);
        
        pf_spatial_flow_vec[i] = normalized_confidenced_disp_filtered;
    }

    sPF_time.toc(" done in: ");


    // Write spatial PF results on disk
    if(cpm_pf_params.write_intermediate_results) {
        CTimer sPF_write_time;
        std::cout << "Writing spatial permeability filter results... " << flush;
        #pragma omp parallel for 
        for (size_t i = 0; i < nb_imgs; ++i) {
            // Remove image file extension
            string img_name1, img_name2; 
            if( i == (nb_imgs-1)) {
                img_name1 = input_images_name_vec[i];
                img_name2 = input_images_name_vec[i-1];
            }
            else {
                img_name1 = input_images_name_vec[i];
                img_name2 = input_images_name_vec[i+1];
            }
            str_replace(img_name1, img_ext, "");
            str_replace(img_name2, img_ext, "");

            Mat1f pf_disp = pf_spatial_flow_vec[i];
            
            if(ang_dir == "ver") // Rotate back to original orientation
                cv::rotate(-pf_disp, pf_disp, cv::ROTATE_90_CLOCKWISE);

            string disp_file = cpm_pf_params.output_PF_dir +  "/PF_spatial__" + img_name1 + "__TO__" + img_name2 + ".pfm";
            WriteFilePFM(-pf_disp, disp_file.c_str(), 1/255.0);
        }
        sPF_write_time.toc(" done in: ");
    }

    // temporal filter
    CTimer tPF_time;
    std::cout << "Running temporal permeability filter... " << flush;
    Mat2f l_prev = Mat2f::zeros(height, width);
    Mat2f l_normal_prev = Mat2f::zeros(height, width);
    Mat2f disp_t0_XYT, disp_t1_XYT;
    vector<Mat2f> disp_t1_XYT_vector;
    vector<Mat1f> pf_temporal_flow_vec(nb_imgs);

    pf_temporal_flow_vec[0] = pf_spatial_flow_vec[0];

    for (size_t i = 1; i < nb_imgs; ++i)
    {
        PF.set_I_T(input_RGB_images_vec[i-1], input_RGB_images_vec[i]);
        Mat2f disp_t0_XY = pf_spatial_flow_vec[i - 1];
        Mat2f disp_t1_XY = pf_spatial_flow_vec[i];

        if(i == 1) {
            PF.computeTemporalPermeability(disp_t0_XY, disp_t1_XY);
            disp_t1_XYT_vector = PF.filterT<Vec2f>(disp_t1_XY, disp_t0_XY, disp_t0_XY,  l_prev, l_normal_prev);
        }
        else {
            PF.computeTemporalPermeability(disp_t0_XYT, disp_t1_XY);
            disp_t1_XYT_vector = PF.filterT<Vec2f>(disp_t0_XY, disp_t1_XY, disp_t0_XYT, l_prev, l_normal_prev);
        }

        disp_t1_XYT = disp_t1_XYT_vector[2];
        l_prev = disp_t1_XYT_vector[0];
        l_normal_prev = disp_t1_XYT_vector[1];

        disp_t0_XYT = disp_t1_XYT;

        pf_temporal_flow_vec[i] = disp_t1_XYT;
    }
    tPF_time.toc(" done in: ");


    // Write temporal PF results on disk
    if(cpm_pf_params.write_intermediate_results) {
        CTimer tPF_write_time;
        std::cout << "Writing temporal permeability filter results... " << flush;
        #pragma omp parallel for 
        for (size_t i = 0; i < nb_imgs; ++i) {
            // Remove image file extension
            string img_name1, img_name2;
            if( i == (nb_imgs-1)) {
                img_name1 = input_images_name_vec[i];
                img_name2 = input_images_name_vec[i-1];
            }
            else {
                img_name1 = input_images_name_vec[i];
                img_name2 = input_images_name_vec[i+1];
            }
            str_replace(img_name1, img_ext, "");
            str_replace(img_name2, img_ext, "");

            // // Forward matching only
            // string flow_fwd_name = "PF_temporal__" + img_name1 + "__TO__" + img_name2;
            // string flow_fwd_file = cpm_pf_params.output_PF_dir + "/" + flow_fwd_name + ".flo";
            // WriteFlowFile(pf_temporal_flow_vec[i], flow_fwd_file.c_str());
            
            // flow_fwd_file = cpm_pf_params.output_PF_dir + "/" + flow_fwd_name + ".png";
            // WriteFlowAsImage(pf_temporal_flow_vec[i], flow_fwd_file.c_str(), -1);

            // Convert to disparity
            vector<Mat1f> pf_temporal_flow_split;
            cv::split(pf_temporal_flow_vec[i], pf_temporal_flow_split);

            if(ang_dir == "ver") // Rotate back to original orientation
                cv::rotate(-pf_temporal_flow_split[0], pf_temporal_flow_split[0], cv::ROTATE_90_CLOCKWISE);

            string disp_file = cpm_pf_params.output_PF_dir +  "/PF_temporal__" + img_name1 + "__TO__" + img_name2 + ".pfm";
            WriteFilePFM(-pf_temporal_flow_split[0], disp_file.c_str(), 1/255.0);
        }
        tPF_write_time.toc(" done in: ");
    }


    /* ---------------- RUN VARIATIONAL REFINEMENT  --------------------------- */
    CTimer var_time;
    std::cout << "Running variational refinement... " << flush;
    color_image_t *im1 = color_image_new(width, height);
    color_image_t *im2 = color_image_new(width, height);
    vector<Mat2f> vr_flow_vec(nb_imgs);
    Mat2f pf_flow;

    for (size_t i = 0; i < nb_imgs ; ++i) {
        if(i == (nb_imgs-1)){
            Mat3f2color_image_t(input_RGB_images_vec[i], im1);
            Mat3f2color_image_t(input_RGB_images_vec[i - 1], im2);
            pf_flow = -pf_temporal_flow_vec[i];
        }
        else {
            Mat3f2color_image_t(input_RGB_images_vec[i], im1);
            Mat3f2color_image_t(input_RGB_images_vec[i + 1], im2);
            pf_flow = pf_temporal_flow_vec[i];
        }
        
        variational_params_t flow_params;
        cpm_pf_params.to_variational_params(&flow_params);
        image_t *flow_x = image_new(width, height), *flow_y = image_new(width, height);
        Mat2f2image_t_uv(pf_flow, flow_x, flow_y);

        variational(flow_x, flow_y, im1, im2, &flow_params);

        Mat2f vr_flow(height, width);
        image_t_uv2Mat2f(vr_flow, flow_x, flow_y);

        vr_flow_vec[i] = vr_flow;
    }
    color_image_delete(im1);
    color_image_delete(im2);
    var_time.toc(" done in: ");


    // Write variational refinement results on disk
    CTimer vr_write_time;
    std::cout << "Writing variational refinement results... " << flush;
    #pragma omp parallel for 
    for (size_t i = 0; i < nb_imgs; ++i) {
        // Remove image file extension
        string img_name1, img_name2;
        if( i == (nb_imgs-1)) {
            img_name1 = input_images_name_vec[i];
            img_name2 = input_images_name_vec[i-1];
            vr_flow_vec[i] = -vr_flow_vec[i]; // Forces all flows to have the same direction
        }
        else {
            img_name1 = input_images_name_vec[i];
            img_name2 = input_images_name_vec[i+1];
        }
        str_replace(img_name1, img_ext, "");
        str_replace(img_name2, img_ext, "");

        // // Forward matching only
        // string flow_fwd_name = "VR__" + img_name1 + "__TO__" + img_name2;
        // string flow_fwd_file = cpm_pf_params.output_VR_dir + "/" + flow_fwd_name + ".flo";
        // WriteFlowFile(vr_flow_vec[i], flow_fwd_file.c_str());
        
        // flow_fwd_file = cpm_pf_params.output_VR_dir + "/" + flow_fwd_name + ".png";
        // WriteFlowAsImage(vr_flow_vec[i], flow_fwd_file.c_str(), -1);

        // Convert to disparity
        vector<Mat1f> vr_flow_split;
        cv::split(vr_flow_vec[i], vr_flow_split);

        if(ang_dir == "ver") // Rotate back to original orientation
                cv::rotate(-vr_flow_split[0], vr_flow_split[0], cv::ROTATE_90_CLOCKWISE);

        string disp_file = cpm_pf_params.output_VR_dir +  "/VR__" + img_name1 + "__TO__" + img_name2 + ".pfm";
        WriteFilePFM(-vr_flow_split[0], disp_file.c_str(), 1/255.0);
    }
    vr_write_time.toc(" done in: ");
    
    total_time.toc("\nTotal elapsed time: ");

    return EXIT_SUCCESS;
}
