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
#include "utils.h"
#include "cpmpf_parameters.h"



void Usage()
{
    std::cout<< endl
        << "Usage:" << endl
        << "  ./CPMPF <input_image_folder> <img_pre> <img_ext> <start_idx> <nb_imgs> [options]" << endl
        << "Options:" << endl
        << "    -h, -help                                  print this message" << endl
        << "  Additional image naming options:" << endl
        << "    -img_idx_width                             length of the image index number" << endl
        << "    -img_skip                                  index number skip" << endl
        << "    -img_suf                                   suffix to add before image format extension" << endl
        << "  Output result folders:" << endl
        << "    -o, -output_VR                             set the final output folder (after variational refinement), default is <input_image_folder>" << endl
        << "    -save_intermediate                         use this flag to save results from CPM and PF steps, use the following flages to set the output folders" << endl
        << "    -write_color_png                           write results as color png files using optical flow convention" << endl
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
        << "    -Sintel                                    parameters for the MPI-Sintel dataset" << endl
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
        else if( isarg("-write_color_png") )
            cpm_pf_params.write_color_png = true;
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
    if (argc < 6){
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

    cpmpf_parameters cpm_pf_params("Sintel"); // Initiates default params
    
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
        input_RGB_images_vec[i] = tmp_img3f;
    }
    CPM_input_time.toc(" done in: ");


    /* ---------------- RUN COARSE-TO-FINE PATCHMATCH --------------------------- */
    CTimer CPM_time;
    std::cout << "Running CPM... " << flush;
    
    CPM cpm;
    cpm_pf_params.to_CPM_params(cpm);

    vector<Mat2f> cpm_flow_fwd(nb_imgs-1), cpm_flow_bwd(nb_imgs-1);
    
    #pragma omp parallel for 
    for (size_t i = 0; i < nb_imgs - 1; ++i) {
        FImage img1(width, height, nch);
        FImage img2(width, height, nch);
        
        Mat3f2FImage(input_RGB_images_vec[i],   img1);
        Mat3f2FImage(input_RGB_images_vec[i+1], img2);

        // Forward flow
        FImage matches;
        cpm.Matching(img1, img2, matches);

        Mat2f flow_fwd(height, width, kMOVEMENT_UNKNOWN);
        Match2Flow(matches, flow_fwd);
        cpm_flow_fwd[i] = flow_fwd;

        // Backward flow
        matches.clear();
        cpm.Matching(img2, img1, matches);

        Mat2f flow_bwd(height, width, kMOVEMENT_UNKNOWN);
        Match2Flow(matches, flow_bwd);
        cpm_flow_bwd[i] = flow_bwd;
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
            string flow_fwd_name = "CPM__" + img_name1 + "__TO__" + img_name2;
            string flow_fwd_file = cpm_pf_params.output_CPM_dir + "/" + flow_fwd_name + ".flo";
            WriteFlowFile(cpm_flow_fwd[i], flow_fwd_file.c_str());
            
            if(cpm_pf_params.write_color_png) {
                flow_fwd_file = cpm_pf_params.output_CPM_dir + "/" + flow_fwd_name + ".png";
                WriteFlowAsImage(cpm_flow_fwd[i], flow_fwd_file.c_str(), -1);
            }


            // Backward matching
            string flow_bwd_name = "CPM__" + img_name2 + "__TO__" + img_name1;
            string flow_bwd_file = cpm_pf_params.output_CPM_dir + "/" + flow_bwd_name + ".flo";
            WriteFlowFile(cpm_flow_bwd[i], flow_bwd_file.c_str());
            
            if(cpm_pf_params.write_color_png) {
                flow_bwd_file = cpm_pf_params.output_CPM_dir + "/" + flow_bwd_name + ".png";
                WriteFlowAsImage(cpm_flow_bwd[i], flow_bwd_file.c_str(), -1);
            }
        }
        CPM_write_time.toc(" done in: ");
    }


    /* ---------------- RUN PERMEABILITY FILTER --------------------------- */
    PermeabilityFilter<Vec3f> PF;
    cpm_pf_params.to_PF_params<Vec3f>(PF);

    // spatial filter
    CTimer sPF_time;
    std::cout << "Running spatial permeability filter... " << flush;
    vector<Mat2f> pf_spatial_flow_vec(nb_imgs);

    for (size_t i = 0; i < nb_imgs - 1; ++i) {
        PF.set_I_XY(input_RGB_images_vec[i]); // Set guide image
        Mat2f flow_forward  = cpm_flow_fwd[i];
        Mat2f flow_backward = cpm_flow_bwd[i];

        // compute flow confidence map
        Mat1f flow_confidence = getFlowConfidence(flow_forward, flow_backward);

        // Apply spatial permeability filter on confidence
        PF.computeSpatialPermeabilityMaps();
        Mat1f flow_confidence_filtered = PF.filterXY(flow_confidence);


        // multiply initial confidence and sparse flow
        Mat2f confidenced_flow = Mat2f::zeros(flow_confidence.rows,flow_confidence.cols);
        for(int y = 0; y < confidenced_flow.rows; y++) {
            for(int x = 0; x < confidenced_flow.cols; x++) {
                for(int c = 0; c < confidenced_flow.channels(); c++) {
                    confidenced_flow(y,x)[c] = flow_forward(y,x)[c] * flow_confidence(y,x);
                }
            }
        }

        //filter confidenced sparse flow
        Mat2f confidenced_flow_XY = PF.filterXY<Vec2f>(confidenced_flow);

        // compute normalized spatial filtered flow FXY by division
        Mat2f normalized_confidenced_flow_filtered = Mat2f::zeros(confidenced_flow_XY.rows,confidenced_flow_XY.cols);
        for(int y = 0; y < confidenced_flow_XY.rows; y++) {
            for(int x = 0; x < confidenced_flow_XY.cols; x++) {
                for(int c = 0; c < confidenced_flow_XY.channels(); c++) {
                    normalized_confidenced_flow_filtered(y,x)[c] = confidenced_flow_XY(y,x)[c] / flow_confidence_filtered(y,x);
                }
            }
        }

        pf_spatial_flow_vec[i] = normalized_confidenced_flow_filtered;
    }

    // For the last image, associate the backward flow with a minus sign
    PF.set_I_XY(input_RGB_images_vec[nb_imgs - 1]); // Set guide image
    Mat2f flow_forward  = cpm_flow_fwd[nb_imgs - 2];
    Mat2f flow_backward = cpm_flow_bwd[nb_imgs - 2];

    // compute flow confidence map
    Mat1f flow_confidence = getFlowConfidence(flow_backward, flow_forward);

    // Apply spatial permeability filter on confidence
    PF.computeSpatialPermeabilityMaps();
    Mat1f flow_confidence_filtered = PF.filterXY(flow_confidence);


    // multiply initial confidence and sparse flow
    Mat2f confidenced_flow = Mat2f::zeros(flow_confidence.rows,flow_confidence.cols);
    for(int y = 0; y < confidenced_flow.rows; y++) {
        for(int x = 0; x < confidenced_flow.cols; x++) {
            for(int c = 0; c < confidenced_flow.channels(); c++) {
                confidenced_flow(y,x)[c] = -flow_backward(y,x)[c] * flow_confidence(y,x);
            }
        }
    }

    //filter confidenced sparse flow
    Mat2f confidenced_flow_XY = PF.filterXY<Vec2f>(confidenced_flow);

    // compute normalized spatial filtered flow FXY by division
    Mat2f normalized_confidenced_flow_filtered = Mat2f::zeros(confidenced_flow_XY.rows,confidenced_flow_XY.cols);
    for(int y = 0; y < confidenced_flow_XY.rows; y++) {
        for(int x = 0; x < confidenced_flow_XY.cols; x++) {
            for(int c = 0; c < confidenced_flow_XY.channels(); c++) {
                normalized_confidenced_flow_filtered(y,x)[c] = confidenced_flow_XY(y,x)[c] / flow_confidence_filtered(y,x);
            }
        }
    }

    pf_spatial_flow_vec[nb_imgs - 1] = normalized_confidenced_flow_filtered;
    

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

            // Forward matching only
            string flow_fwd_name = "PF_spatial__" + img_name1 + "__TO__" + img_name2;
            string flow_fwd_file = cpm_pf_params.output_PF_dir + "/" + flow_fwd_name + ".flo";
            WriteFlowFile(pf_spatial_flow_vec[i], flow_fwd_file.c_str());
            
            if(cpm_pf_params.write_color_png) {
                flow_fwd_file = cpm_pf_params.output_PF_dir + "/" + flow_fwd_name + ".png";
                WriteFlowAsImage(pf_spatial_flow_vec[i], flow_fwd_file.c_str(), -1);
            }
        }
        sPF_write_time.toc(" done in: ");
    }

    // temporal filter
    CTimer tPF_time;
    std::cout << "Running temporal permeability filter... " << flush;
        
    vector<Mat2f> pf_temporal_flow_vec(nb_imgs);
    pf_temporal_flow_vec[0] = pf_spatial_flow_vec[0];

    // PF.init_T<Vec2f>(height, width); // Initializes PF internal accumulated buffer before loop
    // Init PF internal accumulated buffer before loop here, a bit ugly but best solution to handle multiple types for J
    Mat2f l_flow_t0_num = Mat2f::zeros(height, width);
    Mat2f l_flow_t0_den = Mat2f::zeros(height, width);
    for (size_t i = 1; i < nb_imgs; ++i)
    {
        PF.set_I_T(input_RGB_images_vec[i-1], input_RGB_images_vec[i]);
        

        Mat2f flow_t0_XY = pf_temporal_flow_vec[i - 1];
        Mat2f flow_t1_XY = pf_spatial_flow_vec[i];
        PF.set_flow_T(flow_t0_XY, flow_t1_XY);

        PF.computeTemporalPermeability();
        pf_temporal_flow_vec[i] = PF.filterT<Vec2f>(flow_t0_XY, flow_t1_XY, l_flow_t0_num, l_flow_t0_den);
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

            // Forward matching only
            string flow_fwd_name = "PF_temporal__" + img_name1 + "__TO__" + img_name2;
            string flow_fwd_file = cpm_pf_params.output_PF_dir + "/" + flow_fwd_name + ".flo";
            WriteFlowFile(pf_temporal_flow_vec[i], flow_fwd_file.c_str());
            
            if(cpm_pf_params.write_color_png) {
                flow_fwd_file = cpm_pf_params.output_PF_dir + "/" + flow_fwd_name + ".png";
                WriteFlowAsImage(pf_temporal_flow_vec[i], flow_fwd_file.c_str(), -1);
            }
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

        image_delete(flow_x);
        image_delete(flow_y);
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

        // Forward matching only
        string flow_fwd_name = "VR__" + img_name1 + "__TO__" + img_name2;
        string flow_fwd_file = cpm_pf_params.output_VR_dir + "/" + flow_fwd_name + ".flo";
        WriteFlowFile(vr_flow_vec[i], flow_fwd_file.c_str());
        
        if(cpm_pf_params.write_color_png) {
            flow_fwd_file = cpm_pf_params.output_VR_dir + "/" + flow_fwd_name + ".png";
            WriteFlowAsImage(vr_flow_vec[i], flow_fwd_file.c_str(), -1);
        }
    }
    vr_write_time.toc(" done in: ");
    
    total_time.toc("\nTotal elapsed time: ");

    return 0;
}
