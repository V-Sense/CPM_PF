#include "CPM/CPM.h"
#include "CPM/OpticFlowIO.h"
#include "PFilter/PermeabilityFilter.h"
#include "flowIO.h"
#include "utils.h"
extern "C" {
#include "Variational_refinement/variational.h"
}


void Usage()
{
    cout<< "Example use of CPM_PF" << endl
        << "C++ implementation." << endl
        << endl
        << "Usage:" << endl
        << "  ./CPMPF <input_image_folder> <CPM_match_folder> <CPMPF_flow_folder> <refined_CPMPF_flow_folder> [options]" << endl
        << "options:" << endl
        << "    -h help                                     print this message" << endl
        << "  CPM parameters:" << endl
        << "    -m, -max                                    outlier handling maxdisplacement threshold" << endl
        << "    -t, -th                                     froward and backward consistency threshold" << endl
        << "    -c, -cth                                    matching cost check threshold" <<endl
        << "  PF parameters:" << endl
        << "    -i, -iter                                   number of iterantions for spatial permeability filter" << endl
        << "    -l, -lambda                                 lambda para for spatial permeability filter" << endl
        << "    -d, -delta                                  delta para for spatial permeability filter" << endl
        << "    -a, -alpha                                  alpha para for spatial permeability filter" << endl
        << "  predefined parameters:" << endl
        << "    -sintel                                     set the parameters to the one optimized on (a subset of) the MPI-Sintel dataset" << endl
        << "    -hcilf                                      set the parameters to the one optimized on (a subset of) the HCI light field dataset" << endl
        << endl;
}

int main(int argc, char** argv)
{
    if (argc < 5){
        if (argc > 1) fprintf(stderr, "Error: not enough arguments\n");
        Usage();
        exit(1);
	}

    CTimer total_time;

    // load inputs
    char* input_images_folder = argv[1];
    char* CPM_matches_folder = argv[2];
    char* CPMPF_flows_folder = argv[3];
    char* refined_CPMPF_flow_folder = argv[4];
    
    // prepare variables
    cpm_pf_params_t params;
    cpm_pf_params_t &cpm_pf_params = params;

    // load options
    #define isarg(key)  !strcmp(a,key)
    int current_arg = 5;
    while(current_arg < argc){
        const char* a = argv[current_arg++];
        if( isarg("-h") || isarg("-help") )
            Usage();
        else if( isarg("-m") || isarg("-max") )
            cpm_pf_params.max_displacement_input_int = atoi(argv[current_arg++]);
        else if( isarg("-t") || isarg("-th") )
            cpm_pf_params.check_threshold_input_int = atoi(argv[current_arg++]);
        else if( isarg("-c") || isarg("-cth") )
            cpm_pf_params.cost_threshold_input_int = atoi(argv[current_arg++]);
        else if( isarg("-i") || isarg("-iter") )
            cpm_pf_params.iterations_input_int = atof(argv[current_arg++]);
        else if( isarg("-l") || isarg("-lambda") )
            cpm_pf_params.lambda_XY_input_float = atof(argv[current_arg++]);
        else if( isarg("-d") || isarg("-delta") )
            cpm_pf_params.delta_XY_input_float = atof(argv[current_arg++]);
        else if( isarg("-a") || isarg("-alpha") )
            cpm_pf_params.alpha_XY_input_float = atof(argv[current_arg++]);
        else if( isarg("-sintel") ) {
            cpm_pf_params.max_displacement_input_int = 400;
            cpm_pf_params.check_threshold_input_int = 1;
            cpm_pf_params.cost_threshold_input_int = 1880;
            cpm_pf_params.iterations_input_int = 5;
            cpm_pf_params.lambda_XY_input_float = 0;
            cpm_pf_params.delta_XY_input_float = 0.017;
            cpm_pf_params.alpha_XY_input_float = 2;
        }
        else if( isarg("-hcilf") ) {
            cpm_pf_params.max_displacement_input_int = 4;
            cpm_pf_params.check_threshold_input_int = 1;
            cpm_pf_params.cost_threshold_input_int = 1880;
            cpm_pf_params.iterations_input_int = 5;
            cpm_pf_params.lambda_XY_input_float = 0;
            cpm_pf_params.delta_XY_input_float = 0.017;
            cpm_pf_params.alpha_XY_input_float = 2;
        }
        else {
            fprintf(stderr, "unknown argument %s", a);
            Usage();
            exit(1);
        }
    }

    // prepare inputs/outputs folders
    string input_images_folder_string = input_images_folder;
    string CPM_matches_folder_string = CPM_matches_folder;
    string CPMPF_flows_folder_string = CPMPF_flows_folder;
    string refined_CPMPF_flow_folder_string = refined_CPMPF_flow_folder;

    // input RGB images names
    string img_pre("SAI_08_"), img_suf(".png"), img_ext(".png");
    size_t start_idx = 7, nb_imgs = 5;

    vector<string> input_images_name_vec(nb_imgs);
    
    
    /* ---------------- READ INPUT RBG IMAGES --------------------------- */
    CTimer CPM_input_time;
    cout << "Reading input RGB images... " << flush;
    vector<Mat3f> input_RGB_images_vec(nb_imgs);
    int width = -1 , height = -1, nch = -1;

    for (size_t i = 0; i < nb_imgs; i++) {
        std::stringstream ss_idx;
        ss_idx << std::setw(2) << std::setfill('0') << start_idx + i;
		std::string s_idx = ss_idx.str();
        String img_name = img_pre + s_idx + img_suf;
        input_images_name_vec[i] = img_name;
        
        Mat tmp_img = imread(input_images_folder_string + "/" + img_name);
        if ( tmp_img.empty() ) {
            cout << input_images_name_vec[i] << " is invalid!" << endl;
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
                cout << "All input images should have the same size; size of image " << input_images_name_vec[i] << " is different from previous image" << endl;
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
    cout << "Running CPM... " << flush;
    int step = 3;
    CPM cpm(cpm_pf_params);
    cpm.SetStep(step);
    vector<Mat2f> cpm_flow_fwd(nb_imgs-1), cpm_flow_bwd(nb_imgs-1);
    
    // #pragma omp parallel for 
    for (size_t i = 0; i < nb_imgs - 1; ++i) {
        FImage img1(width, height, nch);
        FImage img2(width, height, nch);
        
        Mat3f2FImage(input_RGB_images_vec[i], img1);
        Mat3f2FImage(input_RGB_images_vec[i+1], img2);

        FImage matches;
        

        // Forward flow
        cpm.Matching(img1, img2, matches);

        Mat2f flow_fwd(height, width, UNKNOWN_FLOW);
        // flow_fwd = UNKNOWN_FLOW;
        Match2Mat2f(matches, flow_fwd);
        cpm_flow_fwd[i] = flow_fwd;

        // Backward flow
        matches.clear();
        cpm.Matching(img2, img1, matches);

        Mat2f flow_bwd(height, width, UNKNOWN_FLOW);
        Match2Mat2f(matches, flow_bwd);
        cpm_flow_bwd[i] = flow_bwd;
    }
    CPM_time.toc(" done in: ");

    // Write results on disk
    ostringstream cpm_matches_name_builder, cpm_matches_name_flo_builder, cpm_matches_name_png_builder, cpm_matches_name_txt_builder;

    CTimer CPM_write_time;
    cout << "Writing CPM results... " << flush;
    #pragma omp parallel for 
    for (size_t i = 0; i < nb_imgs - 1; ++i) {
        // Remove image file extension
        string img_name1 = input_images_name_vec[i];
        string img_name2 = input_images_name_vec[i+1];
        str_replace(img_name1, img_ext, "");
        str_replace(img_name2, img_ext, "");

        // Forward matching
        string flow_fwd_name = "CPM__" + img_name1 + "__TO__" + img_name2;
        string flow_fwd_file = CPM_matches_folder_string + "/" + flow_fwd_name + ".flo";
        WriteFlowFile(cpm_flow_fwd[i], flow_fwd_file.c_str());
        
        flow_fwd_file = CPM_matches_folder_string + "/" + flow_fwd_name + ".png";
        WriteFlowAsImage(cpm_flow_fwd[i], flow_fwd_file.c_str(), -1);


        // Backward matching
        string flow_bwd_name = "CPM__" + img_name2 + "__TO__" + img_name1;
        string flow_bwd_file = CPM_matches_folder_string + "/" + flow_bwd_name + ".flo";
        WriteFlowFile(cpm_flow_bwd[i], flow_bwd_file.c_str());
        
        flow_bwd_file = CPM_matches_folder_string + "/" + flow_bwd_name + ".png";
        WriteFlowAsImage(cpm_flow_bwd[i], flow_bwd_file.c_str(), -1);
    }
    CPM_write_time.toc(" done in: ");

    /* ---------------- RUN PERMEABILITY FILTER --------------------------- */
    // spatial filter
    CTimer sPF_time;
    cout << "Running spatial permeability filter... " << flush;
    vector<Mat2f> pf_spatial_flow_vec;

    for (size_t i = 0; i < cpm_flow_fwd.size(); ++i) {
        Mat3f target_img = input_RGB_images_vec[i];
        Mat2f flow_forward  = cpm_flow_fwd[i];
        Mat2f flow_backward = cpm_flow_bwd[i];

        // compute flow confidence map
        Mat1f flow_confidence = getFlowConfidence(flow_forward, flow_backward);
        // Mat flow_confidence_mat;
        // flow_confidence.convertTo(flow_confidence_mat, CV_32FC1);

        // Apply spatial permeability filter on confidence
        Mat1f flow_confidence_filtered = filterXY<Vec3f>(target_img, flow_confidence, cpm_pf_params);
        
        // // Mat flow_confidence_filtered_mat;
        // flow_confidence_filtered.convertTo(flow_confidence_filtered_mat, CV_32FC1);
        // WriteFilePFM(flow_confidence_filtered_mat, format("00%d_flow_confidence_filtered_XY_mat_1_channels.pfm", i), 1/255.);

        // multiply initial confidence and sparse flow
        Mat2f confidenced_flow = Mat2f::zeros(flow_confidence.rows,flow_confidence.cols);
        for(int y = 0; y < confidenced_flow.rows; y++) {
            for(int x = 0; x < confidenced_flow.cols; x++) {
                for(int c = 0; c < confidenced_flow.channels(); c++) {
                    confidenced_flow(y,x)[c] = flow_forward(y,x)[c] * flow_confidence(y,x);
                }
            }
        }
        //string temp_str1 = format("00%d_confidenced_flow.flo", i);
        //WriteFlowFile(confidenced_flow, temp_str1.c_str());

        //filter confidenced sparse flow
        Mat2f confidenced_flow_XY = filterXY<Vec3f, Vec2f>(target_img, confidenced_flow, cpm_pf_params);
        //string temp_str2 = format("00%d_confidenced_flow_XY.flo", i);
        //WriteFlowFile(confidenced_flow_XY, temp_str2.c_str());

        // compute normalized spatial filtered flow FXY by division
        Mat2f normalized_confidenced_flow_filtered = Mat2f::zeros(target_img.rows,target_img.cols);
        for(int y = 0; y < confidenced_flow_XY.rows; y++) {
            for(int x = 0; x < confidenced_flow_XY.cols; x++) {
                for(int c = 0; c < confidenced_flow_XY.channels(); c++) {
                    normalized_confidenced_flow_filtered(y,x)[c] = confidenced_flow_XY(y,x)[c] / flow_confidence_filtered(y,x);
                }
            }
        }

        pf_spatial_flow_vec.push_back(normalized_confidenced_flow_filtered);
        
        ostringstream temp_str3_builder;
        temp_str3_builder << CPMPF_flows_folder_string << setw(4) << setfill('0') << i + 1 << "_Normalized_Flow_XY.flo";
        string temp_str3 = temp_str3_builder.str();
        WriteFlowFile(normalized_confidenced_flow_filtered, temp_str3.c_str());

        vector<Mat1f> normalized_confidenced_flow_filtered_vec;
        split(normalized_confidenced_flow_filtered, normalized_confidenced_flow_filtered_vec);
        
        Mat normalized_confidenced_flow_filtered_mat;
        normalized_confidenced_flow_filtered_vec[0].convertTo(normalized_confidenced_flow_filtered_mat, CV_32FC1);
        temp_str3_builder.str("");
        temp_str3_builder.clear();
        temp_str3_builder << CPMPF_flows_folder_string << setw(4) << setfill('0') << i + 1 << "_Normalized_Flow_XY.pfm";
        temp_str3 = temp_str3_builder.str();
        WriteFilePFM(normalized_confidenced_flow_filtered_mat, temp_str3.c_str(), 1/255.);

    }
    sPF_time.toc(" done in: ");

    // temporal filter
    CTimer tPF_time;
    cout << "Running temporal permeability filter... " << flush;
    Mat2f l_prev = Mat2f::zeros(height, width);
    Mat2f l_normal_prev = Mat2f::zeros(height, width);
    Mat2f It0_XYT, It1_XYT;
    vector<Mat2f> It1_XYT_vector;
    vector<Mat2f> pf_temporal_flow_vec;

    for (size_t i = 1; i < pf_spatial_flow_vec.size(); ++i)
    {
        Mat3f It0 = input_RGB_images_vec[i - 1];
        Mat3f It1 = input_RGB_images_vec[i];
        Mat2f It0_XY = pf_spatial_flow_vec[i - 1];
        Mat2f It1_XY = pf_spatial_flow_vec[i];

        if(i == 1) {
            It1_XYT_vector = filterT<Vec3f, Vec2f>(It1, It0, It1_XY, It0_XY, It1_XY, It0_XY,  l_prev, l_normal_prev);
        }
        else {
            It1_XYT_vector = filterT<Vec3f, Vec2f>(It1, It0, It1_XY, It0_XY, It1_XY, It0_XYT, l_prev, l_normal_prev);
        }

        It1_XYT = It1_XYT_vector[2];
        l_prev = It1_XYT_vector[0];
        l_normal_prev = It1_XYT_vector[1];

        pf_temporal_flow_vec.push_back(It1_XYT);

        ostringstream flowXYT1_name_builder;
        flowXYT1_name_builder << CPMPF_flows_folder_string << setw(4) << setfill('0') << i + 1 << "_XYT.flo";
        string flowXYT1_name = flowXYT1_name_builder.str();
        WriteFlowFile(It1_XYT, flowXYT1_name.c_str());
        It0_XYT = It1_XYT;
    }
    tPF_time.toc(" done in: ");


    /* ---------------- RUN VARIATIONAL REFINEMENT  --------------------------- */
    CTimer var_time;
    cout << "Running variational refinement... " << flush;
    color_image_t *im1 = color_image_new(width, height);
    color_image_t *im2 = color_image_new(width, height);
    for (size_t i = 0; i < pf_temporal_flow_vec.size(); ++i) {
        
        ostringstream refined_cpmpf_flows_name_builder;
        if (i == 0) {
            Mat3f2color_image_t(input_RGB_images_vec[i], im1);
            Mat3f2color_image_t(input_RGB_images_vec[i + 1], im2);
            refined_cpmpf_flows_name_builder << setw(4) << setfill('0') << i + 1 << "_refined_Normalized_Flow_XY";
        }
        else {
            if (i % 2 != 0) {
                Mat3f2color_image_t(input_RGB_images_vec[(i + 1) / 2], im1);
                Mat3f2color_image_t(input_RGB_images_vec[(i + 1) / 2 + 1], im2);
                refined_cpmpf_flows_name_builder << setw(4) << setfill('0') << (i + 1) / 2 + 1 << "_refined_Normalized_Flow_XY";
            }
            else {
                Mat3f2color_image_t(input_RGB_images_vec[i / 2], im1);
                Mat3f2color_image_t(input_RGB_images_vec[i / 2 + 1], im2);
                refined_cpmpf_flows_name_builder << setw(4) << setfill('0') << i / 2 + 1 << "_refined_XYT";
            }
        }

        Mat2f pf_flow = pf_temporal_flow_vec[i];

        variational_params_t flow_params;
        variational_params_default(&flow_params);
        image_t *flow_x = image_new(width, height), *flow_y = image_new(width, height);

        Mat2f2image_t_uv(pf_flow, flow_x, flow_y);

        variational(flow_x, flow_y, im1, im2, &flow_params);


        ostringstream refined_cpmpf_flows_name_flo_builder, refined_cpmpf_flows_name_png_builder;
        refined_cpmpf_flows_name_flo_builder << refined_CPMPF_flow_folder_string << refined_cpmpf_flows_name_builder.str() << ".flo";
        refined_cpmpf_flows_name_png_builder << refined_CPMPF_flow_folder_string << refined_cpmpf_flows_name_builder.str() << ".png";

        string refined_cpm_matches_name_flo = refined_cpmpf_flows_name_flo_builder.str();
        string refined_cpm_matches_name_png = refined_cpmpf_flows_name_png_builder.str();

        vector<Mat1f> flow_vec;
        flow_vec.push_back(Mat1f::zeros(height, width));
        flow_vec.push_back(Mat1f::zeros(height, width));
        image_t2Mat(flow_x, flow_vec[0]);
        image_t2Mat(flow_y, flow_vec[1]);

        Mat2f flow_mat;
        merge(flow_vec, flow_mat);
        WriteFlowFile(flow_mat, refined_cpm_matches_name_flo.c_str());

        // Convert final result to pfm
        // CTimer convert2disp_time;
        Mat disp(height, width, CV_32F);
        image_t2Mat(flow_x, disp);
        disp = -disp;
        // convert2disp_time.toc("Conversion of image_t to Mat: ");

        ostringstream refined_cpmpf_flows_name_pfm_builder;
        refined_cpmpf_flows_name_pfm_builder << refined_CPMPF_flow_folder_string << refined_cpmpf_flows_name_builder.str() << ".pfm";
        string refined_cpm_matches_name_pfm = refined_cpmpf_flows_name_pfm_builder.str();
        WriteFilePFM(disp, refined_cpmpf_flows_name_pfm_builder.str(), 1/255.0);
    }
    color_image_delete(im1);
    color_image_delete(im2);
    var_time.toc(" done in: ");

    
    total_time.toc("\nTotal elapsed time: ");

    return 0;
}
