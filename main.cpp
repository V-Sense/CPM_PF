
#include "CPM_Tip2017Mod/CPM.h"
#include "CPM_Tip2017Mod/OpticFlowIO.h"
#include "PFilter/PermeabilityFilter.h"
#include "flowIO.h"
extern "C" {
#include "PFilter/variational/variational.h"
#include "PFilter/variational/io.h"
}


// draw each match as a 3x3 color block
void Match2Flow(FImage& inMat, FImage& ou, FImage& ov, int w, int h)
{
	if (!ou.matchDimension(w, h, 1)){
		ou.allocate(w, h, 1);
	}
	if (!ov.matchDimension(w, h, 1)){
		ov.allocate(w, h, 1);
	}
	ou.setValue(UNKNOWN_FLOW);
	ov.setValue(UNKNOWN_FLOW);
	int cnt = inMat.height();
	for (int i = 0; i < cnt; i++){
		float* p = inMat.rowPtr(i);
		float x = p[0];
		float y = p[1];
		float u = p[2] - p[0];
		float v = p[3] - p[1];
/*
        if (i < 10) {
            printf("x is %f\n", x);
            printf("y is %f\n", y);
            printf("u is %f\n", u);
            printf("v is %f\n\n", v);
        }
*/

		for (int di = -1; di <= 1; di++){
            for (int dj = -1; dj <= 1; dj++){
				int tx = ImageProcessing::EnforceRange(x + dj, w);
				int ty = ImageProcessing::EnforceRange(y + di, h);
				ou[ty*w + tx] = u;
				ov[ty*w + tx] = v;
			}
		}

        /*
        for (int di = 0; di <= 0; di++){
            for (int dj = 0; dj <= 0; dj++){
                int tx = ImageProcessing::EnforceRange(x + dj, w);
                int ty = ImageProcessing::EnforceRange(y + di, h);
                ou[ty*w + tx] = u;
                ov[ty*w + tx] = v;
            }
        }
        */
	}
}

void WriteMatches(const char *filename, FImage& inMat)
{
	int len = inMat.height();
	FILE *fid = fopen(filename, "w");
	for (int i = 0; i < len; i++){
		float x1 = inMat[4 * i + 0];
		float y1 = inMat[4 * i + 1];
		float x2 = inMat[4 * i + 2];
		float y2 = inMat[4 * i + 3];
		fprintf(fid, "%.0f %.0f %.0f %.0f\n", x1, y1, x2, y2);
		//fprintf(fid, "%.3f %.3f %.3f %.3f 1 100\n", x1, y1, x2, y2);
	}
	fclose(fid);
}

void FImage2image_t(FImage mu, image_t* wx) {
    //image_t *tmp = image_new(mu->width(), mu->height());
    //FImage tmp2(*mu);
    if(!(mu.nelements() > 0))
        cout<<"Can't convert from empty image!"<<endl;
    memcpy(wx->data, mu.pData, sizeof(float) * mu.nelements());
    //tmp->data = tmp2.pData;
    //wx = image_cpy(tmp);
}

void Mat2f2image_t_uv(Mat2f flow, image_t* wx, image_t* wy) {
    vector<Mat1f> flow_ch1_vec;
    split(flow, flow_ch1_vec);

    Mat1f u_flow_ch1 = flow_ch1_vec[0];
    Mat1f v_flow_ch1 = flow_ch1_vec[1];


    //memcpy(wx->data, u_flow_ch1.data, u_flow_ch1.total());
    //memcpy(wy->data, v_flow_ch1.data, v_flow_ch1.total());


    for (size_t y = 0; y < u_flow_ch1.rows; ++y) {
        for (size_t x = 0; x < u_flow_ch1.cols; ++x) {
    //for (size_t y = 0; y < 4; ++y) {
        //for (size_t x = 0; x < 4; ++x) {
            wx->data[y * wx->stride + x] = u_flow_ch1(y, x);
            wy->data[y * wy->stride + x] = v_flow_ch1(y, x);
            //cout<<"y is "<< y <<endl;
            //cout<<"x is "<< x <<endl;
            //cout<<"stride is "<<wx->stride<<endl;
            //cout<<"y * wx->stride + x is "<<y * wx->stride + x<<endl;
            //cout<<"wx->data[y * wx->stride + x] is "<< wx->data[y * wx->stride + x] <<endl;
            //cout<<"u_flow_ch1(y,x) is "<< u_flow_ch1(y, x) <<endl;
        }
    }
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

void run_CPM(FImage img1, FImage img2, int seq_num_of_img1, bool is_forward_matching, cpm_pf_params_t &cpm_pf_params, string output_matches_folder)
{
    int step = 3;
    int w = img1.width();
    int h = img1.height();

    CTimer totalT;
    FImage matches;

    CPM cpm(cpm_pf_params);
    cpm.SetStep(step);
    cpm.Matching(img1, img2, matches);

    totalT.toc("CPM total time: ");

    //FImage u,v;
    //char tmpName[256];

    ostringstream cpm_matches_name_builder;
    if ( is_forward_matching ) {
        cpm_matches_name_builder << setw(4) << setfill('0') << seq_num_of_img1 << '_' << setw(4) << setfill('0') << seq_num_of_img1 + 1;
    }
    else {
        cpm_matches_name_builder << setw(4) << setfill('0') << seq_num_of_img1 << '_' << setw(4) << setfill('0') << seq_num_of_img1 - 1;
    }

    ostringstream cpm_matches_name_flo_builder, cpm_matches_name_png_builder, cpm_matches_name_txt_builder;
    cpm_matches_name_flo_builder << output_matches_folder << cpm_matches_name_builder.str() << ".flo";
    cpm_matches_name_png_builder << output_matches_folder << cpm_matches_name_builder.str() << ".png";
    cpm_matches_name_txt_builder << output_matches_folder << cpm_matches_name_builder.str() << ".txt";
    string cpm_matches_name_flo = cpm_matches_name_flo_builder.str();
    string cpm_matches_name_png = cpm_matches_name_png_builder.str();
    string cpm_matches_name_txt = cpm_matches_name_txt_builder.str();

    FImage u, v;
    Match2Flow(matches, u, v, w, h);
    OpticFlowIO::WriteFlowFile(u.pData, v.pData, w, h, cpm_matches_name_flo.c_str());
    OpticFlowIO::SaveFlowAsImage(cpm_matches_name_png.c_str(), u.pData, v.pData, w, h);
    WriteMatches(cpm_matches_name_txt.c_str(), matches);


    //vector<FImage> output_vec;
    //output_vec.push_back(u);
    //output_vec.push_back(v);

    //return output_vec;
}

//void run_PF(Mat3f target_img, Mat2f flow_forward, Mat2f flow_backward)
//{


//}


int main(int argc, char** argv)
{
    if (argc < 5){
        if (argc > 1) fprintf(stderr, "Error: not enough arguments\n");
        Usage();
        exit(1);
	}

    // load inputs
    char* input_images_folder = argv[1];
    char* CPM_matches_folder = argv[2];
    //char* refined_CPM_matches_folder = argv[3];
    char* CPMPF_flows_folder = argv[3];
    char* refined_CPMPF_flow_folder = argv[4];
    //int max_displacement_input_int, check_threshold_input_int, cost_threshold_input_int, iterations_input_int;
    //float lambda_XY_input_float, delta_XY_input_float, alpha_XY_input_float;
//  float delta_photo_input_float, delta_grad_input_float, alpha_photo_input_float, alpha_grad_input_float;

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

    // perpare inputs/outputs folders
    String input_images_folder_string = input_images_folder;
    String CPM_matches_folder_string = CPM_matches_folder;
    String CPMPF_flows_folder_string = CPMPF_flows_folder;
    String refined_CPMPF_flow_folder_string = refined_CPMPF_flow_folder;

    vector<String> input_images_name_vec;
    glob(input_images_folder_string, input_images_name_vec);

    // preapre inputs/outputs for CPM part and var part
    vector<FImage> cpm_input_images_vec;
    for (size_t i = 0; i < input_images_name_vec.size(); i++) {
        FImage tmp_img;
        tmp_img.imread(input_images_name_vec[i].c_str());
        if ( tmp_img.IsEmpty() ) {
            cout << input_images_name_vec[i] << " is invalid!" << endl;
            continue;
        }
        cpm_input_images_vec.push_back( tmp_img );
        tmp_img.imwrite(input_images_name_vec[i].c_str());
    }


    // run CPM part and var part
    for (size_t i = 0; i < cpm_input_images_vec.size() - 1; ++i) {
        FImage img1, img2;
        img1.copy(cpm_input_images_vec[i]);
        img2.copy(cpm_input_images_vec[i+1]);

        int w = img1.width();
        int h = img1.height();
        if (img2.width() != w || img2.height() != h) {
            printf("CPM can only handle images with the same dimension!\n");
            return -1;
        }

        //vector<FImage> img1_uv_vec, img2_uv_vec;
        //img1_uv_vec = run_CPM(img1, img2, i + 1, true, cpm_pf_params, CPM_matches_folder_string);
        //img2_uv_vec = run_CPM(img2, img1, i + 2, false, cpm_pf_params, CPM_matches_folder_string);

        run_CPM(img1, img2, i + 1, true, cpm_pf_params, CPM_matches_folder_string);
        run_CPM(img2, img1, i + 2, false, cpm_pf_params, CPM_matches_folder_string);
    }


    // perpare inputs/outputs for PF part
    vector<Mat3f> pf_input_images_vec;
    vector<Mat2f> pf_input_matches_vec; // same as cpm_output_matches

    vector<String> pf_input_matches_name_vec;
    glob(CPM_matches_folder_string, pf_input_matches_name_vec);

    for (size_t i = 0; i < input_images_name_vec.size(); i++) {
        Mat tmp_img = imread(input_images_name_vec[i]);
        if ( tmp_img.empty() ) {
            cout << input_images_name_vec[i] << " is invalid!" << endl;
            continue;
        }

        Mat3f tmp_img3f;
        tmp_img.convertTo(tmp_img3f, CV_32F, 1/255.);
        pf_input_images_vec.push_back( tmp_img3f.clone() );
    }

    for (size_t i = 0; i < pf_input_matches_name_vec.size(); i++) {
        Mat2f tmp_flo;
        string tmp_flo_name = pf_input_matches_name_vec[i];
        const char* tmp_flo_name_char = tmp_flo_name.c_str();
        ReadFlowFile(tmp_flo, tmp_flo_name_char);
        if( tmp_flo.empty() ) {
            cout<< pf_input_matches_name_vec[i] << " is invalid!" << endl;
            continue;
        }
        pf_input_matches_vec.push_back( tmp_flo.clone() );
    }


    // run PF part
    // spatial filter
    for (size_t i = 0; i * 2 < pf_input_matches_vec.size(); ++i) {
        Mat3f target_img = pf_input_images_vec[i];
        Mat2f flow_forward = pf_input_matches_vec[i * 2];
        Mat2f flow_backward = pf_input_matches_vec[i * 2 + 1];

        // compute flow confidence map
        Mat1f flow_confidence = getFlowConfidence(flow_forward, flow_backward);
        Mat flow_confidence_mat;
        flow_confidence.convertTo(flow_confidence_mat, CV_32FC1);

        // start of filtering flow confidence map by copy 1 channel to 2 channel
        vector<Mat1f> flow_confidence_2chs_vec;
        flow_confidence_2chs_vec.push_back(flow_confidence);
        flow_confidence_2chs_vec.push_back(flow_confidence);
        Mat2f flow_confidence_2chs;
        merge(flow_confidence_2chs_vec, flow_confidence_2chs);

        Mat2f flow_confidence_2chs_filtered = filterXY<Vec3f, Vec2f>(target_img, flow_confidence_2chs, cpm_pf_params);

        vector<Mat1f> flow_confidence_2chs_filtered_vec;
        split(flow_confidence_2chs_filtered, flow_confidence_2chs_filtered_vec);
        Mat1f flow_confidence_filtered = flow_confidence_2chs_filtered_vec[0];
        Mat flow_confidence_filtered_mat;
        flow_confidence_filtered.convertTo(flow_confidence_filtered_mat, CV_32FC1);
        //WriteFilePFM(flow_confidence_filtered_mat, format("00%d_flow_confidence_filtered_XY_mat.pfm", i), 1/255.);
        //end of filtering flow confidence map by copy 1 channel to 2 channel


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
        //string temp_str3 = format("00%d_Normalized_Flow_XY.flo", i);
        ostringstream temp_str3_builder;
        temp_str3_builder << CPMPF_flows_folder_string << setw(4) << setfill('0') << i + 1 << "_Normalized_Flow_XY.flo";
        string temp_str3 = temp_str3_builder.str();
        WriteFlowFile(normalized_confidenced_flow_filtered, temp_str3.c_str());
    }

    // temporal filter
    Mat2f l_prev = Mat2f::zeros(pf_input_images_vec[0].rows,pf_input_images_vec[0].cols);
    Mat2f l_normal_prev = Mat2f::zeros(pf_input_images_vec[0].rows,pf_input_images_vec[0].cols);
    Mat2f It0_XYT, It1_XYT;
    vector<Mat2f> It1_XYT_vector;
    for (size_t i = 1; i * 2 < pf_input_matches_vec.size(); ++i)
    {
        Mat3f It0 = pf_input_images_vec[i - 1];
        Mat3f It1 = pf_input_images_vec[i];
        Mat2f It0_XY, It1_XY;

        //string flowXY0_name = format("00%d_Normalized_Flow_XY.flo", i - 1);
        ostringstream flowXY0_name_builder;
        flowXY0_name_builder << CPMPF_flows_folder_string << setw(4) << setfill('0') << i << "_Normalized_Flow_XY.flo";
        string flowXY0_name = flowXY0_name_builder.str();
        ReadFlowFile(It0_XY, flowXY0_name.c_str());

        //string flowXY1_name = format("00%d_Normalized_Flow_XY.flo", i);
        ostringstream flowXY1_name_builder;
        flowXY1_name_builder << CPMPF_flows_folder_string << setw(4) << setfill('0') << i + 1 << "_Normalized_Flow_XY.flo";
        string flowXY1_name = flowXY1_name_builder.str();
        ReadFlowFile(It1_XY, flowXY1_name.c_str());

        if(i == 1) {
            //printf("i==1");
            It1_XYT_vector = filterT<Vec3f, Vec2f>(It1, It0, It1_XY, It0_XY, It1_XY, It0_XY, l_prev, l_normal_prev);
        }
        else {
            //printf("i!=1");
            It1_XYT_vector = filterT<Vec3f, Vec2f>(It1, It0, It1_XY, It0_XY, It1_XY, It0_XYT, l_prev, l_normal_prev);
        }

        It1_XYT = It1_XYT_vector[2];
        l_prev = It1_XYT_vector[0];
        l_normal_prev = It1_XYT_vector[1];

        //string flowXYT1_name = format("00%d_XYT.flo", i);
        ostringstream flowXYT1_name_builder;
        flowXYT1_name_builder << CPMPF_flows_folder_string << setw(4) << setfill('0') << i + 1 << "_XYT.flo";
        string flowXYT1_name = flowXYT1_name_builder.str();
        WriteFlowFile(It1_XYT, flowXYT1_name.c_str());
        It0_XYT = It1_XYT;
    }

    //preapre inputs for var part
    vector<color_image_t*> var_input_images_vec;
    vector<Mat2f> var_input_flows_vec;

    vector<String> var_input_flows_name_vec;
    glob(CPMPF_flows_folder_string, var_input_flows_name_vec);

    for (size_t i = 0; i < input_images_name_vec.size(); i++) {
        color_image_t *tmp_img;
        tmp_img = color_image_load(input_images_name_vec[i].c_str());
        if ( tmp_img->stride == 0 ) {
            cout << input_images_name_vec[i] << " is invalid!" << endl;
            continue;
        }
        var_input_images_vec.push_back( tmp_img );
    }


    for (size_t i = 0; i < var_input_flows_name_vec.size(); i++) {
        Mat2f tmp_flo;
        string tmp_flo_name = var_input_flows_name_vec[i];
        const char* tmp_flo_name_char = tmp_flo_name.c_str();
        ReadFlowFile(tmp_flo, tmp_flo_name_char);
        if( tmp_flo.empty() ) {
            cout<< var_input_flows_name_vec[i] << " is invalid!" << endl;
            continue;
        }
        var_input_flows_vec.push_back( tmp_flo.clone() );
    }

    //run var part
    for (size_t i = 0; i < var_input_flows_vec.size(); ++i) {
        color_image_t *im1, *im2;
        Mat2f flo;
        ostringstream refined_cpmpf_flows_name_builder;

        if (i == 0) {
            im1 = color_image_cpy(var_input_images_vec[i]);
            im2 = color_image_cpy(var_input_images_vec[i + 1]);
            refined_cpmpf_flows_name_builder << setw(4) << setfill('0') << i + 1 << "_Normalized_Flow_XY";
        }
        else {
            if (i % 2 != 0) {
                im1 = color_image_cpy(var_input_images_vec[(i + 1) / 2]);
                im2 = color_image_cpy(var_input_images_vec[(i + 1) / 2 + 1]);
                refined_cpmpf_flows_name_builder << setw(4) << setfill('0') << (i + 1) / 2 + 1 << "_Normalized_Flow_XY";
            }
            else {
                im1 = color_image_cpy(var_input_images_vec[i / 2]);
                im2 = color_image_cpy(var_input_images_vec[i / 2 + 1]);
                refined_cpmpf_flows_name_builder << setw(4) << setfill('0') << i / 2 + 1 << "_XYT";
            }
        }

        flo = var_input_flows_vec[i];

        variational_params_t flow_params;
        variational_params_default(&flow_params);
        image_t *wx = image_new(im1->width, im1->height), *wy = image_new(im1->width, im1->height);
        //FImage2image_t(img1_uv_vec[0], wx);
        //FImage2image_t(img1_uv_vec[1], wy);

        Mat2f2image_t_uv(flo, wx, wy);

        variational(wx, wy, im1, im2, &flow_params);


        ostringstream refined_cpmpf_flows_name_flo_builder, refined_cpmpf_flows_name_png_builder;
        refined_cpmpf_flows_name_flo_builder << refined_CPMPF_flow_folder_string << refined_cpmpf_flows_name_builder.str() << ".flo";
        refined_cpmpf_flows_name_png_builder << refined_CPMPF_flow_folder_string << refined_cpmpf_flows_name_builder.str() << ".png";

        string refined_cpm_matches_name_flo = refined_cpmpf_flows_name_flo_builder.str();
        string refined_cpm_matches_name_png = refined_cpmpf_flows_name_png_builder.str();


        writeFlowFile(refined_cpm_matches_name_flo.c_str(), wx, wy);

    }

    printf("Hello World!");
    return 0;
}
