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

#include "utils.h"

/* ---------------- CONVERSION BETWEEN IMAGE TYPES --------------------------- */

void Match2Flow(FImage& inMat, FImage& ou, FImage& ov, int w, int h)
{
	if (!ou.matchDimension(w, h, 1)){
		ou.allocate(w, h, 1);
	}
	if (!ov.matchDimension(w, h, 1)){
		ov.allocate(w, h, 1);
	}
	ou.setValue(kMOVEMENT_UNKNOWN);
	ov.setValue(kMOVEMENT_UNKNOWN);
	int cnt = inMat.height();
	for (int i = 0; i < cnt; i++){
		float* p = inMat.rowPtr(i);
		float x = p[0];
		float y = p[1];
		float u = p[2] - p[0];
		float v = p[3] - p[1];

		for (int di = -1; di <= 1; di++){
            for (int dj = -1; dj <= 1; dj++){
				int tx = img_boundary_check(x + (float)dj, w);
				int ty = img_boundary_check(y + (float)di, h);
				ou[ty*w + tx] = u;
				ov[ty*w + tx] = v;
			}
		}

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
	}
	fclose(fid);
}

/* ---------------- FImage <-> Mat --------------------------- */
void Mat3f2FImage(Mat3f in, FImage &out) {
    for (size_t y = 0; y < in.rows; ++y) {
        for (size_t x = 0; x < in.cols; ++x) {
                float pixBGR[3] = {in(y, x)[0], in(y, x)[1], in(y, x)[2]};
                out.setPixel(y, x, pixBGR);
        }
    }
}

void Match2Flow(FImage matches, Mat2f &flow)
{
	int h = flow.rows;
    int w = flow.cols;
	int cnt = matches.height();
	for (int i = 0; i < cnt; i++){
		float* p = matches.rowPtr(i);
		float x = p[0];
		float y = p[1];
		float u = p[2] - p[0];
		float v = p[3] - p[1];

		for (int di = -1; di <= 1; di++){
            for (int dj = -1; dj <= 1; dj++){
				int tx = img_boundary_check(x + (float)dj, w);
				int ty = img_boundary_check(y + (float)di, h);
				flow(ty, tx)[0] = u;
				flow(ty, tx)[1] = v;
			}
		}

	}
}

void Match2Disp(FImage matches, Mat1f &disp, string parallax)
{
	if(parallax != "ver" && parallax != "vertical" && parallax != "hor" && parallax != "horizontal") 
	{
		cerr << "Wrong parallax direction in conversion from matches to disparity, should be ver or vertical or hor or horizontal" << endl;
		exit(EXIT_FAILURE);
	}

	int h = disp.rows;
    int w = disp.cols;
	int cnt = matches.height();
	for (int i = 0; i < cnt; i++){
		float* p = matches.rowPtr(i);
		float x = p[0];
		float y = p[1];
		float d;
		
		if(parallax == "hor" || parallax == "horizontal")
			d = p[2] - p[0];
		else if(parallax == "ver" || parallax == "vertical")
			d = p[3] - p[1];

		for (int di = -1; di <= 1; di++){
            for (int dj = -1; dj <= 1; dj++){
				int tx = img_boundary_check(x + (float)dj, w);
				int ty = img_boundary_check(y + (float)di, h);
				disp(ty, tx) = d;
			}
		}
	}
}




/* ---------------- FImage <-> image_t --------------------------- */
void FImage2image_t(FImage mu, image_t* flow_x) {
    if(!(mu.nelements() > 0))
        cout<<"Can't convert from empty image!"<<endl;
    memcpy(flow_x->data, mu.pData, sizeof(float) * mu.nelements());
}


/* ---------------- image_t <-> Mat --------------------------- */
void Mat3f2color_image_t(Mat3f in, color_image_t *out) {
    for (size_t y = 0; y < in.rows; ++y) {
        for (size_t x = 0; x < in.cols; ++x) {
                out->c1[y * out->stride + x] = in(y, x)[2] * 255.; // Opencv stores pixels in BGR order
                out->c2[y * out->stride + x] = in(y, x)[1] * 255.;
                out->c3[y * out->stride + x] = in(y, x)[0] * 255.;
        }
    }
}

void Mat1f2image_t(Mat1f disp_in, image_t* disp_out) {
for (size_t y = 0; y < disp_in.rows; ++y) {
        for (size_t x = 0; x < disp_in.cols; ++x) {
            disp_out->data[y * disp_out->stride + x] = disp_in(y, x);
        }
    }
}

void Mat2f2image_t_uv(Mat2f flow, image_t* flow_x, image_t* flow_y) {
for (size_t y = 0; y < flow.rows; ++y) {
        for (size_t x = 0; x < flow.cols; ++x) {
            flow_x->data[y * flow_x->stride + x] = flow(y, x)[0];
            flow_y->data[y * flow_y->stride + x] = flow(y, x)[1];
        }
    }
}

void image_t2Mat1f(image_t* in, Mat1f &out) {
    for (size_t y = 0; y < out.rows; ++y) {
        for (size_t x = 0; x < out.cols; ++x) {
            out(y, x) = in->data[y * in->stride + x];
        }
    }
}

void image_t_uv2Mat2f(Mat2f &flow, image_t* flow_x, image_t* flow_y) {
for (size_t y = 0; y < flow.rows; ++y) {
        for (size_t x = 0; x < flow.cols; ++x) {
            flow(y, x)[0] = flow_x->data[y * flow_x->stride + x];
            flow(y, x)[1] = flow_y->data[y * flow_y->stride + x];
        }
    }
}

/* ---------------- OPERATIONS ON STRING FOR FILE NAMING --------------------------- */
bool str_replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

