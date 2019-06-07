/*  Name:
 *      utils.cpp
 *
 *  Description:
 *      Miscellaneous utilities, mainly conversion between different image formats / data types
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
	ou.setValue(UNKNOWN_FLOW);
	ov.setValue(UNKNOWN_FLOW);
	int cnt = inMat.height();
	for (int i = 0; i < cnt; i++){
		float* p = inMat.rowPtr(i);
		float x = p[0];
		float y = p[1];
		float u = p[2] - p[0];
		float v = p[3] - p[1];

		for (int di = -1; di <= 1; di++){
            for (int dj = -1; dj <= 1; dj++){
				int tx = ImageProcessing::EnforceRange(x + dj, w);
				int ty = ImageProcessing::EnforceRange(y + di, h);
				ou[ty*w + tx] = u;
				ov[ty*w + tx] = v;
			}
		}

	}
}

void Match2Mat2f(FImage matches, Mat2f &flow)
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
				int tx = ImageProcessing::EnforceRange(x + dj, w);
				int ty = ImageProcessing::EnforceRange(y + di, h);
				flow(ty, tx)[0] = u;
				flow(ty, tx)[1] = v;
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

void FImage2image_t(FImage mu, image_t* flow_x) {
    if(!(mu.nelements() > 0))
        cout<<"Can't convert from empty image!"<<endl;
    memcpy(flow_x->data, mu.pData, sizeof(float) * mu.nelements());
}

void Mat2f2image_t_uv(Mat2f flow, image_t* flow_x, image_t* flow_y) {
for (size_t y = 0; y < flow.rows; ++y) {
        for (size_t x = 0; x < flow.cols; ++x) {
            flow_x->data[y * flow_x->stride + x] = flow(y, x)[0];
            flow_y->data[y * flow_y->stride + x] = flow(y, x)[1];
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

void Mat3f2FImage(Mat3f in, FImage &out) {
    for (size_t y = 0; y < in.rows; ++y) {
        for (size_t x = 0; x < in.cols; ++x) {
                float pixBGR[3] = {in(y, x)[0], in(y, x)[1], in(y, x)[2]};
                out.setPixel(y, x, pixBGR);
        }
    }
}

void Mat3f2color_image_t(Mat3f in, color_image_t *out) {
    for (size_t y = 0; y < in.rows; ++y) {
        for (size_t x = 0; x < in.cols; ++x) {
                out->c1[y * out->stride + x] = in(y, x)[2] * 255.; // Opencv stores pixels in BGR order
                out->c2[y * out->stride + x] = in(y, x)[1] * 255.;
                out->c3[y * out->stride + x] = in(y, x)[0] * 255.;
        }
    }
}

void image_t2Mat(image_t* in, Mat &out) {
    for (size_t y = 0; y < out.rows; ++y) {
        for (size_t x = 0; x < out.cols; ++x) {
            out.at<float>(y, x) = in->data[y * in->stride + x];
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

