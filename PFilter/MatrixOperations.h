// Copyright Joan Sol Roo and Christian Richardt (7 August 2014)
//
// solroo@gmail.com
// christian@richardt.name
//
// This software is a computer program whose purpose is to compute temporally
// coherent stereoscopic 3D videos from anaglyph 3D input videos.
//
// This software is governed by the CeCILL license under French law and
// abiding by the rules of distribution of free software. You can use, modify
// and/or redistribute the software under the terms of the CeCILL license as
// circulated by CEA, CNRS and INRIA at "http://www.cecill.info".
//
// As a counterpart to the access to the source code and rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty and the software's author, the holder of the
// economic rights, and the successive licensors have only limited liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading, using, modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean that it is complicated to manipulate, and that also
// therefore means that it is reserved for developers and experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and, more generally, to use and operate it in the same
// conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL license and that you accept its terms.

#pragma once

#include <opencv2/opencv.hpp>


extern const int kVERTICAL;
extern const int kHORIZONTAL;
extern const int kFindNotFound;


// Overloaded function: adds to each channel of the matrix m3 the matching value from a single channel matrix m1
inline cv::Mat3f add(cv::Mat3f m3, cv::Mat1f m1);

// Overloaded funtion: computes the exponential function
//     R(X,Y) = a ^ I(X,Y)
inline cv::Mat pow(float a, cv::Mat src);

// Computes the distance between neigbours in the indicated direction, as an aproximation to the directional derivative
//     direction == kHORIZONTAL : R(X,Y) = I(X,Y) - I(X+1,Y)
//     direction == kVERTICAL :   R(X,Y) = I(X,Y) - I(X,Y+1)
// NOTE: The advanced filter option is based on the one provided at the Transform Domain implementation in Matlab, but was never used
template <class TMat>
inline TMat diff(TMat src, int direction = kHORIZONTAL, bool advanced_filter = false);

// Computes the cummulative sum in indicated direction, as an aproximation of the integral
//     direction == kHORIZONTAL : R(X,Y) = I(0,Y) + ... + I(X,Y)
//     direction == kVERTICAL :   R(X,Y) = I(X,0) + ... + I(X,Y)
template <class TField>
inline cv::Mat_<TField> cumsum(cv::Mat_<TField> src, int direction, int offset = 1);

// Just inline wrappers for OpenCV functions: Function(input,output) -> output = Function(input)
template <class TMat>
inline TMat copyMakeBorder(TMat src, int border_top, int border_bottom, int border_left, int border_right, int border_type = cv::BORDER_REPLICATE);
inline cv::Mat1b grayscale(const cv::Mat3f& image);
inline cv::Mat3f im2float(cv::InputArray src);
inline cv::Mat3b float2im(cv::InputArray src);

// Searchs for the first occurrence GREATER THAN the indicated value V in the given collection C , within the indicaded range
//     R = min(i) / C(i) > V
//     if no such value is found, the default value (kFindNotFound) is returned
// NOTE: this function is a simple implementation of the Matlab equivalent, required by the Domain transform
template <class TValue>
inline int findIndexOf(TValue value, TValue* collection, int range_min, int range_max, int offset = 0, int default_value = kFindNotFound);


// Applies a given mask M to an input image I, asigning 
//     M(X,Y) >  threshold : R(X,Y) = I(X,Y) * M(X,Y)
//     M(X,Y) <= threshold : R(X,Y) = default_value
// NOTE: This method is not optimized
template <class TValue>
cv::Mat_<TValue> applyMask(const cv::Mat_<TValue>&  image, const cv::Mat1f&  mask, float threshold = 0.0001f, TValue default_value = TValue());

// Divides each channel j of the image J by the corresponding value of the confidence image G
//     R(x,y)[j] = J(x,y)[j] / G(x,y)
template<class TValue>
inline cv::Mat_<TValue> normalize(cv::Mat_<TValue> J, cv::Mat1f G);

//////////////////////////////////////////////////////////////////////////
// template / inline functions must be in the same file of their headers
/////////////////////////////////////////////////////////////////////////

inline cv::Mat1b grayscale(const cv::Mat3f& image)
{
	cv::Mat1b resultb;
	cv::cvtColor(float2im(image), resultb, cv::COLOR_BGR2GRAY);
	return resultb;
}


inline cv::Mat3f add(cv::Mat3f m3, cv::Mat1f m1)
{
	int num_channels = m3.channels();
	std::vector<cv::Mat1f> channels(num_channels);
	split(m3, channels);

	for(int c = 0; c < num_channels; c++)
	{
		channels[c] = channels[c] + m1;
	}

	cv::Mat3f F;
	cv::merge(channels, F);
	return F;
}

inline cv::Mat pow(float a, cv::Mat src)
{
	cv::Mat dst;
	cv::exp(log(a) * src, dst);
	return dst;
}

template <class TField>
inline cv::Mat_<TField> cumsum(cv::Mat_<TField> src, int direction, int offset)
{
	int h = src.rows;
	int w = src.cols;

	if(direction == kHORIZONTAL)
	{
		cv::Mat_<TField> dst = cv::Mat_<TField>::zeros(h, w + offset);
		for (int y = 0; y < h; y++)
		{
			for (int x = offset; x < w + offset; x++)
			{
				dst.template at<TField>(y, x) = src.template at<TField>(y, x - offset) + dst.template at<TField>(y, x - offset);
			}
		}
		return dst;
	}

	else if(direction == kVERTICAL)
	{
		cv::Mat_<TField> dst = cv::Mat_<TField>::zeros(h + offset, w);
		for (int y = offset; y < h + offset; y++)
		{
			for (int x = 0; x < w; x++)
			{
				dst.template at<TField>(y, x) = src.template at<TField>(y - offset, x) + dst.template at<TField>(y - offset, x);
			}
		}
		return dst;
	}

	return cv::Mat_<TField>();
}

template <class TMat>
inline TMat diff(TMat src, int direction, bool advanced_filter)
{
	TMat result;

	float delta = 0;
	int ddepth = -1;

	if(direction == kHORIZONTAL)
	{
		cv::Point anchor = cv::Point(-1, -1);
		cv::Mat1f kernel;
		if(!advanced_filter)
		{
			kernel = cv::Mat1f::ones(1, 3);
			kernel(0, 0) = 0;
			kernel(0, 1) = 1;
			kernel(0, 2) = -1;
		}
		else
		{
			kernel = cv::Mat1f(1,5);
			kernel(0, 0) = 1;
			kernel(0, 1) = -8;
			kernel(0, 2) = 0;
			kernel(0, 3) = 8;
			kernel(0, 4) = -1;
			kernel /= 12;
		}

		filter2D(src, result, ddepth, kernel, anchor, delta, cv::BORDER_DEFAULT);
	}
	else if(direction == kVERTICAL)
	{
		cv::Point anchor = cv::Point(-1, 0);
		cv::Mat1f kernel;

		if(!advanced_filter)
		{
			kernel = cv::Mat1f::ones(2, 1);
			kernel(0, 0) = -1;
			kernel(1, 0) = 1;
		}
		else
		{
			kernel = cv::Mat1f(5,1);
			kernel(0, 0) = 1;
			kernel(1, 0) = -8;
			kernel(2, 0) = 0;
			kernel(3, 0) = 8;
			kernel(4, 0) = -1;
			kernel /= 12;
		}
		
		filter2D(src, result, ddepth , kernel, anchor, delta, cv::BORDER_DEFAULT);
	}

	return result;
}

template <class TValue>
inline int findIndexOf(TValue value, TValue * collection, int range_min, int range_max, int offset, int default_value)
{
	for(int i = range_min; i < range_max; i++)
	{
		if(collection[i + offset] > value)
		{
			return i;
		}
	}
	return default_value;
}

template <class TMat>
inline TMat copyMakeBorder(TMat src, int border_top, int border_bottom, int border_left, int border_right, int border_type)
{
	TMat result;
	copyMakeBorder(src, result, border_top, border_bottom, border_left, border_right, border_type);
	return result;
}


inline cv::Mat3f im2float(cv::InputArray src)
{
	cv::Mat3f result;
	src.getMat().convertTo(result, CV_32F, 1 / 255.);
	return result;
}

inline cv::Mat3b float2im(cv::InputArray src)
{
	cv::Mat3b result;
	src.getMat().convertTo(result, result.type(), 255.0, -0.5);
	return result;
}

template <class TValue>
inline cv::Mat_<TValue> manhattan(cv::Mat m, bool normalize)
{
	cv::Mat_<TValue> result;
	int num_channels = m.channels();
	std::vector< cv::Mat_<TValue> > channels(num_channels);
	cv::split(m, channels);
	result = abs(channels[0].clone());

	for(int c = 1; c < num_channels; c++)
		result += abs(channels[c]);

	if(normalize)
		result /= num_channels;

	return result;
}

// Applies a given mask M to an input image I, asigning 
//     M(X,Y)  > threshold : R(X,Y) = I(X,Y) * M(X,Y)
//     M(X,Y) <= threshold : R(X,Y) = default_value
template <class TValue>
cv::Mat_<TValue> applyMask(const cv::Mat_<TValue>&  image, const cv::Mat1f&  mask, float threshold, TValue default_value)
{
	cv::Mat_<TValue> selection(image.rows, image.cols);
	//std::cout<<"bpx1.5"<<std::endl;
	#pragma omp parallel for
	for(int y = 0; y < image.rows; ++y)
	{
		for(int x = 0; x < image.cols; ++x)
		{
			if(mask.at<float>(y, x) < threshold){
				//std::cout<<"bpx1.6"<<std::endl;
				selection.template at<TValue>(y, x) = default_value;
			}
			else{
				//std::cout<<"bpx1.7"<<std::endl;
				selection.template at<TValue>(y, x) = image.template at<TValue>(y, x) * mask.at<float>(y, x);
			}
			//std::cout<<"bpx1.8"<<std::endl;
		}
		//std::cout<<"bpx1.9"<<std::endl;
	}
	//std::cout<<"bpx1.10"<<std::endl;
	return selection;
}

// Divides each channel j of the image J by the corresponding value of the confidence image G
//     R(x,y)[j] = J(x,y)[j] / G(x,y)
template<class TValue>
inline cv::Mat_<TValue> normalize(cv::Mat_<TValue> J, cv::Mat1f G)
{
	int num_channels = J.channels();

	// Apply the normalization function.
	std::vector<cv::Mat1f> channels(num_channels);
	cv::split(J, channels);
	for(int c = 0; c < num_channels; c++)
	{
		cv::divide(channels[c], G, channels[c]);
	}

	// Merge the result matrix.
	cv::Mat_<TValue> F;
	cv::merge(channels, F);

	return F;
}
