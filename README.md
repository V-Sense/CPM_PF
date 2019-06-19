# Optical Flow based Depth Map Estimation from Light Fields

This code implements the coarse to fine patch match [1] followed by the permeability filter described in [2]. It can be used on video sequences to estimate optical flow as well as sequence of images extracted from light fields to estimate disparity [3].

If you use or adapt this code in your work (either as a stand-alone tool or as a component of any algorithm), you need to cite the appropriate papers [1,2,3].

## Compilation 

The code was tested to work in Linux (Ubuntu 16.04). Clone and compile this repository with:

```
git clone git@github.com:V-Sense/CPM_PF.git
cd CPM_PF
mkdir build
cd build
cmake ..
make -j4
```

### Dependencies

- GCC 5.4
- CMake 3.10.2
- LAPACK
- OpenCV 3.4.1 with opencv_contrib repo

This program was tested on 64 bit Ubuntu 16.04 LTS with Intel(R) Core(TM) i7-6700K CPU @ 4.00GHz.

## Usage

The compilation of this code will produce two executables, CMPPF_FLOW and CPMPF_DISP, dedicated to video optical flow estimation and light field disparity estimation respectively. Note the for the light field disparity estimation, only a row or a column of the light field is expected as input, as described in [3]. It is also assumed that the input images are indexed from left to right in a row, and top to bottom in a column, if not the sign of the disparity will have to be corrected.

Below is an excerpt of the manual displayed when using the `-h` or `-help` flag. See the testing section for mode detailed usage examples.

```
  ./CPMPF_* <input_image_folder> <img_pre> <img_ext> <start_idx> <nb_imgs> <ang_dir> [options]
Mandatory parameters:
    <input_image_folder>                       path to input images
    <img_pre>                                  image name prefix. If there is no image prefix use none or ''
    <img_ext>                                  extension of the image file format, can be anything supported by OpenCV
    <start_idx>                                index of the first image to process. Note that the image index is expected to be located in between the image prefix and the image optional suffix of extension
    <nb_imgs>                                  number of images to process
Options:
  Additional image naming options:
    -img_idx_width                             length of the image index number
    -img_skip                                  index number skip to use if the images indices are not consecutive (i.e. different from 1)
    -img_suf                                   suffix to add before image format extension
  Output result folders:
    -o, -output_VR                             set the final output folder (after variational refinement), default is <input_image_folder>
    -save_intermediate                         use this flag to save results from CPM and PF steps, use the following flages to set the output folders
    -output_CPM                                set the output folder for the Coarse-to-fine Patchmatch step, default is <input_image_folder>
    -output_PF                                 set the output folder for the Permeability Filter steps (both spatial and temporal, default is <input_image_folder>
  Predefined parameters:
    -Sintel                                    parameters for the MPI-Sintel dataset
    -HCI                                       parameters for the HCI synthetic light field dataset
    -Stanford                                  parameters for the Stanford gantry light field dataset
    -TCH                                       parameters for the Technicolor camera array light field dataset
```



## Testing

To test the executables, we provide test data and scripts in folder ./testing.
The testing/Sintel folder can be used to test the CPMPF_FLOW, while the other folders can be used to test CPMPF_DISP.

Simply unzip the test images in the folder where they are located, and run the .sh scripts. Note that the different command lines in the .sh scripts also showcase usecases corresponding to different image naming conventions.

If you use these test data in your own work, you need to cite the appropriate papers provided in the dataset_reference.txt files.

## References

1. Hu, Y., Song, R. and Li, Y., 2016. **Efficient coarse-to-fine patchmatch for large displacement optical flow.** In *Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition* (pp. 5704-5712).
2. Schaffner, M., Scheidegger, F., Cavigelli, L., Kaeslin, H., Benini, L. and Smolic, A., 2018. **Towards Edge-Aware Spatio-Temporal Filtering in Real-Time.** In *IEEE Transactions on Image Processing*, *27*(1), pp.265-280.
3. Chen, Y., Alain, M. and Smolic, A., 2017. **Fast and Accurate Optical Flow based Depth Map Estimation from Light Fields.** In *Irish Machine Vision and Image Processing Conference (IMVIP)*.
4. HCI 4D Light Field Dataset http://hci-lightfield.iwr.uni-heidelberg.de/
5. MPI Sintel Flow Dataset http://sintel.is.tue.mpg.de/
6. Middlebury Optical Flow Dataset http://vision.middlebury.edu/flow/