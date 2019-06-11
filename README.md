# Optical Flow based Depth Map Estimation from Light Fields

This code implements the coarse to fine patch match [1] followed by the permeability filter described in [2]. It can be used on video sequences to estimate optical flow as well as sequence of images extracted from light fields to estimate disparity [3].

If you use or adapt this code in your work (either as a stand-alone tool or as a component of any algorithm), you need to cite the appropriate papers [1,2,3].

## Usage

```
  ./CPMPF <input_image_folder> <img_pre> <img_ext> <start_idx> <nb_imgs> <ang_dir> [options]
Options:
    -h, -help                                  print this message
  Additional image naming options:
    -img_idx_width                             length of the image index number
    -img_skip                                  index number skip
    -img_suf                                   suffix to add before image format extension
  Output result folders:
    -o, -output_VR                             set the final output folder (after variational refinement), default is <input_image_folder>
    -save_intermediate                         use this flag to save results from CPM and PF steps, use the following flages to set the output folders
    -output_CPM                                set the output folder for the Coarse-to-fine Patchmatch step, default is <input_image_folder>
    -output_PF                                 set the output folder for the Permeability Filter steps (both spatial and temporal, default is <input_image_folder>
  CPM parameters:
    -CPM_max                                   outlier handling maxdisplacement threshold
    -CPM_fbth                                  forward and backward consistency threshold
    -CPM_cth                                   matching cost check threshold
    -CPM_stereo                                stereo flag
    -CPM_nstep                                 number of step giving the final result resolution
  PF:
    Spatial parameters:
    -PF_iter                                   number of iterations
    -PF_lambda                                 lagrangian factor to balance fidelity to the input data
    -PF_delta                                  transition point of the edge-stopping function
    -PF_alpha                                  falloff rate of the edge-stopping function
  VR parameters:
    -VR_alpha                                  smoothness weight
    -VR_gamma                                  gradient constancy assumption weight
    -VR_delta                                  color constancy assumption weight
    -VR_sigma                                  presmoothing of the images
    -VR_niter_outer                            number of outer fixed point iterations
    -VR_niter_inner                            number of inner fixed point iterations
    -VR_niter_solver                           number of solver iterations 
    -VR_sor_omega                              omega parameter of sor method
  Predefined parameters:
    -Sintel                                    parameters for the MPI-Sintel dataset
    -HCI                                       parameters for the HCI synthetic light field dataset
    -Stanford                                  parameters for the Stanford gantry light field dataset
    -TCH                                       parameters for the Technicolor camera array light field dataset
```




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

## Dependencies

- GCC 5.4
- CMake 3.10.2
- OpenCV 3.4.1 with opencv_contrib repo

This program is tested on 64 bit Ubuntu 16.04 LTS with Intel(R) Core(TM) i7-6700K CPU @ 4.00GHz.

## References

1. Hu, Y., Song, R. and Li, Y., 2016. **Efficient coarse-to-fine patchmatch for large displacement optical flow.** In *Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition* (pp. 5704-5712).
2. Schaffner, M., Scheidegger, F., Cavigelli, L., Kaeslin, H., Benini, L. and Smolic, A., 2018. **Towards Edge-Aware Spatio-Temporal Filtering in Real-Time.** In *IEEE Transactions on Image Processing*, *27*(1), pp.265-280.
3. Chen, Y., Alain, M. and Smolic, A., 2017. **Fast and Accurate Optical Flow based Depth Map Estimation from Light Fields.** In *Irish Machine Vision and Image Processing Conference (IMVIP)*.
4. HCI 4D Light Field Dataset http://hci-lightfield.iwr.uni-heidelberg.de/
5. MPI Sintel Flow Dataset http://sintel.is.tue.mpg.de/
6. Middlebury Optical Flow Dataset http://vision.middlebury.edu/flow/