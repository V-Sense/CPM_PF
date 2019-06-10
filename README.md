# Optical Flow based Depth Map Estimation from Light Fields

This code implements the coarse to fine patch match [1] followed by the permeability filter described in [2]. It can be used on video sequences to estimate optical flow as well as sequence of images extracted from light fields to estimate disparity [3].

If you use or adapt this code in your work (either as a stand-alone tool or as a component of any algorithm), you need to cite the appropriate papers [1,2,3].

## Usage

```
  ./CPMPF <input_image_folder> <img_pre> <img_suf> <img_ext> <start_idx> <nb_imgs> <ang_dir> [options]
Options:
    -h, -help                                  print this message
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
  PF parameters:
    -PF_iter                                   number of iterantions for spatial permeability filter
    -PF_lambda                                 lambda para for spatial permeability filter
    -PF_delta                                  delta para for spatial permeability filter
    -PF_alpha                                  alpha para for spatial permeability filter
  Predefined parameters:
    -Sintel                                    set the parameters to the one optimized on (a subset of) the MPI-Sintel dataset
    -HCI                                       set the parameters to the one optimized on (a subset of) the HCI light field dataset
```



## Folder Structure 

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
2. Schaffner, M., Scheidegger, F., Cavigelli, L., Kaeslin, H., Benini, L. and Smolic, A., 2018. **Towards Edge-Aware Spatio-Temporal Filtering in Real-Time.** *IEEE Transactions on Image Processing*, *27*(1), pp.265-280.
3. Chen, Y., Alain, M. and Smolic, A., 2017. **Fast and Accurate Optical Flow based Depth Map Estimation from Light Fields.** In *Irish Machine Vision and Image Processing Conference (IMVIP)*.
4. HCI 4D Light Field Dataset http://hci-lightfield.iwr.uni-heidelberg.de/
5. MPI Sintel Flow Dataset http://sintel.is.tue.mpg.de/
6. Middlebury Optical Flow Dataset http://vision.middlebury.edu/flow/