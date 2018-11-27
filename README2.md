# CPM_PF
Coarse to fine Patch Match + Permeability Filter: https://github.com/V-Sense/CPM_PF.

This program can be used on a certain sequence of images, for example, one angular diminsion of a Light Field image.

## Usage


```
./CPMPF <input_image_folder> <CPM_match_folder> <CPMPF_flow_folder> <refined_CPMPF_flow_folder> [options]
    options:
    	-h, -help                 print this message
   	  
   	  CPM parameters:
        -m, -max                  outlier handling maxdisplacement threshold
        -t, -th                   froward and backward consistency threshold
        -c, -cth                  matching cost check threshold
      
      PF parameters:
        -i, -iter                 number of iterantions for spatial permeability filter
        -l, -lambda               lambda para for spatial permeability filter
        -d, -delta                delta para for spatial permeability filter
        -a, -alpha                alpha para for spatial permeability filter
      
      predefined parameters:
        -sintel                   set the parameters to the one optimized on (a subset of) the MPI-Sintel dataset
        -hcilf                    set the parameters to the one optimized on (a subset of) the HCI light field dataset
```



## Folder Structure 

The code was tested to work in Linux (Ubuntu 16.04). Clone and compile this repository with:

```
git clone git@github.com:chenyangjamie/CPM_PF.git
cd CPM_PF
mkdir build
cd build
cmake ..
make -j4
```

*Please note*: Folders inside "bin" have to be created manually before running following the usage:

```
CPM_PF
 |--CPM_Tip2017Mod
 |--PFilter
 |--#other files#
 |--(build)
      |--#other folders and files#
      |--bin
          |--(inputImages)
          |--(outputMatches)
          |--(outputFlows)
          |--(refinedOutputFlows)
          |--(utils)
```

*inputImages* contains the input image folders;

*outputMatches* contains the intermedia results from CPM;

*outputFlows* contains the dense flow results from CPMPF (similiar to tip17's results);

*refineOutput* contains the CPMPF results filtered by our IMVIP final step filtering.

*utils* contains tools to read/write and render .flo/.pfm files, and generate disparity/depth for HCI data.

## Visulise results

To visualize .flo & .pfm files, add "utils" folder path to matlab first

```
addpath(genpath('./bin/utils/'));
```

Then use "gen_visual_flo.m" "gen_disp.m" "gen_depth.m" to generate visualised results.

To generate depth or compute MSE/RMSE for HCI dataset, disparity groundtruth "gt_disp_lowres.pfm" and parameter file "LF.mat" has to be placed porperly.

## Dataset

Several sequences from HCI light fields and uploaded them with some Sintel sequences used in tip'17 paper. Results from our implementation are also provided.  Please find the links below:

| Sintel & HCI Light Field |
| :-------------: |
|       [3.9 GB](https://v-sense.scss.tcd.ie/Datasets/Sintel_HCI.zip)            |

## Dependency

- GCC 5.4
- CMake 3.10.2
- OpenCV 3.4.1 with opencv_contrib repo

This program is tested on 64 bit Ubuntu 16.04 LTS with Intel(R) Core(TM) i7-6700K CPU @ 4.00GHz.

## Known Issues

* **[Solved]** Variational C functions will meet problem with *uint16* images. Use matlab function **im2uint8()** to convert them to *uint8*. 

## References

When using this program, please cite following papers:

Schaffner, M., Scheidegger, F., Cavigelli, L., Kaeslin, H., Benini, L. and Smolic, A., 2018. **Towards Edge-Aware Spatio-Temporal Filtering in Real-Time.** *IEEE Transactions on Image Processing*, *27*(1), pp.265-280.

Chen, Y., Alain, M. and Smolic, A., 2017. **Fast and Accurate Optical Flow based Depth Map Estimation from Light Fields.** In *Irish Machine Vision and Image Processing Conference (IMVIP)*.

