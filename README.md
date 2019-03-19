# CPM_PF
Coarse to fine Patch Match + Permeability Filter: [Github](https://github.com/V-Sense/CPM_PF)

This code implements the coarse to fine patch match [1] followed by permeability filter described in [2]. It can be used on video sequences to estimate optical flow as well as sequence of images extracted from light fields to estimate disparity [3].

When using this program, please cite following papers:

1. Hu, Y., Song, R. and Li, Y., 2016. **Efficient coarse-to-fine patchmatch for large displacement optical flow.** In *Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition* (pp. 5704-5712).
2. Schaffner, M., Scheidegger, F., Cavigelli, L., Kaeslin, H., Benini, L. and Smolic, A., 2018. **Towards Edge-Aware Spatio-Temporal Filtering in Real-Time.** *IEEE Transactions on Image Processing*, *27*(1), pp.265-280.
3. Chen, Y., Alain, M. and Smolic, A., 2017. **Fast and Accurate Optical Flow based Depth Map Estimation from Light Fields.** In *Irish Machine Vision and Image Processing Conference (IMVIP)*.

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
git clone git@github.com:V-Sense/CPM_PF.git
cd CPM_PF
mkdir build
cd build
cmake ..
make -j4
```

*Please note*: Folders inside "bin" have to be created manually and copy "utils" folder from the **Datasets** downloaded below before running following the usage:

```
CPM_PF
 |--CPM_Tip2017Mod
 |--PFilter
 |--#other files#
 |--(build)
      |--#other folders and files#
      |--bin
          |--CPMPF
          |--(inputImages)
          |--(outputMatches)
          |--(outputFlows)
          |--(refinedOutputFlows)
          |--(utils)
```

*inputImages* contains the input image folders;

*outputMatches* contains the intermedia results from CPM;

*outputFlows* contains the dense flow results from CPMPF (similiar to [2]s results);

*refineOutput* contains the CPMPF results filtered by our IMVIP final step filtering.

*utils* contains tools to read/write and render .flo/.pfm files, and generate disparity/depth for HCI data.

## Visualise results

To visualize .flo & .pfm files, add "utils" folder path to matlab first

```
addpath(genpath('./bin/utils/'));
```

Then use "gen_visual_flo.m" "gen_disp.m" "gen_depth.m" to generate visualised results.

To generate depth or compute MSE/RMSE for HCI dataset, disparity groundtruth "gt_disp_lowres.pfm" and parameter file "LF.mat" has to be placed porperly.

## Dataset

Several sequences from HCI 4D Light Field Dataset [4] and uploaded them with some sequences from MPI Sintel Dataset [5] and were used in [2]. The flow visualising matlab tool in the "utils" folder is from [6]. Please cite proper papers if using related resources. Results from our implementation are also provided.  Please find the links below:

| Sintel & HCI Light Field |
| :-------------: |
|       [3.9 GB](https://v-sense.scss.tcd.ie/Datasets/Sintel_HCI.zip)            |

Run the program with Sintel/ambush_3 images for example:

```./CPMPF &lt;input_image_folder&gt; &lt;CPM_match_folder&gt; &lt;CPMPF_flow_folder&gt; &lt;refined_CPMPF_flow_folder&gt; [options]
./CPMPF inputImages/Sintel/ambush_3/ outputMatches/Sintel/ambush_3/ outputFlows/Sintel/ambush_3/ refinedOutputFlows/Sintel/ambush_3 -sintel
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