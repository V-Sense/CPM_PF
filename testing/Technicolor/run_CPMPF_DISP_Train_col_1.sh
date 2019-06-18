#!/bin/bash
test_dir="Train"
output_dir=$test_dir"/Results_CPMPF_DISP_center_col/"
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi
../../build/bin/CPMPF_DISP $test_dir/Input_images Train_pr_00001_ .ppm 1 4 ver -img_skip 4 -o $output_dir -save_intermediate -output_CPM $output_dir -output_PF $output_dir -TCH -img_idx_width 2

