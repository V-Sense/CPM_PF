#!/bin/bash
test_dir="Cotton"
output_dir=$test_dir"/Results_CPMPF_DISP_center_col/"
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi
../../build/bin/CPMPF_DISP $test_dir/Input_images input_Cam .png 4 9 ver -img_skip 9 -o $output_dir -save_intermediate -output_CPM $output_dir -output_PF $output_dir -HCI -img_idx_width 3
