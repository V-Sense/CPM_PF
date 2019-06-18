#!/bin/bash
test_dir="LegoKnights"
output_dir=$test_dir"/Results_CPMPF_DISP_center_col/"
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi
../../build/bin/CPMPF_DISP $test_dir/Input_images SAI_ .png 7 5 ver -img_suf _09 -o $output_dir -save_intermediate -output_CPM $output_dir -output_PF $output_dir -Stanford -img_idx_width 2
