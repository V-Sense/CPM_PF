#!/bin/bash
test_dir="1"
output_dir=$test_dir"/Results_CPMPF_FLOW/"
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi
../../build/bin/CPMPF_FLOW $test_dir/Input_images frame_ .png 1 4 -o $output_dir -save_intermediate -write_color_png -output_CPM $output_dir -output_PF $output_dir -Sintel
