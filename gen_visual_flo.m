clear all;
close all;

addpath(genpath('~/Developer/MATLAB/epicflow_with_Matlab_wrapper/epicflow_with_Matlab_wrapper/utils/flow-code-matlab'));
addpath('~/Developer/Datasets/LF_dataset/HCI_4D_Light_Field_Benchmark_tools/matlab_tools/lib');

currentDirFileList = dir( pwd );

for i = 1 : length( currentDirFileList )
	if( isequal( currentDirFileList( i ).name, '.' )||...
        isequal( currentDirFileList( i ).name, '..')||...
	    currentDirFileList( i ).isdir)     %如果是文件夹则跳过
        continue;
    end

    floFilePath = fullfile(pwd, currentDirFileList( i ).name);
    flow = flowToColor(readFlowFile(floFilePath));
    saveImgPath = fullfile(pwd, [currentDirFileList( i ).name '.png']);
    imwrite(flow,saveImgPath);
end
