function FakePix2Pix = realtimePix2Pix()
clc
% Define conda activate and python commands
condaPath = 'C:\Users\menon\miniconda3\Scripts\activate';
envName = 'pix2pix';
%processScriptPath = 'C:\Users\menon\pytorch-CycleGAN-and-pix2pix\realtime_pix2pix.py';
%pythonScriptPath = 'C:\Users\menon\pytorch-CycleGAN-and-pix2pix\test.py';
pythonScriptPath = 'C:\Users\menon\pytorch-CycleGAN-and-pix2pix\realtime_test.py';
pythonScriptArgs = ['--dataroot C:\Users\menon\pytorch-CycleGAN-and-pix2pix\datasets\realtime ' ...
    '--name C:\Users\menon\pytorch-CycleGAN-and-pix2pix\checkpoints\combined_plant_0bp_noBlurryRemoval ' ... 
    '--model pix2pix --direction BtoA --input_nc 1 --output_nc 1 --num_test 1 --epoch latest --gpu_ids -1'];
%20230708_minlap_gan_pro ' ...
% Create the full command
cmd1 = ['"' condaPath '" ' envName ' && python "' processScriptPath '"'];
cmd2 = ['"' condaPath '" ' envName ' && python "' pythonScriptPath '" ' pythonScriptArgs];

% Execute the processing command
[status1, cmdout1] = system(cmd1, '-echo');

% Check if command execution was successful
if status1 == 0
    disp('Processing command executed successfully.')
    disp(cmdout1)
else
    disp('Processing command execution failed.')
    disp(cmdout1)
end

% Execute the testing command
[status2, cmdout2] = system(cmd2, '-echo');

% Check if command execution was successful
if status2 == 0
    disp('Testing command executed successfully.')
    disp(cmdout2)
else
    disp('Testing command execution failed.')
    disp(cmdout2)
end

% Define fake image file path
imgFilePath = 'C:\Users\menon\pytorch-CycleGAN-and-pix2pix\checkpoints\combined_plant_0bp_noBlurryRemoval\test_latest\images\1_fake_B.png';

% Read the image file into a matrix
FakePix2Pix = rgb2gray(imread(imgFilePath));
