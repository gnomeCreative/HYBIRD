%% This scripts write the job file for the HPC at the Sheffield University
close all; clear all; clc;

hybird_name = './hybird.exe';       % name of the .exe file
job_file_name = 'jobfile.sh';       % job file name
job_name = 'deg28F';                % job name  on the hpc
cpus_per_task = 4;                  % number of processor
OMP_NUM_THREADS = cpus_per_task;    % must be equal to cpus
mem = '20G';                        % memory required overall
path = './results/flume/deg28_3D/'; % results and error path
output = [path,'output.txt'];       
error = [path,'error.txt'];         
time = '96:00:00';                  % time for running before being stopped (96 hours is the upper limit)
mail = 'a.pasqua@sheffield.ac.uk';  
config_file = 'configFlume_3D_deg28.cfg';   % configuration file of hybird
results_dir = './results/flume';            % result directory
run_name = 'deg28_3D';                      % name of the result folder


% Open the file for writing
fileID = fopen(job_file_name, 'w');

% Write the job file
fprintf(fileID, '#!/bin/bash\n');
fprintf(fileID, '#SBATCH --job-name=%s\n', job_name);
fprintf(fileID, '#SBATCH --nodes=1\n');
fprintf(fileID, '#SBATCH --ntasks-per-node=1\n');
fprintf(fileID, '#SBATCH --ntasks=1\n');
fprintf(fileID, '#SBATCH --cpus-per-task=%i\n', cpus_per_task);
fprintf(fileID, '#SBATCH --mem=%s\n', mem);
fprintf(fileID, '#SBATCH --output=%s\n', output);
fprintf(fileID, '#SBATCH --error=%s\n', error);
fprintf(fileID, '#SBATCH --time=%s\n', time);
fprintf(fileID, '#SBATCH --mail-user=%s\n', mail);
fprintf(fileID, 'export OMP_NUM_THREADS=%i\n', OMP_NUM_THREADS);
fprintf(fileID, '%s -c %s -d %s -n %s\n', hybird_name, config_file, results_dir, run_name);

% Close the file
fclose(fileID);

