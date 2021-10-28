restoredefaultpath; clear all; close all; clc;

%% settings

xGP = 5; % Dimension of training points for the GP
yGP = 2; % # of target of the GP

data_file = "black_box_ip_sim"; % placed in GPR_PY/data     % black_box_ip_sim,grey_box_ip_sim
res_folder = "blackbox_s01"; % placed in GPR_PY/results_GP_ID/  % blackbox, greybox

sigma_n = 0.01;


%% creating command and run regression procedure

py_string = "GP_ID.py";
py_string = py_string + " -dim_GP_input "+num2str(xGP);
py_string = py_string + " -dim_GP_output "+num2str(yGP);
py_string = py_string + " -csv_data data/"+data_file+".csv";
py_string = py_string + " -saving_path results_GP_ID/"+res_folder+"/";
py_string = py_string + " -flg_norm True -sigma_n_num "+num2str(sigma_n);


cd GPR_PY
disp('Running GP Regression..')
pyrunfile(py_string)
disp('.. done!')
cd ..