"""
Author: Alberto Dalla Libera
GPR-based identification script
"""

import torch
import torch.utils.data
import numpy as np
import gpr_lib.Utils as utils
import gpr_lib.GP_prior.Stationary_GP as Stat_GP
import gpr_lib.Loss.Gaussian_likelihood as loss
import matplotlib.pyplot as plt
import pandas as pd
import sys
import argparse

p = argparse.ArgumentParser('FP GP ID')
p.add_argument('-dim_GP_input',
               type=int,
               default=None,
               help='Dimension of the GP input')
p.add_argument('-dim_GP_output',
               type=int,
               default=None,
               help='Dimension of the GP output')
p.add_argument('-csv_data',
               type=str,
               default='',
               help='Path to the csv file with training data')
p.add_argument('-saving_path',
               type=str,
               default='',
               help='saving path')
p.add_argument('-flg_norm',
               type=bool,
               default=False,
               help='If true normalize data')
p.add_argument('-downsampling_rate',
               type=int,
               default=1,
               help='Downsampling rate')
p.add_argument('-flg_SOD',
               type=bool,
               default=False,
               help='If true apply SOD')
p.add_argument('-threshold_SOD',
               type=float,
               default=1.,
               help='SOD threshold')
p.add_argument('-flg_cuda',
               type=bool,
               default=False,
               help='If true use cuda')
p.add_argument('-sigma_n_num',
               type=float,
               default=None,
               help='Manual regularization (if None not applied)')
p.add_argument('-num_data',
               type=int,
               default=-1,
               help='Select the first num_data (if None use all the data)')
p.add_argument('-N_epoch',
               type=int,
               default=1501,
               help='Number of training epoch')
p.add_argument('-N_epoch_print',
               type=int,
               default=250,
               help='Epochs between prints during training')
p.add_argument('-batch_size',
               type=int,
               default=-1,
               help='Size of the training batch')
p.add_argument('-lr',
               type=float,
               default=0.01,
               help='Learning rate')
p.add_argument('-flg_plot',
               type=bool,
               default=False,
               help='If true plot estimates')
locals().update(vars(p.parse_known_args()[0]))

# ---- Set properties ----
print('\n---- Set properties ----')
# Set type
dtype = torch.float64
# set device
if flg_cuda:
    device = torch.device('cuda')  
else:
    device = torch.device('cpu')


# ---- Load data ----
print('\n---- Load data ----')
# Read data
data = pd.read_csv(csv_data)
if num_data==-1:
    num_data=data.shape[0]
# Get the I/O data
x_names = ['x_'+str(GP_input_dim) for GP_input_dim in range(1,dim_GP_input+1)]
y_names = ['y_'+str(GP_output_dim) for GP_output_dim in range(1,dim_GP_output+1)]
X = data[x_names].values[:num_data,:]
Y = data[y_names].values[:num_data,:]
# Normalize output
norm_coef = np.max(np.abs(Y),0)
num_samples_tr, input_dim = X.shape
num_output = Y.shape[1]
print('num data tr: ', num_samples_tr)
print('input dim: ', input_dim)
print('num_output: ', num_output)
# sys.exit()


# ---- Train models ----
print('\n---- Train models ----')
m_list = []
X_tc = torch.tensor(X, dtype=dtype, device=device)
Y_tc = torch.tensor(Y/norm_coef, dtype=dtype, device=device)
for gp_model in range(dim_GP_output):
    print('\n---- Train GP model '+str(gp_model+1)+' ----')
    # Initialize the model
    print('Initialize the model...')
    active_dims = torch.tensor(range(dim_GP_input))
    # active_dims = np.arange(dim_GP_input)
    lengthscales_init = np.ones(input_dim)
    sigma_n_init = 1.
    m = Stat_GP.RBF(active_dims=active_dims,
                    lengthscales_init=lengthscales_init, flg_train_lengthscales=True,
                    sigma_n_init=sigma_n_init, flg_train_sigma_n=True,
                    mean_init=None, flg_train_mean = False,
                    scale_init=np.ones(1), flg_train_scale=False,
                    name='GP_'+str(gp_model+1), dtype=dtype, sigma_n_num=sigma_n_num, device=device)
    # Train the model
    print('Train the model...')
    N_epoch = N_epoch
    N_epoch_print = N_epoch_print
    if batch_size==-1:
        batch_size = num_data
    # move data in torch
    tr_dataset = torch.utils.data.TensorDataset(X_tc, Y_tc[:, gp_model:gp_model+1])
    trainloader = torch.utils.data.DataLoader(tr_dataset, batch_size=batch_size, shuffle=True)
    # Set optimizer
    optimizer = torch.optim.Adam(m.parameters(), lr=lr)
    criterion = loss.Marginal_log_likelihood()
    m.fit_model(trainloader=trainloader, 
                optimizer=optimizer, criterion=criterion,
                N_epoch=N_epoch, N_epoch_print=N_epoch_print,
                f_saving_model=None, f_print=None)
    m_list.append(m)


# ---- Compute predictions ----
print('\n---- Compute predictions ----')
if flg_SOD:
    print('---- SOD reduction ----')
    # compute inducing inputs
    inducing_inputs_indices_list = []
    for gp_model in range(dim_GP_output):
        inducing_inputs_indices_list.append(m_list[gp_model].get_SOD(X_tc, Y_tc[:, gp_model:gp_model+1],
                                                                     threshold=threshold_SOD*torch.sqrt(m_list[gp_model].get_sigma_n_2()),
                                                                     flg_permutation=False))
    inducing_inputs_indices = list(set(sum(inducing_inputs_indices_list,[])))
    inducing_inputs_indices.sort()
    X_tc_tr = X_tc[inducing_inputs_indices,:]
    Y_tc_tr = Y_tc[inducing_inputs_indices,:]
    print('\nX_tc_tr shape:', X_tc_tr.shape)
    print('Y_tc_tr shape:', Y_tc_tr.shape)
else:
    # get the posterior
    X_tc_tr = X_tc[::downsampling_rate,:]
    Y_tc_tr = Y_tc[::downsampling_rate,:]
np.savetxt(saving_path+'X.csv', X_tc_tr.detach().cpu().numpy(), delimiter=',')
# test models
for gp_model in range(dim_GP_output):
    print('Test GP '+str(gp_model+1)+'...')
    y_hat_tc, _, alpha, m_X, K_X_inv = m_list[gp_model].get_estimate(X_tc_tr, Y_tc_tr[:,gp_model:gp_model+1], X_tc, flg_return_K_X_inv=True)
    # Save the model
    print('Save the model...')
    np.savetxt(saving_path+'alpha_'+str(gp_model+1)+'.csv', norm_coef[gp_model]*alpha.detach().cpu().numpy(), delimiter=',')
    np.savetxt(saving_path+'l_'+str(gp_model+1)+'.csv', np.exp(m_list[gp_model].log_lengthscales_par.detach().cpu().numpy()), delimiter=',')
    np.savetxt(saving_path+'y_hat_'+str(gp_model+1)+'.csv', norm_coef[gp_model]*y_hat_tc.detach().cpu().numpy(), delimiter=',')
    if flg_plot:
        # plot results
        plt.figure()
        plt.grid()
        plt.plot(Y[:,gp_model], label='$y_'+str(gp_model+1)+'$')
        plt.plot(norm_coef[gp_model]*y_hat_tc.detach().cpu().numpy(), label='$\\hat{y}_'+str(gp_model+1)+'$')
        plt.legend()
plt.show()