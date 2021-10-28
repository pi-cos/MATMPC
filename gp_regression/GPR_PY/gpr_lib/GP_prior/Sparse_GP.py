"""
Author: Alberto Dalla Libera
Definition of linear and polynomial GP.
"""



import torch
import numpy as np
from . import GP_prior




class Linear_GP(GP_prior.GP_prior):
    """
    Implementation of the GP with linear kernel (dot product covariance).
    f(X) = phi(X)*w, where the w is defined as a gaussian variable (mean_w and Sigma_w)
    """


    def __init__(self, active_dims,
                 mean_init=None, flg_mean_trainable=False, flg_no_mean=False,
                 sigma_n_init=None, flg_train_sigma_n=False,
                 Sigma_function=None, Sigma_f_additional_par_list=None, 
                 Sigma_pos_par_init=None, flg_train_Sigma_pos_par=True,
                 Sigma_free_par_init=None, flg_train_Sigma_free_par=True,
                 flg_offset=False,
                 scale_init=np.ones(1), flg_train_scale=False,
                 name='', dtype=torch.float64, sigma_n_num=None, device=None):
        """
        Initialization of the object:
        - mean_init and Sigma_function define the prior on the model parameters w
        - Sigma_pos_par Sigma_free par are, respectively, the positive and free parameters of the w prior
        - if flg_offset is true a constant feature equal to 1 is added to the input X
        - active_dims define the input dimension interested i.e. the phi function
        """
        # initilize the GP object
        super(Linear_GP, self).__init__(active_dims,
                                        sigma_n_init=sigma_n_init, flg_train_sigma_n=flg_train_sigma_n,
                                        scale_init=scale_init, flg_train_scale=flg_train_scale,
                                        name=name, dtype=dtype, sigma_n_num=sigma_n_num, device=device)
        # check active dims
        if active_dims is None:
            raise RuntimeError("Active_dims are needed")
        self.num_features = active_dims.size
        # save flg_offset (flg_offset=True => ones added to the phi)
        self.flg_offset = flg_offset
        # check mean init
        self.check_mean(mean_init, flg_mean_trainable, flg_no_mean)
        # check Sigma init
        self.check_sigma_function(Sigma_function, Sigma_f_additional_par_list,
                                  Sigma_pos_par_init, flg_train_Sigma_pos_par,
                                  Sigma_free_par_init, flg_train_Sigma_free_par)
    

    def check_mean(self, mean_init, flg_mean_trainable, flg_no_mean):
        if mean_init is None:
            mean_init = np.zeros(1)
            flg_no_mean = True
        self.flg_no_mean = flg_no_mean
        self.mean_par = torch.nn.Parameter(torch.tensor(mean_init, dtype=self.dtype, device=self.device), requires_grad=flg_mean_trainable)


    def check_sigma_function(self, Sigma_function, Sigma_f_additional_par_list,
                             Sigma_pos_par_init, flg_train_Sigma_pos_par,
                             Sigma_free_par_init, flg_train_Sigma_free_par):
        if Sigma_function is None:
            raise RuntimeError("Specify a Sigma function")
        self.Sigma_function = Sigma_function
        self.Sigma_f_additional_par_list = Sigma_f_additional_par_list
        if Sigma_pos_par_init is None:
            self.Sigma_pos_par = None
        else:
            self.Sigma_pos_par = torch.nn.Parameter(torch.tensor(np.log(Sigma_pos_par_init), dtype=self.dtype, device=self.device), requires_grad=flg_train_Sigma_pos_par)
        if Sigma_free_par_init is None:
            self.Sigma_free_par = None
        else:
            self.Sigma_free_par = torch.nn.Parameter(torch.tensor(Sigma_free_par_init, dtype=self.dtype, device=self.device), requires_grad=flg_train_Sigma_free_par)


    def get_phi(self,X):
        """
        Returns the regression matrix associated to the inputs X
        """
        num_samples = X.shape[0]
        if self.flg_offset:
            return torch.cat([X[:,self.active_dims],torch.ones(num_samples,1, dtype=self.dtype, device=self.device)],1)
        else:
            return X[:,self.active_dims]


    def get_Sigma(self):
        """
        Computes the Sigma matrix
        """
        if self.Sigma_pos_par is None:
            return self.Sigma_function( self.Sigma_pos_par, self.Sigma_free_par, *self.Sigma_f_additional_par_list)
        else:
            return self.Sigma_function( torch.exp(self.Sigma_pos_par), self.Sigma_free_par, *self.Sigma_f_additional_par_list)


    def get_mean(self,X):
        """
        Returns phi(X)*w_mean
        """
        if self.flg_no_mean:
            N = X.size()[0]
            return torch.zeros(N, 1, dtype=self.dtype, device=self.device)
        else:
            return torch.matmul(self.get_phi(X), self.mean_par)


    def get_covariance(self, X1, X2=None, flg_noise=False):
        """
        Returns phi(X)^T*Sigma*phi(X)
        """
        # get the parameters variance
        Sigma = self.get_Sigma()
        # get the covariance
        phi_X1 = self.get_phi(X1)
        if X2 is None:
            K_X =  torch.exp(self.scale_log)*torch.matmul(phi_X1, torch.matmul(Sigma, phi_X1.transpose(0,1)))
            # check if we need to add the noise
            if flg_noise & self.GP_with_noise:
                N = X1.size()[0]
                return K_X + self.get_sigma_n_2()*torch.eye(N, dtype=self.dtype, device=self.device)
            else:
                return K_X
        else:
            return torch.exp(self.scale_log)*torch.matmul(phi_X1, torch.matmul(Sigma, self.get_phi(X2).transpose(0,1)))


    def get_diag_covariance(self, X, flg_noise=False):
        """
        Returns the diag of the cov matrix
        """
        # Get the parameters the variance
        Sigma = self.get_Sigma()
        # get the diag of the covariance
        phi_X = self.get_phi(X)
        diag = torch.exp(self.scale_log)*torch.sum(torch.matmul(phi_X, Sigma)*(phi_X), dim=1)
        if flg_noise & self.GP_with_noise:
            return diag + self.get_sigma_n_2()
        else:
            return diag


    def get_parameters(self, X, Y, flg_print=False):
        """
        Returns the posterior estimate of w, the parameters of the regression.
        NB: the parameters returned are correct only if this is the only GP
        """
        # get the mean and the inverse of the kernel matrix using forward
        m_X, _, K_X_inv, _ = self.forward(X)
        Y = Y - m_X
        # get sigma and phi
        Sigma = self.get_Sigma()
        phi_X_T = torch.transpose(self.get_phi(X), 0, 1)
        # get the parameters
        w_hat = torch.matmul(Sigma, torch.matmul(phi_X_T, torch.matmul(K_X_inv, Y)))
        if flg_print:
            print(self.name+' linear parameters estimated: ', w_hat.data)
        return w_hat


    def get_parameters_inv_lemma(self, X, Y, flg_print=False):
        """
        Returns the estimate of w, the parameters of the regression.
        This implementation exploit the sparsity of the covariance and uses 
        the matrix inversion lemma.
        NB: the parameters returned are correct only if this is the only GP
        """
        # get the mean and the phi and the covariance of the parameters
        m_X = self.get_mean(X)
        Y = Y - m_X
        Phi_X = self.get_phi(X)
        Sigma = self.get_Sigma()
        sigma_n_square = self.get_sigma_n_2()
        # return the parameters
        cov = torch.inverse(Sigma) + sigma_n_square*torch.matmul(Phi_X.transpose(0,1), Phi_X)
        cov = torch.inverse(cov)
        w_hat = sigma_n_square*torch.matmul(torch.matmul(cov, Phi_X.transpose(0,1)), Y)
        if flg_print:
            print(self.name+' linear parameters estimated: ', w_hat.data)
        return w_hat




class Poly_GP(Linear_GP):
    """
    GP with polynomial kernel. Implemented extending the Linear_GP
    """
    def __init__(self, active_dims, poly_deg,
                 sigma_n_init=None, flg_train_sigma_n=True,
                 Sigma_function=None, Sigma_f_additional_par_list=None, 
                 Sigma_pos_par_init=None, flg_train_Sigma_pos_par=True,
                 Sigma_free_par_init=None, flg_train_Sigma_free_par=True,
                 flg_offset=True,
                 scale_init=np.ones(1), flg_train_scale=False,
                 name='', dtype=torch.float64, sigma_n_num=None, device=None):
        # initilize the mean parameters (no mean considered)
        mean_init=None
        flg_mean_trainable=False
        flg_no_mean=True
        # initialize the linear model
        super(Poly_GP, self).__init__(active_dims=active_dims,
                                      mean_init=mean_init, flg_mean_trainable=flg_mean_trainable, flg_no_mean=flg_no_mean,
                                      sigma_n_init=sigma_n_init, flg_train_sigma_n=flg_train_sigma_n,
                                      Sigma_function=Sigma_function, Sigma_f_additional_par_list=Sigma_f_additional_par_list, 
                                      Sigma_pos_par_init=Sigma_pos_par_init, flg_train_Sigma_pos_par=flg_train_Sigma_pos_par,
                                      Sigma_free_par_init=Sigma_free_par_init, flg_train_Sigma_free_par=flg_train_Sigma_free_par,
                                      flg_offset=flg_offset,
                                      scale_init=np.ones(1), flg_train_scale=False,
                                      name=name, dtype=dtype, sigma_n_num=sigma_n_num, device=device)
        # save the poly deg
        self.poly_deg = poly_deg
        # save the scaling parameter of the poly transformation
        self.scale_log = torch.nn.Parameter(torch.tensor(np.log(scale_init), dtype=self.dtype, device=self.device),
                                            requires_grad=flg_train_scale)
    

    def get_covariance(self, X1, X2=None, flg_noise=False):
        """
        Returns the linear covariance raised to self.poly_deg
        """
        return torch.exp(self.scale_log)*(super(Poly_GP, self).get_covariance(X1, X2, flg_noise))**self.poly_deg


    def get_diag_covariance(self, X, flg_noise=False):
        return torch.exp(self.scale_log)*(super(Poly_GP, self).get_diag_covariance(X,flg_noise=flg_noise))**self.poly_deg


    def get_parameters(self, X, Y, flg_print=False):
        raise NotImplementedError()


    def get_parameters_inv_lemma(self, X, Y, flg_print=False):
        raise NotImplementedError()




def get_SOR_GP(exact_GP_object):
    """
    function that returns a Subset of Regressors GP, given a GP object.
    This model is a low-rank approximation of an exact GP model.
    Given a GP model with kernel function k(x_1,x_2), SOR approximate its covariance
    with k_SOR(x_1,x_2) = k(x_1,U) K(U,U)^-1 k(U,x_2), where:
    - U = {set of inducing inputs}
    - K(U,U) = kernel matrix associated to the inducing inputs 
    """


    # create the SOR_GP class dynamically
    class SOR_GP(type(exact_GP_object)):
        """SOR_GP"""
        def __init__(self, exact_GP_object):
            """Initialize the object inheriting all the exact_GP_object parameters"""
            # initialize the GP object randomly
            GP_prior.GP_prior.__init__(self,
                                       active_dims=[0],
                                       sigma_n_init=None,
                                       flg_train_sigma_n=False,
                                       scale_init=np.ones(1),
                                       flg_train_scale=False,
                                       name='',
                                       dtype=torch.float64,
                                       sigma_n_num=None,
                                       device=torch.device('cpu'))
            # assign all the variables of the exact_GP_object
            self.__dict__ = exact_GP_object.__dict__
            self.name = 'SOR_GP '+self.name


        def init_inducing_inputs(self, inducing_inputs, flg_train_inducing_inputs=False):
            """set the set U, a matrix of dimension (num_inducing_inputs, num_feature)
               which represent the initial value of the inducing inputs.
               If flg_train_inducing_inputs=True U is considered as a trainable hyperparameter
            """
            self.U = torch.nn.Parameter(torch.tensor(inducing_inputs, dtype=self.dtype, device=self.device),
                                        requires_grad=flg_train_inducing_inputs)


        def set_inducing_inputs_from_data(self, X, Y, threshold, flg_trainable):
            """
            Set the inducing inputs with an online procedure
            """
            print('\nSelection of the inducing inputs...')
            # get number of samples
            num_samples = X.shape[0]
            # init the set of inducing inputs with the first sample
            self.U = torch.nn.Parameter(X[0:1,:], requires_grad=flg_trainable)
            inducing_inputs_indices = [0]
            # iterate all the samples
            for sample_index in range(2,num_samples):
                # get the estimate 
                # _, var, _ = self.get_SOR_estimate(X[:sample_index-1,:], Y[:sample_index-1,:], X[sample_index:sample_index+1,:])
                _, var, _ = self.get_estimate(X[inducing_inputs_indices,:], Y[inducing_inputs_indices,:], X[sample_index:sample_index+1,:])
                # check 
                # print('torch.sqrt(var)',torch.sqrt(var))
                if torch.sqrt(var)>threshold:
                    self.U.data = torch.cat([self.U.data,X[sample_index:sample_index+1,:]],0)
                    inducing_inputs_indices.append(sample_index)
            print('Shape of the inducing inputs selected:',self.U.shape)
            return inducing_inputs_indices


        def get_SOR_alpha(self, X, Y):
            """
            Returns the coefficients that defines the SOR posterior distribution
            
            inputs:
            - X = training inputs
            - Y = training outputs

            outputs:
            - alpha = vector defining the SOR posterior distribution
            - m_X = prior mean of X
            - Sigma = inverse of (K_UU + sigma_n^-2*K_UX*K_XU)
            """
            # get the mean and the phi and the covariance of the parameters
            m_X = self.get_mean(X)
            Y = Y - m_X
            K_XU = self.get_covariance(X,self.U)
            K_UU = self.get_covariance(self.U)
            sigma_n_square_inv = 1/self.get_sigma_n_2()
            sigma_n_square = self.get_sigma_n_2()
            # return the parameters
            Sigma_inv = K_UU + sigma_n_square_inv*torch.matmul(K_XU.transpose(0,1), K_XU)
            U_Sigma_inv = torch.cholesky(Sigma_inv, upper=True)
            Sigma = torch.cholesky_inverse(U_Sigma_inv, upper=True)
            
            SOR_alpha = sigma_n_square_inv*torch.matmul(torch.matmul(Sigma, K_XU.transpose(0,1)), Y)
            return SOR_alpha, m_X, Sigma


        def get_SOR_estimate_from_alpha(self, X, X_test, SOR_alpha, m_X, Sigma=None):
            """
            Compute the SOR posterior distribution in X_test, given the alpha vector.
            
            input:
            - X = training input locations (used to compute alpha)
            - X_test = test input locations
            - SOR_alpha = vector of coefficients defining the SOR posterior
            - m_X = prior mean of X
            - Sigma = inverse of (K_UU + sigma_n^-2*K_UX*K_XU)
    
            output:
            - Y_hat = posterior mean
            - var = diagonal elements of the posterior variance
            
            If Sigma_inv is given the method returns also the confidence intervals (variance of the gaussian)
            """
            # get covariance and prior mean
            K_X_test_U = self.get_covariance(X_test, self.U)
            m_X_test = self.get_mean(X_test)
            # get the estimate 
            Y_hat = m_X_test + torch.matmul(K_X_test_U, SOR_alpha)
            # if Sigma_inv is given compute the confidence intervals
            if Sigma is not None:
                var = torch.sum(torch.matmul(K_X_test_U, Sigma)*(K_X_test_U), dim=1)
                # var = self.get_sigma_n_2()*torch.sum(torch.matmul(K_X_test_U, Sigma)*(K_X_test_U), dim=1)
            return Y_hat, var
         

        def get_SOR_estimate(self, X, Y, X_test, flg_return_Sigma=False):
            """
            Returns the SOR posterior distribution in X_test, given the training samples X Y and the inducing inputs self.U.
            
            input:
            - X = training input
            - Y = training output
            - X_test = test input
    
            output:
            - Y_hat = mean of the test posterior
            - var = diagonal elements of the variance posterior
            - SOR_alpha = coefficients defining the SOR posterior
            - m_X = prior mean of the training samples
            - Sigma = (K_UU + sigma_n^-2*K_UX*K_XU)
            """
            # get the coefficent and the mean
            SOR_alpha, m_X, Sigma = self.get_SOR_alpha(X, Y)
            # get the estimate and the confidence intervals
            Y_hat, var = self.get_SOR_estimate_from_alpha(X, X_test, SOR_alpha, m_X, Sigma=Sigma)
            # return the opportune values
            if flg_return_Sigma:
                return Y_hat, var, SOR_alpha, m_X, Sigma
            else:
                return Y_hat, var, SOR_alpha


        def SOR_forward(self, X):
            """
            Returns the elements that define the likelihood distribution of the model
            when considering the SOR approximation
            
            input:
            - X = training inputs (X has dimension [num_samples, num_features])
            
            output:
            - m_X = mean
            - K_X = None
            - K_X_inv = inverse of (K_XU*K_UU^-1*K_UX+sigma_n^2)
            - log_det = (K_XU*K_UU^-1*K_UX+sigma_n^2)
            """
            # get the mean
            m_X = self.get_mean(X)
            # get kernel matrices
            K_X = None
            N = X.shape[0]
            K_UU = self.get_covariance(self.U)
            K_XU = self.get_covariance(X, self.U)
            sigma_n_square_inv = 1/self.get_sigma_n_2()
            # compute the K_UU^-1 logdet
            U_K_UU = torch.cholesky(K_UU, upper=True)
            K_UU_inv_log_det = -2*torch.sum(torch.log(torch.diag(U_K_UU)))
            # compute Sigma
            Sigma_inv = K_UU + sigma_n_square_inv*torch.matmul(K_XU.transpose(0,1), K_XU)
            # compute Sigma inverse and logdet
            U_Sigma_inv = torch.cholesky(Sigma_inv, upper=True)
            Sigma_inv_log_det = 2*torch.sum(torch.log(torch.diag(U_Sigma_inv)))
            Sigma = torch.cholesky_inverse(U_Sigma_inv, upper=True)
            # compute K_X_inv
            K_X_inv = sigma_n_square_inv*torch.eye(N, dtype=self.dtype, device=self.device)
            K_X_inv -= sigma_n_square_inv**2*torch.matmul(K_XU, torch.matmul(Sigma, K_XU.transpose(0,1)))
            # compute the log_det
            log_det = N*torch.log(self.get_sigma_n_2()) + K_UU_inv_log_det + Sigma_inv_log_det
            return m_X, K_X, K_X_inv, log_det


        def fit_SOR_model(self,trainloader=None, 
                          optimizer=None, criterion=None,
                          N_epoch=1, N_epoch_print=1,
                          f_saving_model=None, f_print=None):
            """
            Optimize the SOR model hyperparameters
    
            input:
            - trainloader = torch train loader object
            - optimizer = torch optimizer object
            - criterion = loss function
            - N_epoch = number of epochs
            - N_epoch_print = number of epoch between print two prints of the current loss and model parameters
            - f_saving_model = customizable function that save the model
            - f_print_model = customizable function that print the model (eventually with performance)
            """
            # print initial parametes and initial estimates
            print('\nInitial parameters:')
            self.print_model()  
            # iterate over the training data for N_epochs
            for epoch in range(0,N_epoch):
                # initialize loss grad and counter
                running_loss = 0.0
                N_btc = 0
                optimizer.zero_grad()
                #iterate over the training set
                for i, data in enumerate(trainloader, 0):
                    # get the training data
                    inputs, labels = data
                    # zero the parameter gradients
                    optimizer.zero_grad()
                    # forward + backward + optimize
                    out_SOR_GP_priors = self.SOR_forward(inputs)
                    loss = criterion(out_SOR_GP_priors, labels)
                    loss.backward(retain_graph=False)
                    optimizer.step()
                    # update the running loss
                    running_loss = running_loss + loss.item()
                    N_btc = N_btc + 1
                # print statistics and save the model
                if epoch%N_epoch_print==0:
                    print('\nEPOCH:', epoch)
                    self.print_model()
                    print('Running loss:', running_loss/N_btc)
                    if f_saving_model is not None:
                        f_saving_model(epoch)
                    if f_print is not None:
                        f_print()
            # print the final parameters
            print('\nFinal parameters:')
            self.print_model()



    # init the object and return
    return SOR_GP(exact_GP_object)