"""
Author: Alberto Dalla Libera
High level definition of the GP object.
"""

import torch
import numpy as np





class GP_prior(torch.nn.Module):
    """
    Superclass of each GP model (this class extends torch.nn.Module)
    """


    def __init__(self, active_dims,
                 sigma_n_init=None, flg_train_sigma_n=False,
                 scale_init=np.ones(1), flg_train_scale=False,
                 name='', dtype=torch.float64, sigma_n_num=None, device=torch.device('cpu')):
        """
           Initialize the Module object and set flags regarding:
           - noise
            sigma_n_int = initial noise (None = no noise)
            flg_train_sigma_n = set to true to train the noise parameter
           - sigma_n_num = used to compensate numerical noise
           - active dims = indices of the input selected by the GP model to compute mean and covariance
           - data type
           - device
           - name
           - scaling parameters:
             scale_init = initialization of the scaling parameters
             flg_train_scale = set to true to train the scaling parameter
        """
        # initilize the Module object
        super(GP_prior,self).__init__()
        # set name device and type
        self.name = name
        self.dtype = dtype
        self.device = device
        # active dims
        self.active_dims = torch.tensor(active_dims, requires_grad=False, device=device)
        # set sigma_n_log (the log of the noise standard deviation) 
        if sigma_n_init is None:
            self.GP_with_noise = False
        else:
            self.GP_with_noise = True
            self.sigma_n_log = torch.nn.Parameter(torch.tensor(np.log(sigma_n_init), dtype=self.dtype, device=self.device),
                                                  requires_grad=flg_train_sigma_n)
        # scaling parameters
        self.scale_log = torch.nn.Parameter(torch.tensor(np.log(scale_init), dtype=self.dtype, device=self.device),
                                            requires_grad=flg_train_scale)
        # standard deviation of the rounding errors compensation
        if sigma_n_num is not None:
            self.sigma_n_num = torch.tensor(sigma_n_num, dtype=self.dtype, device=self.device)
        else:
            self.sigma_n_num = torch.tensor(0., dtype=self.dtype, device=self.device)
    

    def to(self, dev):
        """
        Set the device and move the parameters
        """
        # set the new device
        super(GP_prior, self).to(dev)
        self.device = dev
        # move the model parameters in the new device
        self.sigma_n_num = self.sigma_n_num.to(dev)
    

    def get_sigma_n_2(self):
        """
        Returns the variance of the noise (measurement noise + rounding noise)
        """
        return torch.exp(self.sigma_n_log)**2 + self.sigma_n_num**2


    def forward(self, X):
        """
        Returns the prior distribution and the inverse anf the log_det of the prior covariace
        
        input:
        - X = training inputs (X has dimension [num_samples, num_features])
        
        output:
        - m_X = prior mean of X
        - K_X = prior covariance of X
        - K_X_inv = inverse of the prior covariance
        - log_det = log det of the prior covariance
        """
        # get the covariance
        N = X.size()[0]
        if self.GP_with_noise:
            K_X = self.get_covariance(X, flg_noise=True) 
        else:
            K_X = self.get_covariance(X)
        # get inverse and log_det with cholesky
        L = torch.linalg.cholesky(K_X)
        log_det = 2*torch.sum(torch.log(torch.diag(L)))
        K_X_inv = torch.cholesky_inverse(L)
        # get the mean
        m_X = self.get_mean(X)
        # return the values
        return m_X, K_X, K_X_inv, log_det


    def get_mean(self, X):
        """
        Returns the prior mean of X
        X is assumed to be a tenso of dimension [num_samples, num_features]
        """
        raise NotImplementedError()


    def get_covariance(self, X1, X2=None, flg_noise=False):
        """
        Returns the covariance betweeen the input locations X1 X2. 
        If X2 is None X2 is assumed to be equal to X1
        """
        raise NotImplementedError()


    def get_diag_covariance(self, X, flg_noise=False):
        """
        Returns the diagonal elements of the covariance betweeen the input locations X1
        """
        raise NotImplementedError()


    def get_alpha(self, X, Y):
        """
        Returns alpha, the vector of coefficients defining the posterior distribution

        inputs:
        - X = training input
        - Y = training output

        outputs:
        - alpha = vector defining the posterior distribution
        - m_X = prior mean of X
        - K_X_inv = inverse of the prior covariance of X
        """
        m_X, _, K_X_inv, _ = self.forward(X)
        alpha = torch.matmul(K_X_inv, Y-m_X)
        return alpha, m_X, K_X_inv


    def get_estimate_from_alpha(self, X, X_test, alpha, m_X, K_X_inv=None):
        """
        Compute the posterior distribution in X_test, given the alpha vector.
        
        input:
        - X = training input locations (used to compute alpha)
        - X_test = test input locations
        - alpha = vector of coefficients defining the posterior
        - m_X = prior mean of X
        - K_X_inv = inverse of the prior covariance of X

        output:
        - Y_hat = posterior mean
        - var = diagonal elements of the posterior variance
        
        If K_X_inv is given the method returns also the confidence intervals (variance of the gaussian)
        If Y_test is given the method prints the MSE
        """ 
        # get covariance and prior mean
        K_X_test_X = self.get_covariance(X_test, X)
        m_X_test = self.get_mean(X_test)
        # get the estimate 
        Y_hat = m_X_test + torch.matmul(K_X_test_X, alpha)
        # if K_X_inv is given compute the confidence intervals
        if K_X_inv is not None:
            num_test = X_test.size()[0]
            var = self.get_diag_covariance(X_test) - torch.sum(torch.matmul(K_X_test_X, K_X_inv)*(K_X_test_X), dim=1)
        else:
            var = None
        return Y_hat, var
         

    def get_estimate(self, X, Y, X_test, flg_return_K_X_inv=False):
        """
        Returns the posterior distribution in X_test, given the training samples X Y.
        
        input:
        - X = training input
        - Y = training output
        - X_test = test input

        output:
        - Y_hat = mean of the test posterior
        - var = diagonal elements of the variance posterior
        - alpha = coefficients defining the posterior
        - m_X = prior mean of the training samples
        - K_X_inv = inverse of the training covariance

        The function returns:
           -a vector containing the sigma squared confidence intervals
           -the vector of the coefficient
           -the K_X inverse in case required through flg_return_K_X_inv"""
        # get the coefficent and the mean
        alpha, m_X, K_X_inv = self.get_alpha(X, Y)
        # get the estimate and the confidence intervals
        Y_hat, var = self.get_estimate_from_alpha(X, X_test, alpha, m_X, K_X_inv=K_X_inv)
        #return the opportune values
        if flg_return_K_X_inv:
            return Y_hat, var, alpha, m_X, K_X_inv
        else:
            return Y_hat, var, alpha


    def print_model(self):
        """
        Print the model parameters
        """
        print(self.name+' parameters:')
        for par_name, par_value in self.named_parameters():
            print('-', par_name, ':', par_value.data)


    def fit_model(self,trainloader=None, 
                  optimizer=None, criterion=None,
                  N_epoch=1, N_epoch_print=1,
                  f_saving_model=None, f_print=None):
        """
        Optimize the model hyperparameters

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
                out_GP_priors = self(inputs)
                loss = criterion(out_GP_priors, labels)
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


    def get_SOD(self, X, Y, threshold, flg_permutation=False):
        """
        Returns the SOD points with an online procedure
        SOD: most importants subset of data
        """
        print('\nSelection of the inducing inputs...')
        # get number of samples
        num_samples = X.shape[0]
        # init the set of inducing inputs with the first sample
        SOD = X[0:1,:]
        inducing_inputs_indices = [0]
        # get a permuation of the inputs
        # perm_indices = torch.arange(1,num_samples)
        perm_indices = range(1,num_samples)
        if flg_permutation:
            perm_indices = perm_indices[torch.randperm(num_samples-1)]
        # iterate all the samples
        for sample_index in perm_indices:
            # get the estimate 
            _, var, _ = self.get_estimate(X[inducing_inputs_indices,:], Y[inducing_inputs_indices,:], X[sample_index:sample_index+1,:])
            if torch.sqrt(var)>=threshold:
                SOD = torch.cat([SOD,X[sample_index:sample_index+1,:]],0)
                inducing_inputs_indices.append(sample_index)
            # else:
            #     print('torch.sqrt(var) = ',torch.sqrt(var))
            #     print('threshold = ',threshold)
        print('Shape inducing inputs selected:', SOD.shape)
        return inducing_inputs_indices


    def __add__(self, other_GP_prior):
        """
        Returns a new GP given by the sum of two GP
        """
        return Sum_Independent_GP(self, other_GP_prior)


    def __mul__(self, other_GP_prior):
        """
        Returns a new GP with mean given by the product of the means and covariances
        """
        return Multiply_GP_prior(self, other_GP_prior)




class Combine_GP(GP_prior):
    """
    Class that extend GP_prior and provide common utilities to combine GP
    """


    def __init__(self, *gp_priors_obj):
        """
        Initialize the multiple kernel object
        """
        # initialize a new GP object
        super(Combine_GP, self).__init__(active_dims=[], sigma_n_num=gp_priors_obj[0].sigma_n_num,
                                         scale_init=np.ones(1), flg_train_scale=False,
                                         dtype=gp_priors_obj[0].dtype, device=gp_priors_obj[0].device)
        # build a list with all the models
        self.gp_list = torch.nn.ModuleList(gp_priors_obj)
        # check the noise flag
        GP_with_noise = False
        for gp in self.gp_list:
            GP_with_noise = GP_with_noise or gp.GP_with_noise 
        self.GP_with_noise = GP_with_noise


    def to(self, dev):
        """
        Move all the models to the desired device
        """
        super(Combine_GP, self).to(dev)
        self.device = dev
        for gp in self.gp_list:
            gp.to(dev)


    def print_model(self):
        """
        Print the parameters of all the models in the gp_list
        """
        for gp in self.gp_list:
            gp.print_model()


    def get_sigma_n_2(self):
        """
        Iterate over all the models in the list and returns the global noise variance
        """
        sigma_n_2 = torch.zeros(1, dtype=self.dtype, device=self.device)
        for gp in self.gp_list:
            if gp.GP_with_noise:
                sigma_n_2 += gp.get_sigma_n_2()
        return sigma_n_2




class Sum_Independent_GP(Combine_GP):
    """
    Class that sum GP_priors objects
    """


    def __init__(self, *gp_priors_obj):
        """
        Initialize the gp list
        """
        super(Sum_Independent_GP, self).__init__(*gp_priors_obj)


    def get_mean(self, X):
        """
        Returns the sum of the means returned by the models in gp_list
        """
        N = X.size()[0]
        mean = torch.zeros(N,1, dtype=self.dtype, device=self.device)
        for gp in self.gp_list:
            mean += gp.get_mean(X)
            return mean


    def get_covariance(self, X1, X2=None, flg_noise=False):
        """
        Returns the sum of the covariances returned by the mdoels in gp_list
        """
        # get dimensions
        N1 = X1.size()[0]
        if X2 is None:
            N2 = N1
        else:
            N2 = X2.size()[0]
        # initialize the covariance
        cov = torch.sum(torch.cat([(gp.get_covariance(X1,X2,flg_noise=False)).unsqueeze(0) for gp in self.gp_list],0),0)
        # add the noise
        if flg_noise & self.GP_with_noise:
            cov += self.get_sigma_n_2()*torch.eye(N1, dtype=self.dtype, device=self.device)
        return cov


    def get_diag_covariance(self, X, flg_noise=False):
        """
        Returns the sum of the diagonals of the covariances in gp list
        """
        # initialize the vector
        diag = torch.zeros(X.size()[0], dtype=self.dtype, device=self.device)
        # iterate in the list and sum the diagonals
        for gp in self.gp_list:
            diag += gp.get_diag_covariance(X, flg_noise=False)
        # add the noise
        if flg_noise & self.GP_with_noise:
            diag += self.get_sigma_n_2()
        return diag




class Multiply_GP_prior(Combine_GP):
    """
    Class that generates a GP_prior multiplying GP_priors objects
    """


    def __init__(self, *gp_priors_obj):
        """
        Initilize the GP list
        """
        super(Multiply_GP_prior, self).__init__(*gp_priors_obj)


    def get_mean(self, X):
        """
        Returns the product of the means returned by the models in gp_list
        """
        # initilize the mean vector
        N = X.size()[0]
        mean = torch.ones(N,1, dtype=self.dtype, device=self.device)
        # multiply all the means
        for gp in self.gp_list:
            mean = mean*gp.get_mean(X)
        return mean


    def get_covariance(self, X1, X2=None, flg_noise=False):
        """
        Returns the element-wise product of the covariances returned by the models in gp_list
        """
        # get size
        N1 = X1.size()[0]
        if X2 is None:
            N2 = N1
        else:
            N2 = X2.size()[0]
        cov = torch.prod(torch.cat([(gp.get_covariance(X1,X2,flg_noise=False)).unsqueeze(0) for gp in self.gp_list],0),0)
        if flg_noise & self.GP_with_noise:
            cov = cov + self.get_sigma_n_2()*torch.eye(N1, dtype=self.dtype, device=self.device)
        return cov


    def get_diag_covariance(self, X, flg_noise=False):
        """
        Returns the product of the diagonal elements of the covariance returned by the models in gp_list
        """
        # initilize the diagona
        N = X.size()[0]
        diag = torch.ones(N, dtype=self.dtype, device=self.device)
        # multiply all the diagonals
        for gp in self.gp_list:
            diag *= gp.get_diag_covariance(X, flg_noise=False)
        # add the nosie
        if flg_noise & self.GP_with_noise:
           diag +=self.get_sigma_n_2()
        return diag




def Scale_GP_prior(GP_prior_class, GP_prior_par_dict,
                   f_scale, active_dims_f_scale,
                   pos_par_f_init=None, flg_train_pos_par_f=True,
                   free_par_f_init=None, flg_train_free_par_f=True,
                   additional_par_f_list=[]):
    """
    Funciton that returns a GP_prior scaled. This class implement the following model:
    y(x) = a(x)f(x) + e, where f(x) is a GP and a(x) a deterministic function.
    The function a(x) can be parametrize with respect to a set of trainable prameters.
    This class retuns an instance of a new class defined dynamically in the following
    """

    # define dynamically the new class
    class Scaled_GP(GP_prior_class):
        """
        Class that extends the GP_prior_class with the scaling parameters
        """

        def __init__(self, GP_prior_par_dict,
                     f_scale, active_dims_f_scale,
                     pos_par_f_init, flg_train_pos_par_f,
                     free_par_f_init, flg_train_free_par_f,
                     additional_par_f_list):
            # initialize the object of the superclass
            super(Scaled_GP, self).__init__(**GP_prior_par_dict)
            # save the scaling info 
            self.f_scale = f_scale
            self.active_dims_f_scale = active_dims_f_scale
            self.additional_par_f_list = additional_par_f_list
            if pos_par_f_init is None:
                self.flg_pos_par = False
                self.pos_par_f_log = None
            else:
                self.flg_pos_par = True
                self.pos_par_f_log = torch.nn.Parameter(torch.tensor(np.log(pos_par_f_init), dtype=self.dtype, device=self.device), requires_grad=flg_train_pos_par_f)
            if free_par_f_init is None:
                self.flg_free_par = False
                self.free_par_f = None            
            else:
                self.flg_free_par = True
                self.free_par_f = torch.nn.Parameter(torch.tensor(free_par_f_init, dtype=self.dtype, device=self.device), requires_grad=flg_train_free_par_f)
        

        def get_scaling(self,X):
            """
            Returns the scaling funciton evaluated in X
            """
            if self.flg_pos_par:
                pos_par = torch.exp(self.pos_par_f_log)
            else:
                pos_par = None
            return self.f_scale(X[:, self.active_dims_f_scale], pos_par, self.free_par_f, *self.additional_par_f_list).reshape(-1)


        def get_mean(self, X):
            """
            Calls the get_mean of the superclass and apply the scaling
            """
            # get the supercalss mean and scale the result
            return self.get_scaling(X)*super(Scaled_GP, self).get_mean(X)


        def get_covariance(self, X1, X2=None, flg_noise=False):
            """
            Calls the get_covariance of the superclass and apply the scaling
            """
            # get the scaling functions
            a_X1 = self.get_scaling(X1).reshape([-1,1])
            # if required evaluate the scaling function in X2 and get the covariance
            if X2 is None:
                K = a_X1*super(Scaled_GP, self).get_covariance(X1, X2, flg_noise=False)*(a_X1.transpose(0,1))
            else:
                a_X2 = self.get_scaling(X2).reshape([-1,1])
                K = a_X1*super(Scaled_GP, self).get_covariance(X1, X2, flg_noise=False)*(a_X2.transpose(0,1))
            # if required add the noise and return the covariance
            if flg_noise & self.GP_with_noise:
                N = K.size()[0]
                return K + self.get_sigma_n_2()*torch.eye(N, dtype=self.dtype, device=self.device)
            else:
                return K 


        def get_diag_covariance(self, X, flg_noise=False):
            """
            Calls the get_diag_covariance of the superclass and apply the scaling
            """
            # evaluate the scaling function in X1
            a_X = self.get_scaling(X)
            diag = a_X**2*super(Scaled_GP, self).get_diag_covariance(X, flg_noise=False)
            # if required add the noise and return the covariance
            if flg_noise & self.GP_with_noise:
                return diag + self.get_sigma_n_2()
            else:
                return diag


    # return an object of the new class
    return Scaled_GP(GP_prior_par_dict,
                     f_scale, active_dims_f_scale,
                     pos_par_f_init, flg_train_pos_par_f,
                     free_par_f_init, flg_train_free_par_f,
                     additional_par_f_list)