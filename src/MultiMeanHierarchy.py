#!/usr/bin/env python
"""
Bayesian analysis of multiple populations characterized by J (mean, sample standard 
error mean_j, sigma_j) pairs, or alternatively J sets {yi}j of N_j size raw data sets
using a hierarchical gaussian model characterized by hyper-parameters 
mu, tau govering gaussian prior distribution of population means
using the approach of gelman et al, DBA3 chapter 5, e.g. the famous 8-schools case
"""
import random as rn
import numpy as np
import matplotlib.pyplot as plt
import SciInf_utilities as ku
from math import *
percentsign = '%'
#-------------------------------------
def mu_params(tau,means_data,sterr_data):
  """ calculates  global mean, spread for know tau
    DBA3 eq. 5.20"""
  nset = len(means_data)
  mu_prec = 0.
  mu_av = 0.
  for j in range(nset):
    mu_precj = 1./(sterr_data[j]**2 + tau**2)
    mu_prec += mu_precj
    mu_av += mu_precj*means_data[j]
  mu_av /= mu_prec
  return mu_av, mu_prec

def tau_post(tau,means_data,sterr_data):
  """ calculate posterior marginal for hyperparamter tau conditioned
  on data set means, std. errs according to DBA3 eqs. 5.20, 5.21 """
  nset = len(means_data)
  mu_av, mu_prec = mu_params(tau,means_data,sterr_data)
  mu_stdev = sqrt(1./mu_prec)
  logprob = 0.
  for j in range(nset):
    mu_precj = 1./(sterr_data[j]**2 + tau**2)
    logprob += 0.5*log(mu_precj) - 0.5*(means_data[j] - mu_av)**2*mu_precj
  tau_prob = mu_stdev*exp(logprob)
  return tau_prob
#
#-------------------------------------
print('\nBayesian hierarchical multiple means analysis')
print('using a hierarchical model characterized by hyper-parameters ')
print('mu, tau govering gaussian prior distribution of population means')
print('using the approach of Gelman et al, DBA3 chapter 5, in the 8-schools case\n')
print('\nInput data as \n(1) pairs of mean, std. error \n(2) J sets of raw data {y_i}')
input_type = int(input('\nOption 1 or 2 >> '))
#input_type = 1 # debug
means_data = []
sterr_data = []
if(input_type == 1):
  data_file = input('\nname of file containing one (mean, std. err) per line>> ')
  #data_file = 'EightSchools.dat' # debug
  print('\nreading mean, std. err data from file ',data_file)
  ku.read_xy(means_data,sterr_data,data_file)
  nset = len(means_data)
else:
  data_raw = []
  while(1):
    data_file = input('\nname of file  1 set of raw data, one float per line (exit to end input)>> ')
    if(data_file == 'exit'): break
    yj = []
    ku.read_x(yj,data_file)
    data_raw.append(yj)
  # done with file reads,  convert raw data to mean, std errors
  nset = len(data_raw)
  print('\n# of data sets read: ',nset )
  for j in range(nset):
    muj = ku.average_x(data_raw[j])
    nj = len(data_raw[j])
    mu2j = ku.average_xy(data_raw[j],data_raw[j])
    varj = mu2j - muj*muj
    sterrj = sqrt(varj/(nj - 1))
    means_data.append(muj)
    sterr_data.append(sterrj)
#
print('\ninput data mean    stderr: ')
for j in range(nset):
  print('%12.5f %12.5f ' % (means_data[j],sterr_data[j]))
#
#==============================================================
# compute mean of means, and stdev of means to get overall location, and
# scale range for hyperparameter tau
#==============================================================
mu_mu = ku.average_x(means_data)
mu2 = ku.average_xy(means_data,means_data)
mu_var = mu2 - mu_mu*mu_mu
mu_stdev = sqrt(mu_var)
print('\nglobal mean, stdev of means: %12.5f %12.5f '% (mu_mu,mu_stdev))
tau_range = 8.*mu_stdev
#
#==============================================================
# compute posterior marginal distbn. for hyperparameter mu, tau, the center, spread of means
# assuming a uniform prior for tau and mu
#==============================================================
#NPOINT = 2501 # debug
NPOINT = ku.NPOINT
dtau = tau_range/(NPOINT - 1.)
tau_axis = np.zeros(NPOINT)
tau_prob = np.zeros(NPOINT)
tau_val = 0.
for i in range(NPOINT):
  tau_axis[i] = tau_val
  tau_prob[i] = tau_post(tau_val,means_data,sterr_data)
  tau_val += dtau
tau_prob /= np.max(tau_prob)
tau_cdf = ku.pdf_to_cdf(tau_axis,tau_prob)
ndim = 99
quantile_axis = np.zeros(ndim) # values for inverse cdf of p(tau) for sampling
tau_quantile = np.zeros(ndim) # value of tau for each percentile
for i in range(ndim):
  quantile_axis[i] = 1.*(i+1)
  tau_quantile[i] = ku.quantile(tau_axis,tau_cdf,quantile_axis[i])
tau_up = tau_quantile[95]  
print('tau 95{:1s} limits: ({:12.5f} to {:12.5f})'.format(percentsign,0.,tau_up))
#print(tau_axis)
#print(tau_quantile)

if(ku.MAKEPLOT):
  plt.figure(1)
  plt.title('posterior marginal for hyperparameter tau (spread of means)')
  plt.plot(tau_axis,tau_prob,'g-')
  plt.plot(tau_axis,tau_cdf,'r-')
  plt.xlim(0.,tau_axis[-1])
  plt.ylim(0.,1.1)
  plt.xlabel('tau ')
  plt.ylabel('p(tau|data) ')
  plt.show()
#
"""
if(ku.MAKEPLOT):
  plt.figure(2)
  plt.title(' inverse cdf for p(tau|data')
  plt.plot(quantile_axis,tau_quantile,'g-')
  plt.xlim(0.,100.)
  plt.xlabel('%')
  plt.ylabel('tau ')
  plt.show()
"""
#
# sample global spread paramter tau from p(tau|data)
# then sample global location mu from p(mu|tau,data) = Normal(mu_mu,mu_stdev)
# (DBA3 eq. between 5.19 and 5.20)
# then for each data set, sample its location parameter theta_j 
# from p(theta_j | mu, tau, data) = Normal(theta_av_j, stdev_j)
# (DBA3 eq. 5.17)
rn.seed(123)
tau_sample = []
mu_check = 0.
nsample = 5000
theta_sample = np.zeros((nset,nsample))
#ts = np.zeros((nsample,nset))
print('\nsampling randomly from posterior...')
for i in range(nsample):
  i1 = rn.randint(0,ndim-1)
  tau_val = tau_quantile[i1]
  tau_sample.append(tau_val)
  mu_av, mu_prec = mu_params(tau_val,means_data,sterr_data)
  mu_stdev = sqrt(1./mu_prec)
  mu_val = rn.normalvariate(mu_av,mu_stdev)
  mu_check += mu_val
  #print(i1,tau_val,mu_val)
  for j in range(nset):
    theta_prec = 1./sterr_data[j]**2 + 1./tau_val**2
    theta_stdev = sqrt(1./theta_prec)
    theta_av = (means_data[j]/sterr_data[j]**2 + mu_val/tau_val**2)/theta_prec
    theta_val = rn.normalvariate(theta_av,theta_stdev)
    theta_sample[j][i] = theta_val
    #ts[i][j] = theta_val
mu_check /= nsample
print('mean of means from %d samples: %12.5f   ' % (nsample,mu_check))
#
# extract median, 95% credible interval ranges
# and also produce 'forest plot', basically a very skinny boxplot
# good when there are a lot of sets
yarr = np.zeros(nset)
xq25 = np.zeros((2,nset))
xq95 = np.zeros((2,nset))
xmed = np.zeros(nset)
# indices for various quantiles in sorted sample
q50 = nsample//2
q025 = int(nsample*0.025)
q975 = int(nsample*0.975)
q25 = int(nsample*0.25)
q75 = int(nsample*0.75)
for j in range(nset):
  theta_sample[j].sort()
  yarr[j] = j + 1
  xmed[j] = theta_sample[j][q50]
  #xq25[j] = (theta_sample[j][q75] - theta_sample[j][q25])/2.
  #xq95[j] = (theta_sample[j][q975] - theta_sample[j][q025])/2.
  xq25[0][j] = (xmed[j] - theta_sample[j][q25])
  xq25[1][j] = (theta_sample[j][q75] - xmed[j])
  xq95[0][j] = (xmed[j] - theta_sample[j][q025])
  xq95[1][j] = (theta_sample[j][q975] - xmed[j])
  print('set %3d median %12.5f 95%s CI (%12.5f , %12.5f) ' \
  %(j+1,theta_sample[j][q50],percentsign,theta_sample[j][q025],theta_sample[j][q975]))
#print(theta_sample)
if(ku.MAKEPLOT):
  plt.figure(3)
  plt.errorbar(xmed,yarr,xerr=xq95,elinewidth=1,fmt='.',ecolor='r')
  plt.errorbar(xmed,yarr,xerr=xq25,elinewidth=2,fmt='.',ecolor='b')
  plt.title('posterior credible intervals for mean (blue=50% red=95%)')
  plt.xlabel('Mean')
  plt.ylabel('Data Set')
  plt.show()
