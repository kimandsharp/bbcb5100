#!/usr/bin/env python
"""
implement bayesian estimation of mean of population, using exact (t-distribution) 
and chi-sq posterior pdf of std. dev
20 may 2023- add option to do log normal: replace x by ln x
"""
from math import sqrt, exp,log
import numpy as np
import matplotlib.pyplot as plt
import sys
from SciInf_utilities import *
#import pymol_cgo
#--------------------------------------

print("\nImplement bayesian estimation of location parameter with Gaussian noise")
print("using exact (t-distribution) and chi-sq posterior pdf of magnitude of noise")
print('option to model distribution as log normal by adding 2nd argument "log" after filename\n')
# main
x = []
y = []
logNormal = False
if(len(sys.argv)>1):
  file1 = sys.argv[1]
else:
  file1 = input('data file with one value per line> ')
if(len(sys.argv)>2):
  if(sys.argv[2][0:3] == 'log'): logNormal = True
print('lognormal ',logNormal)
n_x = read_x(x,file1)
if(logNormal):
  if(min(x) > 0.):
    for i in range(n_x):
      x[i] = log(x[i])
  else:
    print('not all data is > 0, cannot use log normal distribution')
    logNormal = false
#
# basic stats
#
min_x = min(x)
max_x = max(x)
av_x = average_x(x)
av_xx = average_xy(x,x)
var_x = av_xx - av_x**2
sigma_x = sqrt(var_x)
sigma_av = sqrt(var_x/n_x)
for i in range(len(x)):
  y.append(0.5)
#
print('\n===========================================================')
print('sample (data) summary')
print('===========================================================')
print(' Min X {:12.5f} Max X    {:12.5f} '.format(min_x,max_x))
print(' Av X  {:12.5f} Var of X {:12.5f} Std.dev of  X  {:12.5f}'.format(av_x,var_x,sigma_x))
print('===========================================================\n')
exponent = n_x/2. # result if use log prior for sigma or prob(sigma) = const./sigma
#exponent = (n_x-1)/2. # result if use flat prior for sigma or prob(sigma) = const
#
#generate posterior pdf and cdf for mean
#
xrange = 4. # sigma range for x-axis
av_min = av_x - xrange*sigma_av
av_incr = 2*xrange*sigma_av/(NPOINT - 1)
av_axis = np.zeros(NPOINT)
log_av_pdf = np.zeros(NPOINT)
for i in range(NPOINT):
  av_axis[i] = av_min + i*av_incr
  log_av_pdf[i] = -1.*exponent*log(1. + (av_axis[i] - av_x)**2/var_x)
av_pdf = np.exp(log_av_pdf)
#
if(logNormal): # transform back
  for i in range(n_x):
    x[i] = exp(x[i])
  for i in range(NPOINT):
    av_axis[i] = exp(av_axis[i])
    av_pdf[i] = av_pdf[i]/av_axis[i] # *  |dy/dx| = dlnx/dx = 1/x
pdf_max = max(av_pdf)
av_pdf = av_pdf/pdf_max
#
av_cdf = pdf_to_cdf(av_axis,av_pdf)
write_pdf_cdf(av_axis,av_pdf,av_cdf,title='Location Parameter pdf cdf',filename='location_pdf_cdf.dat')
#
summarize(av_axis,av_pdf,av_cdf,title='Location Parameter')
#
# plot original data
#
if(MAKEPLOT):
  plt.figure(1)
  plt.subplot(211)
  plt.boxplot(x,notch=0,sym='b+',vert=0,showmeans=True)
  plt.yticks([1],['X 1'],rotation=0,fontsize=12)
  plt.title('Input Data')
  #
  # plot posterior pdf, cdf for mean
  #
  plt.subplot(212)
  plt.plot(av_axis,av_pdf,'g-')
  plt.plot(av_axis,av_cdf,'r-')
  plt.scatter(x,y)
  plt.title('posterior pdf,cdf for Location Parameter X')
  plt.xlabel('X')
  plt.ylabel('p(X)')
  plt.ylim((0.,1.2))
  plt.grid(True)
  plt.show()
#
#generate posterior pdf and cdf for st.dev
#
xrange = 4. # range for x-axis
sd_min = sigma_x/xrange
sd_max = sigma_x*xrange
sd_incr = (sd_max - sd_min)/(NPOINT - 1)
sd_axis = np.zeros(NPOINT)
log_sd_pdf = np.zeros(NPOINT)
for i in range(NPOINT):
  sd_i = sd_min + i*sd_incr
  var_i = sd_i*sd_i
  sd_axis[i] = sd_i
  #sd_pdf[i] = exp(-0.5*n_x*var_x/var_i)/sd_i**n_x
  log_sd_pdf[i] = (-0.5*n_x*var_x/var_i) - n_x*log(sd_i)
pdf_max = max(log_sd_pdf)
log_sd_pdf = log_sd_pdf - pdf_max
sd_pdf = np.exp(log_sd_pdf)
sd_cdf = pdf_to_cdf(sd_axis,sd_pdf)
write_pdf_cdf(sd_axis,sd_pdf,sd_cdf,title='Noise Magnitude pdf cdf',filename='noise_pdf_cdf.dat')
#
summarize(sd_axis,sd_pdf,sd_cdf,title='Noise Magnitude')
#
# plot posterior pdf, cdf of st. dev
#
if(MAKEPLOT):
  plt.figure(1)
  plt.plot(sd_axis,sd_pdf,'g-')
  plt.plot(sd_axis,sd_cdf,'r-')
  plt.title('posterior pdf,cdf for Noise magnitude Z')
  plt.xlabel('Z')
  plt.ylabel('p(Z)')
  plt.grid(True)
  plt.show()
"""
#
# output joint p(mean, stdev) to file for plotting
#
print(av_axis)
print(sd_axis)
fileout = open('meanStd.dat','w')
fileout.write('# data for 3d plot of log p(mean,stdev) from MeanStd.py \n')
ilw = int(NPOINT/2 - 10)
iup = int(NPOINT/2 + 10)
#for i in range(ilw,iup):
for i in range(0,NPOINT,10):
  av_i = av_axis[i]
  #for j in range(ilw,iup):
  for j in range(0,NPOINT,10):
    sd_i = sd_axis[j]
    logProb = -(n_x + 1)*log(sd_i) - n_x*(var_x + (av_i - av_x)**2)/2/sd_i**2
    fileout.write('{:8.3f}  {:8.3f}  {:12.3g}\n'.format(av_i,sd_i,logProb))
fileout.close
"""
