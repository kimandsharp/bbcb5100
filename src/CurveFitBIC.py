#!/usr/bin/env python
"""
python version of curve_fit_bayes.f
using Bayesian information criterion (BIC)
G. Schwartz "Estimating the dimension of a model"
(1978) Ann. Stats. 6:461-464
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
from numpy.linalg import inv
import sys
MAKEPLOT = True
#-------------------------------
# model is linear combination of basis functions-
# 
def base_func(ic,xx):
  # ic is index into which basis function
  # here  model is a polynomial, so basis functions are 
  # simply x^i, so index i just specifies exponent
  fx  = xx**ic
  return fx
#-------------------------------
ncmax = 8 # max number of basis functions (1=constant, 2=linear, 3=quadratic, etc)
print('\nCurve fit to a set of basis functions')
print('using Bayesian information criterion (BIC)')
print('G. Schwartz "Estimating the dimension of a model"')
print('(1978) Ann. Stats. 6:461-464')
print('BIC = log Likelihood - 0.5 k log(n)')
print('n is number of (x,y) data pairs, k is number of parameters ')
print('since basis functions are powers of x  1=constant, 2=linear, 3=quadratic, etc\n')
#
# read in data
#
if(len(sys.argv)>1):
  input_file = sys.argv[1]
else:
   input_file = input("file with one x, y data pair per line> ")
   #input_file = 'xy.dat'
print('input file: ',input_file)
xraw = []
yraw = []
ndata = read_xy(xraw,yraw,input_file)
xx = np.zeros(ndata,'float')
yy = np.zeros(ndata,'float')
print('--------------')
print(' Input Data')
print('--------------')
for i in range(ndata):
  xx[i] = xraw[i]
  yy[i] = yraw[i]
  print('%12.5g %12.5g ' %(xx[i],yy[i]))
#
# basics stats
#
yylw = np.min(yy)
yyup = np.max(yy)
xx_mean = np.mean(xx)
yy_mean = np.mean(yy)
x2 = xx*xx
y2 = yy*yy
xx_var = np.average(x2) - xx_mean**2
yy_var = np.average(y2) - yy_mean**2
print('mean of     x: %12.5g y: %12.5g ' % (xx_mean,yy_mean))
print('variance of x: %12.5g y: %12.5g ' % (xx_var,yy_var))
print('y  min: %12.5g max: %12.5g ' %(yylw,yyup))
#
# 1st term in exponent of gaussian posterior
#
V0 = np.sum(y2)
yyfit = np.zeros((ndata))
dyy = np.zeros(ndata)
#
# loop over increasing # of basis functions/polynomial order
#
ncup = min(ndata-3,ncmax)
nc_best = 1
logBic_best = -100.
done = False
for n in range(0,ncup):
  nc = n + 1
  print('---------------')
  print('# of parameters: ',nc)
  print('---------------')
  Bvec       = np.zeros(nc)
  Cmin       = np.zeros(nc)
  Amat       = np.zeros((nc,nc))
  Amatinv    = np.zeros((nc,nc))
  Covmat     = np.zeros((nc,nc))
  Csigma     = np.zeros(nc)
  #
  # solve for coefficients that minimize ss error
  # ( max likelihood solution
  #
  for i in range(ndata):
    for j in range(nc):
      Bvec[j] += yy[i]*base_func(j,xx[i])
      for l in range(nc):
        Amat[j][l] += base_func(j,xx[i])*base_func(l,xx[i])
  Amatinv = inv(Amat)
  for j in range(nc):
    for l in range(nc):
      Cmin[j] += Amatinv[j][l]*Bvec[l]
  Vmin = 0.
  #
  # residuals and sum sq. err
  #
  for i in range(ndata):
    yyfit[i] = 0.
    for j in range(nc):
      yyfit[i] += Cmin[j]*base_func(j,xx[i])
    dyy[i] = yyfit[i] - yy[i]
    Vmin += dyy[i]**2
  print('Sum Sq. dev: %15.6g '% (Vmin))
  Vmin = V0
  for j in range(nc):
    Vmin -= Bvec[j]*Cmin[j] # alternative way to find SS dev
  #
  # curvature matrix gives sigmas for coefficients
  # noise estimated from rms residuals, or 'errors' see Gull, 1988
  sigma2 = Vmin/(ndata - nc)
  sigma = math.sqrt(abs(sigma2))
  for j in range(nc):
    for l  in range(nc):
      Covmat[j][l] = Amatinv[j][l]*sigma2
  for j in range(nc):
    if(Covmat[j][j] < 0.):
      print('Higher order polynomial fit degenerate - stopping here\n')
      done = True
    if(done): break
    Csigma[j]  = math.sqrt(Covmat[j][j])
  if(done): break
  #
  # posterior values for coefficients- includes shrinkage from prior
  # for coefficients
  # and posterior for # of terms nc
  #
  print('coefficients    +/- 1 sigma')
  for j in range(nc):
    print('%12.5g %12.5g  '%(Cmin[j],Csigma[j]))
  print('|noise| %12.5g     Sum. Sq. dev. fit: %12.5g ' % (sigma,Vmin))
  logBic = - math.log(Vmin) - 0.5*nc*math.log(ndata)
  #logBic = - Vmin - 0.5*nc*math.log(ndata) # BIC alternative 1
  #logBic = - math.log(Vmin) - nc # BIC alternative 2
  if(logBic > logBic_best):
    logBic_best = logBic
    nc_best = nc
  print('Bayesian Info. Criterion log P for # parameters = %1d: %12.5g ' % (nc,logBic))
  
  #
  #-----------------------------------------
  if(MAKEPLOT):
    if(nc == 1):
      ptitle = 'fit for %3d parameter' % (nc)
    else:
      ptitle = 'fit for %3d parameters' % (nc)
    plt.figure()
    plt.subplot(211)
    plt.scatter(xx,yy,color='red',marker='o')
    plt.plot(xx,yyfit)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(ptitle)
    plt.grid(True)
    #
    plt.subplot(212)
    plt.scatter(xx,dyy,color='red',marker='o')
    plt.title('residuals')
    plt.xlabel('x')
    plt.ylabel('dy')
    plt.grid(True)
    plt.show()
  #-----------------------------------------
print('=================')
print('Best # of parameters = %1d with Bayesian Info. Criterion log P %12.5g ' % (nc_best,logBic_best))
print('=================')

sys.exit()
