#!/usr/bin/env python
"""
implement  standard t-test
"""
from math import sqrt, exp, log
import numpy as np
import matplotlib.pyplot as plt
from SciInf_utilities import *
import random
import sys
#--------------------------------------

print("\n t-test to see if difference in 2 population means")
print("is significant\n")
x = []
y = []
"""
#
# dummy up some random data
#
random.seed(1234567)
S1 = 5.
A1 = 0.44
for i in range(20):
  xi = S1*random.random() + A1
  yi = S1*random.random() 
  x.append(xi)
  y.append(yi)
  print('%9.4f  %9.4f '%(xi,yi))
n_x = len(x)
n_y = len(y)
"""
if(len(sys.argv) > 2):
  file1 = sys.argv[1]
  file2 = sys.argv[2]
else:
  file1 = input('first data file with one value per line> ')
  file2 = input('second data file with one value per line> ')
print('\n input file 1: ',file1)
print(' input file 2: ',file2,'\n')
n_x = read_x(x,file1)
n_y = read_x(y,file2)
"""
x = [30.02,29.99,30.11,29.97,30.01,29.99]
n_x = len(x)
y = [29.89,29.93,29.72,29.98,30.02,29.98]
n_y = len(y)
"""
dof = n_x + n_y - 2
av_x = average_x(x)
av_y = average_x(y)
av_xx = average_xy(x,x)
av_yy = average_xy(y,y)
var_x = av_xx - av_x**2
var_y = av_yy - av_y**2
sd_x = sqrt(var_x)
sd_y = sqrt(var_y)
nterm = sqrt(1./n_x + 1./n_y)
print('dof, nterm: ',dof,nterm)
sigma_xy = sqrt((n_x*var_x + n_y*var_y)/dof)
dav = av_y - av_x
tStat = dav/(nterm*sigma_xy)
print('\n===========================================================')
print('sample (data) summary')
print('===========================================================')
print(' Av X1 {:12.5f} Av X2 {:12.5f} Var of X1 {:12.5f} Var of X2 {:12.5f} '.format(av_x,av_y,var_x,var_y))
print(' Av X2 - Av X1 {:12.5f} '.format(dav))
print(' sigma <X2-X1> {:12.5f} '.format(sigma_xy))
print(' tStat {:12.5f} '.format(tStat))
print('===========================================================\n')
#
#
# generate t-distribution
#
t_min = -5.
t_max = +5.
t_incr = (t_max - t_min)/(NPOINT - 1)
t_axis = np.zeros(NPOINT)
t_pdflog = np.zeros(NPOINT)
for i in range(NPOINT):
  t_axis[i] = t_min + i*t_incr
  t_pdflog[i] = -0.5*(dof + 1.)*log(1. + t_axis[i]**2/dof)
pdf_max = max(t_pdflog)
t_pdflog = t_pdflog - pdf_max
t_pdf = np.exp(t_pdflog)
t_cdf = pdf_to_cdf(t_axis,t_pdf)
#write_pdf_cdf(t_axis,t_pdf,t_cdf,title='t-dist pdf cdf',filename='tDist_pdf_cdf.dat')
p05 = quantile(t_axis,t_cdf,2.5)
print('t value for 2 tailed p=0.05: +',p05)
i = 0
while(i<NPOINT):
 if(t_axis[i] < tStat):
   i += 1
 else:
   break
print('1-tailed probability: %9.3f'%(t_cdf[i]))
#sys.exit()
#
# plot original data
#
data_all = [x, y]
MAKEPLOT = True
if(MAKEPLOT):
  plt.figure(1)
  plt.subplot(211)
  plt.boxplot(data_all,notch=0,sym='b+',vert=0,showmeans=True)
  plt.yticks([1,2],['X 1','X 2'],rotation=0,fontsize=12)
  plt.title('input data')
  #plt.xlim(xmin=0.3,xmax=1.)
  #plt.show()
  plt.subplot(212)
  plt.plot(t_axis,t_pdf,'g-')
  plt.plot(t_axis,t_cdf,'r-')
  #plt.title('posterior pdf,cdf for diff. in means')
  plt.ylim((0.,1.2))
  plt.xlabel('t')
  plt.ylabel('p(t)')
  plt.grid(True)
  plt.show()
