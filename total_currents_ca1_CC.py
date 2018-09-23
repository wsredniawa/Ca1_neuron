# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 09:04:55 2018

@author: Wladek
"""

import numpy as np
import os
import pylab as py
from mpl_toolkits.mplot3d import axes3d
import matplotlib.cm as cm

py.close('all')

os.chdir('\\Users\\Wladek\\Dysk Google\\Ca1neuron\\neuron_model\\arrayTomography')

posxyz = np.loadtxt('matRecxyz.txt', skiprows = 1).T
curr = np.loadtxt('matRecData.txt').T
posswc = np.loadtxt('twinApical.swc').T
time = np.loadtxt('Rectime.txt')
inj = np.loadtxt('current_inj.txt', skiprows = 2)

ntime = len(time)

tp_man = 100 #time point to be plotted in 3d
skip = 4 #downsampling of .swc morphology

vs = np.reshape(curr[-1], (np.shape(posxyz)[1], ntime))
ks = np.reshape(curr[-2], (np.shape(posxyz)[1], ntime))
nas = np.reshape(curr[-3], (np.shape(posxyz)[1], ntime))
pas = np.reshape(curr[-4], (np.shape(posxyz)[1], ntime))
cap = np.reshape(curr[0], (np.shape(posxyz)[1], ntime))

for ii in range(ntime):
    ks[:,ii] *= posxyz[3]
    nas[:,ii] *=posxyz[3]
    pas[:,ii] *=posxyz[3]
    cap[:,ii] *=posxyz[3]

posxyz[2]*=1e1
posswc[4]*=1e1

wsp_plot = posxyz
xmin = np.min(wsp_plot[0])#-70
xmax = np.max(wsp_plot[0])#-45
ymin = np.min(wsp_plot[1])#-70
ymax = np.max(wsp_plot[1])#0
zmin = int(np.min(wsp_plot[2]))#-10
zmax = int(np.max(wsp_plot[2]))#60

n_ele = 32
wsp_ele=np.zeros((3,n_ele))
wsp_ele[1]=np.linspace(ymin, ymax, n_ele)
wsp_ele[0] = 100
wsp_ele[2] = 0
pots = np.zeros((ntime,n_ele))

for ele in range(n_ele):
    for n in range(np.shape(posxyz)[1]):
        dist = np.sqrt((wsp_ele[0, ele]- posxyz[0, n])**2+(wsp_ele[1, ele]- posxyz[1, n])**2+(wsp_ele[2, ele]- posxyz[2, n])**2)
        pots[:, ele] += (-ks[n, :]-nas[n, :]-pas[n,:])/dist
sigma= 0.3
pots /= 4*np.pi*sigma

print('injected current:', np.sum(inj[:,1]))
print('sum of all currents all compartments', np.sum(ks+nas+pas+cap)*1e-2)
#%%
def plot_neuron(tp):
    fig2 = py.figure(figsize=(9, 6), dpi=200)
    grid = py.GridSpec(8, 8, hspace=0.1, wspace=0.1)
    ax = fig2.add_subplot(grid[:-1, 1:7], projection='3d')
    ax.set_title('True CSD')
#    ax.scatter(posxyz[0], posxyz[1], posxyz[2], alpha = 0.1, s = 20, 
#               c = vs[:, tp], cmap='PRGn', vmin = -60, vmax = 60, antialiased=True)
#    ax.scatter(posswc[2, ::skip], posswc[3, ::skip], posswc[4, ::skip], alpha = 0.5, s = 8, 
#               c = vs[:1216, tp], cmap='PRGn', vmin = -60, vmax = 60, antialiased=True)
    ax.scatter(posswc[2,::skip], posswc[3,::skip], posswc[4,::skip], alpha = 0.3, s = 12, 
               c = 'grey',  linewidths = 0, antialiased=True)
    ax.scatter(wsp_ele[0], wsp_ele[1], wsp_ele[2], alpha = 1, s = 6, 
               c = pots[tp,:], cmap='PRGn', vmin = -0.2, vmax = 0.2, linewidths = 0.5, 
               edgecolors='black', antialiased=True)
    ax.scatter(posxyz[0], posxyz[1], posxyz[2], alpha = 0.5, s = 15, 
                c = -ks[:, tp]-nas[:,tp]-pas[n,tp], cmap=cm.bwr, vmin = -0.3, vmax = 0.3, antialiased=True)
    ax3 = fig2.add_subplot(grid[-1, 1:-1])

    ax3.plot(time, pots[:,0], label = 'ele bot')
    ax3.plot(time, pots[:,15], label = 'ele mid')
    ax3.plot(time, pots[:,31], label = 'ele top')
    py.legend()
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    
    ax3.set_xlim(0,16)
#    py.ylim(-0.5, 0.5)
    ax3.axvline(time[tp], linestyle = '--')
    ax.view_init(elev = 90, azim=270)
    return fig2

plot_neuron(tp_man)
#%%
'''plots to study currents in different compartments'''
py.figure()
lp = 100
ran=10
py.subplot(311)
py.title('v')
for i in range(ran):
    py.plot(time, vs[i*lp], label = str(i)+':'+ str(np.max(vs[i*lp])))
py.axvline(time[tp_man], linestyle = '--')
py.legend()
py.subplot(312)
py.title('ik')
for i in range(ran):
    py.plot(time, ks[i*lp], label = str(i)+':'+ str(np.max(ks[i*lp])))
py.legend()
py.subplot(313)
py.title('ina')
for i in range(ran):
    py.plot(time, nas[i*lp], label = str(i)+':'+ str(np.max(nas[i*lp])))
py.legend()
#%%