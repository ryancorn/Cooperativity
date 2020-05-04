# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 18:56:53 2020

@author: Ryan
"""


import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import paramiko as pm
from IPython import get_ipython

get_ipython().run_line_magic('matplotlib', 'qt')

def vfront_exact(vg,vs,fres,fcat,r):
    r_crit = fcat - (vg/vs)*fres
    
    if r >= r_crit and r < fcat:
        V = vg*(vg*fres-vs*fcat)**2/(vg*(vg*fres-vs*fcat)*(fres+fcat)+(vg+vs)*(vg*fres+vs*fcat)*r-2*(vg+vs)*np.sqrt(vg*fres*fcat*r*(vg*fres-vs*fcat+vs*r)))
    if r < r_crit:
        V = 0
    if r >= fcat:
        V = vg
        
    return V

def linear(a,b,x):
    y = a*x+b
    return y

#%% Download SCC Data

### Specify Directories

local_home = 'C:/Users/Ryan/Documents/Research/Microtubule/' ###Where to send data on local comp
path='/projectnb/qedksk/ryancorn/Microtubule/' ### Location of folder named fig_path on SCC
fig_path = 'Cubic_Growth_4'

os.chdir(local_home)

if fig_path not in os.listdir():
    os.mkdir(local_home+fig_path)
    os.chdir(local_home+fig_path)
    
else:
    os.chdir(local_home+fig_path)

### Setup SSH

ssh = pm.SSHClient()
ssh.set_missing_host_key_policy(pm.AutoAddPolicy())

try:
    ssh.connect('scc1.bu.edu', username='ryancorn', password='YedasPie44')
except pm.SSHException:
    print("Failed to connect.")
    
t = pm.Transport('scc1.bu.edu')
t.connect(username='ryancorn',password='YedasPie44')
sftp = pm.SFTPClient.from_transport(t)
shell = ssh.invoke_shell()
    
sftp.chdir(path+fig_path) 

dir_list = [f for f in sftp.listdir() if 'C_Star_' in f]

for i in dir_list:
    
    ##Make Local Directory
    if i not in os.listdir():
        print('Making New Directory at '+str(local_home+fig_path+'/'+i))
        os.mkdir(local_home+fig_path+'/'+i)
    
    ##Download files
    sftp.chdir(path+fig_path+'/'+i)
    
    input_files = [f for f in sftp.listdir() if 'Microtubules_C_Star' in f]

    for file in input_files:
        print(file)    
        sftp.get(file, local_home+fig_path+'/'+i+'/'+file)
        
    sftp.chdir(path+fig_path)

### Close SSH 
sftp.close()
shell.close()
ssh.close()

#%% Calculate Spatial Decay Rate

local_home = 'C:/Users/Ryan/Documents/Research/Microtubule/Cubic_Growth_4/'
os.chdir(local_home)

decay = []
r = []
c_star = []
rg_profile_1 = []
rs_profile_1 = []
rg_profile_2 = []
rs_profile_2 = []
rg_profile_3 = []
rs_profile_3 = []

dir_list = [f for f in os.listdir() if 'C_Star' in f]

for direc in dir_list:
    
    r_temp = []
    rg_temp_1 = []
    rs_temp_1 = []
    rg_temp_2 = []
    rs_temp_2 = []    
    rg_temp_3 = []
    rs_temp_3 = []
    decay_temp = []
    c_star_temp = []
    
    print(direc)
    os.chdir(local_home+direc)

    input_files = [f for f in os.listdir() if "Microtubule" in f]
    input_files.sort(key = lambda x: int(x.rsplit('_',2)[1]))
    
    for file in input_files:
        
        print(file)
        
        infile = open(os.getcwd()+"/"+file, 'rb')
        params = pickle.load(infile)
        infile.close()
        
        rg_temp_1.append(params['rg_profile'][0]) 
        rs_temp_1.append(params['rs_profile'][0])
        rg_temp_2.append(params['rg_profile'][1]) 
        rs_temp_2.append(params['rs_profile'][1])        
        rg_temp_3.append(params['rg_profile'][2]) 
        rs_temp_3.append(params['rs_profile'][2])
        rg = params['rg_profile']
        rs = params['rs_profile']
        c_star_temp.append(params['c_star'])
        r_temp.append(params['r'])

        fit_range_rg = [[[] for y in range(len(rg[x]))] for x in range(len(rg))]
        fit_range_rs = [[[] for y in range(len(rs[x]))] for x in range(len(rs))]
    
        for i in range(len(rg)):
            for j in range(len(rg[i])):
                fit_list_rg = [x for x in rg[i][j] if x < 0.2 and x > 10e-5]
                fit_list_rs = [x for x in rs[i][j] if x < 0.2 and x > 10e-5]
                
                if len(fit_list_rg) < 20:
                    print('Potentially underfitting for '+file+', i = '+str(i)+', j = '+str(j))
                    
                fit_range_rg[i][j] = fit_list_rg
                fit_range_rs[i][j] = fit_list_rs
            
        for i in range(len(rg)):
            fits = []
            decay_temp_temp = []
            for j in range(len(rg[i])):
                fit = np.polyfit(np.arange(len(fit_range_rg[i][j])),np.log(fit_range_rg[i][j]),1)                            
                decay_temp_temp.append(fit[0])
        #print(np.average(decay_temp_temp))
        decay_temp.append(np.average(decay_temp_temp))
    
    c_star.append(c_star_temp)
    r.append(r_temp)
    decay.append(decay_temp)
    rg_profile_1.append(rg_temp_1)
    rs_profile_1.append(rs_temp_1)
    rg_profile_2.append(rg_temp_2)
    rs_profile_2.append(rs_temp_2)
    rg_profile_3.append(rg_temp_3)
    rs_profile_3.append(rs_temp_3)
    
#%% Testing Calculating spatial decay rate using only last profile

local_home = 'C:/Users/Ryan/Documents/Research/Microtubule/Cubic_Growth_4/'
os.chdir(local_home)

decay = []
r = []
c_star = []
rg_profile_1 = []
rs_profile_1 = []
rg_profile_2 = []
rs_profile_2 = []
rg_profile_3 = []
rs_profile_3 = []

dir_list = [f for f in os.listdir() if 'C_Star' in f]

for direc in dir_list:
    
    r_temp = []
    rg_temp_1 = []
    rs_temp_1 = []
    rg_temp_2 = []
    rs_temp_2 = []    
    rg_temp_3 = []
    rs_temp_3 = []
    decay_temp = []
    c_star_temp = []
    
    print(direc)
    os.chdir(local_home+direc)

    input_files = [f for f in os.listdir() if "Microtubule" in f]
    input_files.sort(key = lambda x: int(x.rsplit('_',2)[1]))
    
    for file in input_files:
        
        print(file)
        
        infile = open(os.getcwd()+"/"+file, 'rb')
        params = pickle.load(infile)
        infile.close()
        
        rg_temp_1.append(params['rg_profile'][0]) 
        rs_temp_1.append(params['rs_profile'][0])
        rg_temp_2.append(params['rg_profile'][1]) 
        rs_temp_2.append(params['rs_profile'][1])        
        rg_temp_3.append(params['rg_profile'][2]) 
        rs_temp_3.append(params['rs_profile'][2])
        rg = params['rg_profile'][0][-1]
        c_star_temp.append(params['c_star'])
        r_temp.append(params['r'])

        fit_range_rg = [x for x in rg if x < 0.2 and x > 10e-5]
                
        if len(fit_list_rg) < 20:
            print('Potentially underfitting for '+file)
            
        fit = np.polyfit(np.arange(len(fit_range_rg)),np.log(fit_range_rg),1)                            
        decay_temp.append(fit[0])
    
    c_star.append(c_star_temp)
    r.append(r_temp)
    decay.append(decay_temp)
    rg_profile_1.append(rg_temp_1)
    rs_profile_1.append(rs_temp_1)
    rg_profile_2.append(rg_temp_2)
    rs_profile_2.append(rs_temp_2)
    rg_profile_3.append(rg_temp_3)
    rs_profile_3.append(rs_temp_3)
#%% Plot some profiles

plot = 2

if plot == 0:

    for l in range(len(rg_profile_1)):
        fig1, ax1 = plt.subplots(figsize = (25,15))
        for i in range(len(rg_profile_1[l])):
            plt.plot(rg_profile_1[l][i][2],label = 'r = '+str(r[l][i]))
                

        ax1.set_xlim([0,750])
        ax1.set_ylim([10e-10,.2])
        ax1.set_yscale('log')
        ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14) 
        plt.legend( loc='best',prop={'size': 18})

if plot == 1:


    for l in range(len(rg_profile_2)):
        fig2, ax2 = plt.subplots(figsize = (25,15))
        for i in range(len(rg_profile_2[l])):
            plt.plot(rg_profile_2[l][i][0],label = 'r = '+str(r[l][i]))
                

        ax2.set_xlim([0,750])
        ax2.set_ylim([10e-10,.2])
        ax2.set_yscale('log')
        ax2.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14) 
        plt.legend( loc='best',prop={'size': 18})
    
if plot  == 2:
            
# =============================================================================
#     for l in range(len(rg_profile_3)):
#         fig3, ax3 = plt.subplots(figsize = (25,15))
#         for i in range(len(rg_profile_3[l])):
#             plt.plot(rg_profile_3[l][i][0],label = 'r = '+str(r[l][i]))
# =============================================================================
        fig3, ax3 = plt.subplots(figsize = (25,15))
        for i in range(len(rg_profile_3[-1])):        
            plt.plot(rg_profile_3[0][i][4], label = 'r = '+str(r[-1][i]))
    
        fig3.suptitle(r'$\ln(\rho_g)$ for $C^* = $'+str(c_star[0][0]),fontsize = 30)
        ax3.set_xlim([0,750])
        ax3.set_ylim([10e-10,.2])
        ax3.set_yscale('log')
        ax3.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14) 
        plt.legend( loc='best',prop={'size': 18})  
         
        


#%% Plot Spatial Decay Rate

fig1, ax1 = plt.subplots(figsize = (25,15))

for j in range(len(decay)):
    plt.plot(r[j],decay[j], linestyle = '', marker = 'o', markersize = 12, label = r'Spatial Decay Rate for $C^*$ = '+str(c_star[j][0]))

plt.xlabel('Growth Rate : r', fontsize = 30)
plt.ylabel('Spatial Decay Rate', fontsize = 30)

plt.legend( loc='best',prop={'size': 18})
ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)    
