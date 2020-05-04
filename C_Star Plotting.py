# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:58:15 2020

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

def vfront_c_star(vg,vs,fres,fcat,r,c_star):
    r = -1*r*c_star
    
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

def p_c(vg,vs,fres,fcat):
    return (vs*fcat-vg*fres)**2/(vg*vs*(vg+vs))

#%% Download SCC Data

### Specify Directories

local_home = 'C:/Users/Ryan/Documents/Research/Microtubule/' ###Where to send data on local comp
path='/projectnb/qedksk/ryancorn/Microtubule/' ### Location of folder named fig_path on SCC
fig_path = 'Cubic_Growth_12'

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
        if file not in os.listdir(local_home+fig_path+'/'+i):
            print(file)
            sftp.get(file, local_home+fig_path+'/'+i+'/'+file)
        else:
            pass
        
    sftp.chdir(path+fig_path)

### Close SSH 
sftp.close()
shell.close()
ssh.close()

#%% Calculate Front Velocity

def partition(i,n,profile_times,front_pos):
    l = len(front_pos[i])
    fit_time = []
    fit_pos = []
    #fit_time.append(profile_times[i][:l//n])
    #fit_pos.append(front_pos[i][:l//n])
    for j in range(n):
        #print(profile_times[i][(j*l)//n:(j+1)*(l)//n])
        fit_time.append(profile_times[i][(j*l)//n:(j+1)*(l)//n])
        fit_pos.append(front_pos[i][(j*l)//n:(j+1)*(l)//n])
        
    return [fit_time, fit_pos]

local_home = 'C:/Users/Ryan/Documents/Research/Microtubule/Polymer_Growth_7/'
os.chdir(local_home)

convergence = True
n = 50  ##Set number of partitions

c_star = []
r = []
r_conv = []
vfront_third = []
dx = []
dt = []
vfront = []
plot_times = []
cg = []
profile_times = []
break_list = []
error = []
front = []

dir_list = [f for f in os.listdir() if 'C_Star' in f]

for direc in dir_list:
    
    #print(direc)
    os.chdir(local_home+direc)
    
    dt_temp = []
    dx_temp = []
    r_temp = []
    front_pos = []
    front_times = []
    profile_times_temp = []
    r_temp_2 = []
    fit_list = []
    cg_temp = []
    break_list_temp = []
    error_temp = []
    
    counter = 0
    
    input_files = [f for f in os.listdir() if "Microtubule" in f]
    input_files.sort(key = lambda x: int(x.rsplit('_',2)[1]))
    
    for file in input_files:
        
        #print(file)
        
        infile = open(os.getcwd()+"/"+file, 'rb')
        params = pickle.load(infile)
        infile.close()
        
        front_pos.append(params['front_pos'])
        front_times.append(params['front_times'])
        profile_times_temp.append(params['profile_times'])
        dx_temp.append(params['dx'])
        dt_temp.append(params['dt'])
        r_temp.append(params['r'])
        cg_temp.append(params['c_g_array'])
        break_list_temp.append(params['break_param'])
        #error_temp.append(params['error'])
        time = len(np.arange(0,params['T'],params['dt']))
        print('C* = '+str(params['c_star']))
        #print('r = '+str(params['r']))
        #print('Percent Finished = '+str(len(params['front_pos'])/time))
        
        if counter == 0:
            c_star.append(params['c_star'])
            
        counter += 1
    vfront_temp = []
    vfront_second_temp = []
    vfront_third_temp = []
    plot_times_temp = []
    for i in range(len(r_temp)):
        temp_list = []
        temp_times = []
        fit_list.append(partition(i,n,front_times,front_pos))
        for j in range(n):
            #print(fit_list[i][0][j])
            fit = np.polyfit(fit_list[i][0][j], np.multiply(fit_list[i][1][j],dx_temp[i]),1)
            mid_pt = len(fit_list[i][0][j])
            temp_times.append(fit_list[i][0][j][mid_pt//2])
            temp_list.append(fit[0])
            
        vfront_temp.append(temp_list)
        plot_times_temp.append(temp_times)

        if np.abs(temp_list[-1]-temp_list[-2])/temp_list[-2] < 0.01 and convergence == True:
            vfront_third_temp.append(temp_list[-1])
            r_temp_2.append(r_temp[i])
        
        elif break_list_temp[i] == True:
            vfront_third_temp.append(0)
            r_temp_2.append(r_temp[i])
                
        elif break_list_temp[i] == False and np.abs(temp_list[-1]-temp_list[-2])/temp_list[-2] > 0.01:
            vfront_third_temp.append(temp_list[-1])
            r_temp_2.append(r_temp[i])
            print('c* = '+str(c_star[-1])+', r = '+str(r_temp[i])+' did not converge')
    
    front.append(front_pos)
    error.append(error_temp)            
    r_conv.append(r_temp_2)
    r.append(r_temp)
    dt.append(dt_temp)
    dx.append(dx_temp)
    vfront_third.append(vfront_third_temp)
    vfront.append(vfront_temp)
    plot_times.append(plot_times_temp)
    profile_times.append(profile_times_temp)
    cg.append(cg_temp)
    break_list.append(break_list_temp)

#%% Plot Vfront in time (semi-log)

for i in range(len(r)):
    fig1, ax1 = plt.subplots(figsize = (25,15))
    for j in range(len(r[i])):
        ax1.plot(plot_times[i][j],np.abs(np.subtract(vfront[i][j],vfront_c_star(params['vg'],params['vs'],params['fres'],params['fcat'],r[i][j],c_star[i]))), marker = 'o', markersize = 12, label = '$r = $'+str(r[i][j]))
    
        plt.xlabel('Time', fontsize = 30)
        plt.ylabel(r'$|V_{sim} - V_{theory}|$', fontsize = 30)
        plt.ylim([10,10e-4])
        #plt.ylim([0,20])
        plt.title(r'$C^* = $'+str(c_star[i]), fontsize = 35)
        
        plt.legend( loc='best',prop={'size': 18})
        #ax1.set_yscale('log')
        ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)

#%% Plot Vfront in time (log-log)

for i in range(len(r)):
    fig1, ax1 = plt.subplots(figsize = (25,15))
    for j in range(len(r[i])):
        ax1.plot(plot_times[i][j],np.abs(np.subtract(vfront[i][j],vfront_c_star(params['vg'],params['vs'],params['fres'],params['fcat'],r[i][j],c_star[i]))), marker = 'o', markersize = 12, label = '$r = $'+str(r[i][j]))
    
        plt.xlabel('Time', fontsize = 30)
        plt.ylabel(r'$|V_{sim} - V_{theory}|$', fontsize = 30)
        plt.ylim([10,10e-4])
        plt.title(r'$C^* = $'+str(c_star[i]), fontsize = 35)
        
        plt.legend( loc='best',prop={'size': 18})
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)

#%%Plot cg profiles
fig1, ax1 = plt.subplots(figsize = (25,15))
c_choose = 0
r_choose = -1
numPoints = 8
color = iter( plt.cm.jet( np.linspace(0,1,numPoints//1) ) )

for i in range(4,len(cg[c_choose][r_choose])-1,3):
    c = next(color)
    ax1.plot(cg[c_choose][r_choose][i],label = r't = '+str(round(profile_times[c_choose][r_choose][i],3)), color = c)

plt.xlabel(r'Position',fontsize = 30)
plt.ylabel(r'$C_{g}(t,x)$',fontsize = 30)
plt.xlim([-5,14000])
plt.ylim([0,0.03])
plt.legend( loc='best',prop={'size': 18})
#%% Plot vfront in real time
dV = []
for i in range(len(r)):
    fig1, ax1 = plt.subplots(figsize = (25,15))
    dV_temp = []
    for j in range(len(r[i])):
        print(j)
        #ax1.plot(plot_times[i][j],np.abs(np.subtract(vfront[i][j],vfront_c_star(params['vg'],params['vs'],params['fres'],params['fcat'],r[i][j],c_star[i]))), marker = 'o', markersize = 12, label = '$r = $'+str(r[i][j]))
        ax1.plot(plot_times[i][j],vfront[i][j], marker = 'o', markersize = 12, label = '$r = $'+str(r[i][j]))
        dV_temp.append(np.abs(np.subtract(vfront[i][j],vfront_c_star(params['vg'],params['vs'],params['fres'],params['fcat'],r[i][j],c_star[i]))))
        plt.xlabel('Time', fontsize = 30)
        plt.ylabel(r'$|V_{sim} - V_{theory}|$', fontsize = 30)
        plt.title(r'$C^* = $'+str(c_star[i]), fontsize = 35)
        
        plt.legend( loc='best',prop={'size': 18})
        ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)
    dV.append(dV_temp)
    
dV = np.asarray(dV)

#%% Determine exponential decay constant

decay = []
tau = []
log_V = np.log(vfront)

for i in range(len(r)):
    decay_temp = []
    for j in range(len(r[i])):
        l = len(np.where(dV[i][j] > 10e-4)[0])
        try:
            fit = np.polyfit(plot_times[i][j][:l],log_V[i][j][:l],1)
            decay_temp.append(np.abs(fit[0]))
        except:
            decay_temp.append(1)
    decay.append(decay_temp)

for i in range(len(r)):
    tau_temp = []
    for j in range(len(r[i])):
        tau_temp.append(decay[i][j]**(-1))
    tau.append(tau_temp)

#%% Plot front pos and error
fig1, ax1 = plt.subplots(figsize = (25,15))

for i in range(len(error[0])):
    plt.plot(error[0][i], label  = r[0][i])

plt.legend( loc='best',prop={'size': 18})

fig1, ax1 = plt.subplots(figsize = (25,15))

for i in range(len(front[0])):
    plt.plot(front[0][i], label  = r[0][i])

plt.legend( loc='best',prop={'size': 18})
#%% Plot Front Velocity vs Growth Rate
    
growthrange = np.arange(0.1,4.5,0.005)   
    
fig1, ax1 = plt.subplots(figsize = (25,15))

#print('p_c = '+str(p_c(params['vg'],params['vs'],params['fres'],params['fcat'])))
#p_crit = p_c(params['vg'],params['vs'],params['fres'],params['fcat'])
vexact_plot = []

for i in range(len(growthrange)):
    vexact_plot.append(vfront_exact(params['vg'],params['vs'],params['fres'],params['fcat'],growthrange[i]))

c_star_plot = []

for i in range(len(c_star)):
    c_star_plot_temp = []
    for j in range(len(growthrange)):
        c_star_plot_temp.append(vfront_c_star(params['vg'],params['vs'],params['fres'],params['fcat'],growthrange[j],c_star[i]))
    c_star_plot.append(c_star_plot_temp)    
    
#ax1.plot(growthrange, vexact_plot, label = r'Linearized Solution for $C^* = -1.0$')
numPoints = len(vfront_third)

color = iter( plt.cm.jet( np.linspace(0,1,numPoints//1) ) )

plot_points = np.asarray(c_star_plot[0])
idx = np.where(plot_points > 0)[0][0]

#ax1.plot(-1*c_star[0]*growthrange[idx:], plot_points[plot_points > 0], label = r'Exact Solution, $B=0$', color = c)

#plt.plot(np.multiply(r[0],-1*c_star[0]),vfront_third[0], linestyle = '', marker = 'o', markersize = 12, color = c)


for i in range(len(c_star)):
    c = next(color)
    plt.plot(np.multiply(r[i],1),vfront_third[i], linestyle = '', marker = 'o', markersize = 12, color = c, label = 'Simulations')
    #plt.plot(np.multiply(r[i],-1*c_star[i]),vfront_third[i], linestyle = '', marker = 'o', markersize = 12, label = r'Simulation for $B$ = '+str(round(-1/c_star[i],3)), color = c)

plt.axvline(x = p_c(params['vg'],params['vs'],params['fres'],params['fcat']), label = r'$p_c$', linestyle = '--', color = 'k')
plt.plot(np.arange(0,10,0.2),params['vg']*np.ones_like(np.arange(0,10,0.2)), linestyle = '--', color = 'g', label = r'$v_g$')
plt.xlabel(r'Nucleation Rate : $p$', fontsize = 30)
plt.ylabel('Front Velocity : V', fontsize = 30)
plt.ylim([0,params['vg']+0.3])
plt.xlim([0,3])

plt.legend( loc='best',prop={'size': 18})
ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)

# =============================================================================
# fig1, ax1 = plt.subplots(figsize = (25,15))
# ax1.plot(r[0], vfront_third[0], linestyle = '', marker = 'o', markersize = 12, color = 'g', label = 'Simulations')
# plt.axvline(x = p_c(params['vg'],params['vs'],params['fres'],params['fcat']), label = r'$p_c$')
# plt.xlim([0,r[0][-1]+0.2])
# plt.ylim([0,2])
# plt.xlabel('Growth Rate : p', fontsize = 30)
# plt.ylabel('Front Velocity : V', fontsize = 30)
# plt.legend( loc='best',prop={'size': 18})
# =============================================================================

ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)
#%% Plot Vexact - Vsim

V_diff = []

for i in range(len(r)):
    V_diff_temp = []
    for j in range(len(r_conv[i])):
        V_diff_temp.append(np.abs(vfront_c_star(params['vg'],params['vs'],params['fres'],params['fcat'],r[i][j],c_star[i]) - vfront_third[i][j]))#/vfront_c_star(params['vg'],params['vs'],params['fcat'],params['fres'],r[i][j],c_star[i]))
    V_diff.append(V_diff_temp)
    
fig1, ax1 = plt.subplots(figsize = (25,15))

for i in range(len(V_diff)):
    ax1.plot(np.multiply(r[i],-1*c_star[i]),V_diff[i], linestyle = '', marker = 'o', markersize = 12, label = r'$C^* = $'+str(c_star[i]))

plt.xlabel(r'$r$', fontsize = 30)
plt.ylabel(r'$\frac{|V_{linearized} - V_{sim}|}{V_{linearized}}$', fontsize = 30)

plt.legend( loc='best',prop={'size': 18})
ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)

#%% V_theory - V_sim

V_diff = []

for i in range(len(c_star)):
    V_diff_temp = []
    for j in range(len(r[i])):
        if vfront_third[i][j] > 0.1:
            V_diff_temp.append(np.abs(vfront_c_star(params['vg'],params['vs'],params['fres'],params['fcat'],r[i][j],c_star[i]) - vfront_third[i][j])/vfront_c_star(params['vg'],params['vs'],params['fcat'],params['fres'],r[i][j],c_star[i]))
    #print(V_diff_temp)
    V_diff.append(np.average(V_diff_temp))
    
fig1, ax1 = plt.subplots(figsize = (25,15))
ax1.plot(np.divide(-1,c_star),V_diff, linestyle = '', marker = 'o', markersize = 12)

plt.xlabel(r'Cooperativity, $B$', fontsize = 30)
plt.ylabel(r'$\frac{|V(B = 0) - V(B)|}{V(B)}$', fontsize = 30)

#plt.legend( loc='best',prop={'size': 18})
ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)

#%% V_linear  - V_pushed for Polymer case

##we now use B = 1 as linearized prediction

V_diff = []

for i in range(0,len(c_star)):
    V_diff_temp = []
    for j in range(len(r[i])):
        V_diff_temp.append(np.abs(vfront_third[0][j] - vfront_third[i][j])/vfront_third[i][j])
    V_diff.append(np.average(V_diff_temp))
    
fig1, ax1 = plt.subplots(figsize = (25,15))
ax1.plot(np.divide(-1,c_star[1:]),V_diff[1:], linestyle = '', marker = 'o', markersize = 12)

plt.xlabel(r'Cooperativity, $B$', fontsize = 30)
plt.ylabel(r'$\frac{|V(B = 1) - V(B)|}{V(B)}$', fontsize = 30)
plt.ylim([0,0.5])

#plt.legend( loc='best',prop={'size': 18})
ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)
