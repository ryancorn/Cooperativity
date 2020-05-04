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

def vfront_exact(v_p,v_m,f_p,f_m,r_pp):
    if r_pp < f_p:
        V = (2*v_m*np.sqrt(f_p*r_pp)-(f_p+r_pp)*v_m+f_m*v_p)/(f_m+f_p+r_pp-2*np.sqrt(f_p*r_pp))
    if r_pp >= f_p:
        V = v_p
    return V

def vfront_c_star(v_p,v_m,f_p,f_m,r_pp,c_star):
    r_pp = -1*r_pp*c_star
    if r_pp < f_p:
        V = (2*v_m*np.sqrt(f_p*r_pp)-(f_p+r_pp)*v_m+f_m*v_p)/(f_m+f_p+r_pp-2*np.sqrt(f_p*r_pp))
    if r_pp >= f_p:
        V = v_p
    return V

def vfront_d_m_c_star(v_p,v_m,f_p,f_m,r_pp,c_star,d_m):
    r_crit = f_p - (f_p*f_m)/(d_m+f_m)
    r_pp = -1*r_pp*c_star
    #r_pp = r_pp/4*(1-2*c_star)
    c11 = r_pp-f_p
    c22 = -d_m-f_m
    c12 = f_m
    c21 = f_p
    v_m = -v_m
    if r_pp <= r_crit:
        return 0
    if r_pp < f_p and r_pp > r_crit:
        return ((c22*v_p - c11*v_m)+(c22*v_p+c11*v_m)*np.sign(v_p-v_m)*np.sqrt((c12*c21)/(c21*c12 - c11*c22)))/((c22 - c11)+(c22+c11)*np.sign(v_p-v_m)*np.sqrt((c12*c21)/(c21*c12 - c11*c22)))
    if r_pp >= f_p:
        return v_p

def linear(a,b,x):
    y = a*x+b
    return y

#%% Download SCC Data

### Specify Directories

local_home = 'C:/Users/Ryan/Documents/Research/Microtubule/' ###Where to send data on local comp
path='/projectnb/qedksk/ryancorn/Levy_Model/' ### Location of folder named fig_path on SCC
fig_path = 'Levy_Cubic_8'

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
    
    input_files = [f for f in sftp.listdir() if 'Levy_Model_C_Star' in f]
    in_dir = os.listdir(local_home+fig_path+'/'+i)
    
    for file in input_files:
        if file not in in_dir:
            print(file)
            sftp.get(file, local_home+fig_path+'/'+i+'/'+file)
        else: 
            print('File '+str(file)+' already downloaded.')
            pass
        
    sftp.chdir(path+fig_path)

### Close SSH 
sftp.close()
shell.close()
ssh.close()

#%% Calculate Front Velocity

local_home = 'C:/Users/Ryan/Documents/Research/Microtubule/Levy_Cubic_8/'
os.chdir(local_home)

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

c_star = []
r_pp = []
front_loc = []
front_loc_time = []
profile_times = []
cp = []
cm = []
vfront_first = []
vfront_second = []
vfront_third = []
dx = []
dt = []
break_list = []
vfront = []
front = []

dir_list = [f for f in os.listdir() if 'C_Star' in f]
dir_list.sort(key = lambda x: int(x.rsplit('_',1)[1]))

for direc in dir_list:
    
    print(direc)
    os.chdir(local_home+direc)
    
    dt_temp = []
    dx_temp = []
    r_pp_temp = []
    front_pos = []
    front_times = []
    cp_temp = []
    cm_temp = []
    profile_times_temp = []
    r_pp_temp_2 = []
    break_list_temp = []
    
    counter = 0

    input_files = [f for f in os.listdir() if "Levy_Model" in f]
    input_files.sort(key = lambda x: int(x.rsplit('_',2)[1]))
    
    for file in input_files:
        
        print(file)
        
        infile = open(os.getcwd()+"/"+file, 'rb')
        params = pickle.load(infile)
        infile.close()
        front_pos.append(params['front_pos'])
        profile_times_temp.append(params['profile_times'])
        front_times.append(params['front_times'])
        dx_temp.append(params['dx'])
        dt_temp.append(params['dt'])
        r_pp_temp.append(params['r_pp'])
        cp_temp.append(params['cp'])
        cm_temp.append(params['cm'])
        break_list_temp.append(params['break_param'])
        if counter == 0:
            c_star.append(params['c_star'])
        
        counter += 1
        
    vfront_first_temp = []
    vfront_second_temp = []
    vfront_third_temp = []
    vfront_temp = []
    
    print(params['c_star'])
    
    ##Fit first third of points
            
    for i in range(0,len(r_pp_temp)):
        fit = np.polyfit(front_times[i][0:len(front_pos[i])//3], np.multiply(front_pos[i][0:len(front_pos[i])//3],dx_temp[i]), 1)
        vfront_first_temp.append(fit[0])
                    
    ##Fit middle third of points
                
    for i in range(0,len(r_pp_temp)):
        fit = np.polyfit(front_times[i][(len(front_pos[i]))//3:(2*len(front_pos[i]))//3], np.multiply(front_pos[i][(len(front_pos[i]))//3:(2*len(front_pos[i]))//3],dx_temp[i]), 1)
        vfront_second_temp.append(fit[0])
                    
    ##Fit latter third of points
        
    for i in range(0,len(r_pp_temp)):
        
        if break_list_temp[i] == True:
            vfront_third_temp.append(0)
            r_pp_temp_2.append(r_pp_temp[i])
        
        fit_list = []
        temp_list = []
        temp_times = []
        n = 50
        fit_list.append(partition(i,n,front_times,front_pos))
        for j in range(n):
            #print(fit_list[i][0][j])
            fit = np.polyfit(fit_list[0][0][j], np.multiply(fit_list[0][1][j],dx_temp[i]),1)
            mid_pt = len(fit_list[0][0][j])
            temp_times.append(fit_list[0][0][j][mid_pt//2])
            temp_list.append(fit[0])
            
        v1 = np.average(temp_list[-5:])
        v2 = np.average(temp_list[-10:-5])
        vfront_temp.append(temp_list)
        if np.abs(v1-v2)/v2 < 0.01 and break_list_temp[i] == False:
            r_pp_temp_2.append(r_pp_temp[i])
            vfront_third_temp.append(temp_list[-1])
        elif break_list_temp[i] == True:
            r_pp_temp_2.append(r_pp_temp[i])
            vfront_third_temp.append(0)
        elif v1 or v2 == 0:
            r_pp_temp_2.append(r_pp_temp[i])
            vfront_third_temp.append(0)
        elif np.abs(v1-v2)/v2 > 0.01 and break_list_temp[i] == False:
            r_pp_temp_2.append(r_pp_temp[i])
            vfront_third_temp.append(temp_list[-1])
            print('c* = '+str(c_star[-1])+', r = '+str(r_temp[i])+' did not converge')
            

    front.append(front_pos)     
    vfront.append(vfront_temp)        
    r_pp.append(r_pp_temp_2)
    profile_times.append(profile_times_temp)
    front_loc.append(front_pos)
    front_loc_time.append(front_times)
    dt.append(dt_temp)
    dx.append(dx_temp)
    cp.append(cp_temp)
    cm.append(cm_temp)
    vfront_first.append(vfront_first_temp)
    vfront_second.append(vfront_second_temp)
    vfront_third.append(vfront_third_temp)
    break_list.append(break_list_temp)
    
#%% Plot 1st/2nd/3rd Velocity Fits

for i in range(len(vfront_third)):
    fig1, ax1 = plt.subplots(figsize = (25,15))
    
    plt.plot(r[i],vfront_first[i], linestyle = '', marker = 'o', markersize = 12, label = r'First Third for $C^*$ = '+str(c_star[i][0]))
    plt.plot(r[i],vfront_second[i], linestyle = '', marker = 'o', markersize = 12, label = r'Middle Third for $C^*$ = '+str(c_star[i][0]))
    plt.plot(r[i],vfront_third[i], linestyle = '', marker = 'o', markersize = 12, label = r'Last Third for $C^*$ = '+str(c_star[i][0]))

    plt.xlabel('Growth Rate : r', fontsize = 30)
    plt.ylabel('Front Velocity : V', fontsize = 30)

    plt.legend( loc='best',prop={'size': 18})
    ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)

#%% Plot Profiles
fig1, ax1 = plt.subplots(figsize = (25,15))

r_choose = 1

for j in range(len(c_star)):
    c = np.add(cp[j][r_choose][-1],cm[j][r_choose][-1])
    front_pos_idx = np.argmin(np.abs(c-0.1))
    print(c[front_pos_idx])
    cp_plot = cp[j][r_choose][-1][front_pos_idx-500:front_pos_idx+500]
    cm_plot = cm[j][r_choose][-1][front_pos_idx-500:front_pos_idx+500]
    c_plot = np.add(cp_plot, cm_plot)
    #plt.plot(np.divide(cm_plot,c_plot), label = 'C* = '+str(c_star[j]))
    plt.plot(cm_plot, label = 'C* = '+str(c_star[j]))

plt.legend( loc='best',prop={'size': 18})
plt.title(r'$c_-$',fontsize = 30)

fig1, ax1 = plt.subplots(figsize = (25,15))


for j in range(len(c_star)):
    c = np.add(cp[j][r_choose][-1],cm[j][r_choose][-1])
    front_pos_idx = np.argmin(np.abs(c-0.1))
    print(c[front_pos_idx])
    cp_plot = cp[j][r_choose][-1][front_pos_idx-500:front_pos_idx+500]
    cm_plot = cm[j][r_choose][-1][front_pos_idx-500:front_pos_idx+500]
    c_plot = np.add(cp_plot, cm_plot)
    #plt.plot(np.divide(cp_plot,c_plot), label = 'C* = '+str(c_star[j]))
    plt.plot(cp_plot, label = 'C* = '+str(c_star[j]))

plt.legend( loc='best',prop={'size': 18})
plt.title(r'$c_+$',fontsize = 30)

#fig1, ax1 = plt.subplots(figsize = (25,15))

# =============================================================================
# for i in np.arange(len(cp[c_choose][r_choose])):
#     plt.plot(cp[c_choose][r_choose][i])
# 
# plt.title('cp',fontsize = 30)
#     
# fig1, ax1 = plt.subplots(figsize = (25,15))
# 
# for i in np.arange(len(cp[c_choose][r_choose])):
#     plt.plot(cm[c_choose][r_choose][i])
# 
# plt.title('cm', fontsize = 30)
# =============================================================================
#%%Fit front position
c_choose = -1
r_choose = -1

c = []
front = []

for i in range(len(cp[c_choose][r_choose])):
    #c.append(cp[c_choose][r_choose][i]+cm[c_choose][r_choose][i])
    c = np.add(cp[c_choose][r_choose][i],cm[c_choose][r_choose][i])
    c_max = np.argmax(c)
    #print(len(c_fitting)+len(c[:c_max]))
    #print(c_max)
    c_fitting = c[:]
    front_pos_idx_new = np.argmin(np.abs(c_fitting-0.5))
    front.append(front_pos_idx_new)


#%% Plot front position
fig1, ax1 = plt.subplots(figsize = (25,15))
c_choose = -1
r_choose = -1

plt.plot(front_loc_time[c_choose][r_choose],front_loc[c_choose][r_choose],marker = 'o', markersize = 12, linestyle = '')

#%% Plot Front Velocity vs Growth Rate
    
fig1, ax1 = plt.subplots(figsize = (25,15))
growthrange = np.arange(0,15,0.005)   
# =============================================================================
# vexact_plot = []
# 
# for i in range(len(growthrange)):
#     vexact_plot.append(vfront_exact(params['v_p'],params['v_m'],params['f_p'],params['f_m'],growthrange[i]))
# 
# ax1.plot(growthrange, vexact_plot, label = r'Exact Solution $C^* = -1.0$')
# 
# =============================================================================
numPoints = 6
color = iter( plt.cm.jet( np.linspace(0,1,numPoints//1) ) )
v_linear_plot = []
shifts = []

for i in range(0,len(c_star),2):
    c = next( color )
    growthrange = np.arange(0,15,0.005)   
    shift_list = []
    v_linear_plot = []
    for j in range(len(growthrange)):
        v_linear_plot.append(vfront_d_m_c_star(params['v_p'],params['v_m'],params['f_p'],params['f_m'],growthrange[j],c_star[i],params['d_m']))
    
    for l in range(len(r_pp[i])):
        arg = np.argmin(np.abs(np.subtract(v_linear_plot,vfront_third[i][l])))
        shift_list.append(-1*c_star[i]*np.abs(growthrange[arg] - r_pp[i][l]))
        
    shift = np.average(shift_list) #shift is average deviation in r for a given C*
    shifts.append(shift)
    #shift = np.multiply(-1*c_star[i], shift)
    #print(shift)
    #growthrange = growthrange - shift
    #print(shift)
    ax1.plot(np.multiply(r_pp[i],-1*c_star[i]),vfront_third[i], linestyle = '', marker = 'o', markersize = 12, label = r'Simulation for $B$ = '+str(round(-1/c_star[i],3)), color = c)
    #ax1.plot(np.multiply(r_pp[i],-1*c_star[i]),vfront_third[i], linestyle = '', marker = 'o', markersize = 12, label = r'Simulation for $B$ = '+str(round(-1/c_star[i],3)), color = c)

c = next(color)
idx = min(np.where(np.asarray(v_linear_plot) > 0)[0])
growthrange = growthrange[idx:]
v_linear_plot = v_linear_plot[idx:]
ax1.plot(np.multiply(-1*c_star[i],growthrange), v_linear_plot, label = r'Exact Solution, $B=0$', color = c)

plt.plot(np.arange(0.1,3,0.005),params['v_p']*np.ones_like(np.arange(0.1,3,0.005)), linestyle = '--', color = 'k', label = r'$v_{+}$')
plt.xlabel(r'Nucleation Rate : $\tilde{r}_{++}$', fontsize = 30)
plt.ylabel('Front Velocity : V', fontsize = 30)
plt.xlim([0,3])
plt.ylim([0,params['v_p']+0.3])
#plt.legend( loc='best',prop={'size': 18},bbox_to_anchor=(0.82, 0.85))
plt.legend( loc='best',prop={'size': 18})
ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)

#%% Plot shift as function of C*

fig1, ax1 = plt.subplots(figsize = (25,15))
growthrange = np.arange(0.1,17,0.005)

# =============================================================================
# shift_list = []
# line_plot = []
# 
# for i in range(len(r_pp)):
#     v_linear_plot = []
#     shift_list_temp = []
#     for j in range(len(growthrange)):
#         v_linear_plot.append(vfront_d_m_c_star(params['v_p'],params['v_m'],params['f_p'],params['f_m'],growthrange[j],c_star[i],params['d_m']))
#     for j in range(len(r_pp[i])):
#         idx = min(np.where(np.asarray(v_linear_plot) > 0)[0])
#         growthrange = growthrange[idx:]
#         v_linear_plot = v_linear_plot[idx:]
#         arg = np.argmin(np.abs(np.subtract(v_linear_plot,vfront_third[i][j])))
#         #print(growthrange[arg])
#         #print(r_pp[i][j])
#         shift = np.abs(growthrange[arg] - r_pp[i][j])
#         shift_list_temp.append(shift)
#     shift_list.append(np.average(shift_list_temp))
# =============================================================================

# =============================================================================
# fits = np.polyfit(np.log(np.divide(-1,c_star)), np.log(shift_list), 1)
# plot_range = np.arange(np.log(-1/c_star[0]),np.log(-1/c_star[-1]+40),0.01)
# 
# ax1.plot(np.log(np.divide(-1,c_star)),np.log(shift_list), linestyle = '', marker = 'o', markersize = 12, label = r'Simulations', color = 'b')
# ax1.plot(np.log(plot_range), fits[0]*np.log(plot_range) + fits[1])
# =============================================================================

ax1.plot(np.divide(-1,c_star), shifts , color = 'b', markersize = 12, marker = 'o', linestyle = '')
plt.xlabel(r'Cooperativity, $B$', fontsize = 30)
plt.ylabel(r'Average Shift Size, $\langle\Delta \tilde{r}\rangle$', fontsize = 30)
#plt.xlim([-1,0])
#plt.ylim([0,6])
#plt.legend( loc='best',prop={'size': 18})
ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)
#ax1.set_yscale('log')
#ax1.set_xscale('log')

#%% Plot Vexact - Vsim

V_diff = []

for i in range(len(r_pp)):
    V_diff_temp = []
    for j in range(len(r_pp[i])):
        V_diff_temp.append(np.abs(vfront_c_star(params['v_p'],params['v_m'],params['f_p'],params['f_m'],r_pp[i][j],c_star[i]) - vfront_third[i][j])/vfront_c_star(params['v_p'],params['v_m'],params['f_p'],params['f_m'],r_pp[i][j],c_star[i]))
    V_diff.append(V_diff_temp)
    
fig1, ax1 = plt.subplots(figsize = (25,15))

for i in range(len(V_diff)):
    plt.plot(r_pp[i],V_diff[i], linestyle = '', marker = 'o', markersize = 12, label = r'$C^* = $'+str(c_star[i]))

plt.xlabel(r'$r_{++}$', fontsize = 30)
plt.ylabel(r'$\frac{|V_{linearized} - V_{sim}|}{V_{linearized}}$', fontsize = 30)

plt.legend( loc='best',prop={'size': 18})
ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)

#%% Plot Vexact - Vsim for different B
V_diff = []

f_p = params['f_p']
f_m = params['f_m']
v_p = params['v_p']
v_m = params['v_m']
d_m = params['d_m']

def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

##using linearized prediction
for i in range(len(r_pp)):
    V_diff_temp = np.zeros_like(r_pp[i])
    for j in range(len(r_pp[i])):
        V_diff_temp[j] = (np.abs(vfront_d_m_c_star(params['v_p'],params['v_m'],params['f_p'],params['f_m'],r_pp[i][j],c_star[i],params['d_m']) - vfront_third[i][j])/vfront_third[i][j])#np.average(vfront_third[i]))
    #V_fixed = reject_outliers(V_diff_temp)
    V_diff.append(np.mean(V_diff_temp))
    #V_diff.append(V_diff_temp[5])
        
fig1, ax1 = plt.subplots(figsize = (25,15))
plt.plot(np.divide(-1,c_star[:-2]),V_diff[:-2], linestyle = '', marker = 'o', markersize = 12)

plt.xlabel(r'Cooperativity, $B$', fontsize = 30)
plt.ylabel(r'$\frac{|V(B = 0) - V(B)|}{V(B)}$', fontsize = 30)

#plt.legend( loc='best',prop={'size': 18})
ax1.tick_params(axis='both', which='both', labelsize=20, width = 3,length  = 14)