#!/share/pkg.7/python3/3.6.5/install/bin/python3.6
#$ -j y
import numpy as np
import pickle
import argparse
import math
import time

# =============================================================================
# ############# the following block defines the varaibles that can be passed
# parser = argparse.ArgumentParser()
# parser.add_argument("-T", type=float, help="simulation time",default=10.)
# parser.add_argument("-x", type=float, help="max spatial length",default=10.)
# parser.add_argument("--v_p", type=float, help="polymerisation rate",default=1.)
# parser.add_argument("--v_m", type=float, help="depolymerisation rate",default=1.)
# parser.add_argument("--f_p", type=float, help="catastrophe rate",default=1.)
# parser.add_argument("--f_m", type=float, help="rescue rate",default=1.)
# parser.add_argument("--r_pp", type=float, help="nucleation rate",default=1.)
# parser.add_argument("--r_pm", type=float, help="nucleation rate",default=1.)
# parser.add_argument("--r_mp", type=float, help="nucleation rate",default=1.)
# parser.add_argument("--r_mm", type=float, help="nucleation rate",default=1.)
# parser.add_argument("--d_m", type=float, help="nucleation rate",default=1.)
# parser.add_argument("--d_p", type=float, help="nucleation rate",default=1.)
# parser.add_argument("-k", type=float, help="nucleation carrying capacity",default=10.)
# parser.add_argument("--file_suffix", help="file_suffix",default='')
# parser.add_argument("--path", help="filepath",default="/usr4/ugrad/ryancorn/keisuke_model_output/")
# args = parser.parse_args()

# ########## the variables are assigned values, if variable is not passed, default value is given
# 
# T=args.T
# x=args.x
# v_p=args.v_p
# v_m=args.v_m 
# f_p=args.f_p
# f_m=args.f_m
# r_pp=args.r_pp
# r_pm=args.r_pm
# r_mp=args.r_mp
# r_mm=args.r_mm
# d_m=args.d_m
# d_p=args.d_p
# k=args.k
# file_suffix=args.file_suffix
# path = args.path
# path=args.path
# =============================================================================
use_SCC = True

if use_SCC == True:
    parser = argparse.ArgumentParser()
    parser.add_argument("IPT", help = "Input dictionary containing params")
    args = parser.parse_args()
    
    ipt = args.IPT
    infile = open(ipt, 'rb')
    params_dict = pickle.load(infile)
    infile.close()
    
if use_SCC == False:
    
    print('Not Using SCC Methods') 
    
    v_p = 2 #velocity of growing microtubule
    v_m = 1 #velocity of shrinking microtubule
    r_pp = 3
    f_p = 30 #switching rate from growing to shrinking : fcat   
    f_m = 5  #switching rate from shrinking to growing : fres
    r_pm = 0 #nucleation rate for shrinker from grower
    r_mp = 0 #nucleation rate for grower from shrinker
    r_mm = 0 #nucleation rate for shrinker from shrinker
    d_p = 0 #rate growers die
    d_m = 1 #rate shrinkers die
    k = 1 #carrying capacity
    T = 1000 #total amount of time
    x = v_p*T
    num_saves = 50
    figure = 'C_Star_0'
    ## For vg != vs
    ## q*vg = p*vs
    ## vg/vs = p/q
    p = 2
    q = 1
    c_star = -1.0
    
    #print('r = '+str(r[i]))
    dt = 6*0.01/max(f_p,f_m,r_pp)
    dx = 6*v_p*0.01/max(f_p,f_m,r_pp)
        
    params_dict = {'v_p' : v_p,
        'v_m' : v_m,
        'f_p' : f_p,
        'f_m' : f_m,
        'r_pp' : r_pp,
        'r_pm' : r_pm,
        'r_mp' : r_mp,
        'r_mm' : r_mm,
        'd_p' : d_p,
        'd_m' : d_m,
        'k' : k,
        'T' : T,
        'x' : x,
        'figure' : figure,
        'num_saves' : num_saves,
        'p' : p,
        'q' : q,
        'dt' : dt,
        'dx' : dx,
        'c_star' : c_star,
        'idx' : 0}

## For vg != vs
## q*vg = p*vs
## vg/vs = p/q

v_p = params_dict['v_p']
v_m = params_dict['v_m']
f_p = params_dict['f_p']
f_m = params_dict['f_m']
r_pp = params_dict['r_pp']
r_pm = params_dict['r_pm']
r_mp = params_dict['r_mp']
r_mm = params_dict['r_mm']
d_m = params_dict['d_m']
d_p = params_dict['d_p']
k = params_dict['k']
T = params_dict['T']
x = params_dict['x']
figure = params_dict['figure']
idx = params_dict['idx']
p = params_dict['p']
q = params_dict['q']
dt = params_dict['dt']
dx = params_dict['dx']
c_star = params_dict['c_star']
num_saves = params_dict['num_saves']

def partition(n,profile_times,front_pos):
    l = len(front_pos)
    fit_time = []
    fit_pos = []
    #fit_time.append(profile_times[i][:l//n])
    #fit_pos.append(front_pos[i][:l//n])
    for j in range(n):
        #print(profile_times[i][(j*l)//n:(j+1)*(l)//n])
        fit_time.append(profile_times[(j*l)//n:(j+1)*(l)//n])
        fit_pos.append(front_pos[(j*l)//n:(j+1)*(l)//n])
        
    return [fit_time, fit_pos]

def error_calc(v1,v2):
    if v1 > 0:
        return np.abs(v1 - v2/v1)
    else:
        return 0

def simulate(v_p,v_m,f_p,f_m,r_pp,r_pm,r_mp,r_mm,d_p,d_m,k,figure,idx,dt,dx,p,q,T,x,c_star,num_saves):
       
    x_array = np.arange(0,x,dx)
    front_times = []
    profile_times = []
    front_pos = []
    cp_profile = []
    cm_profile = []
    cp_array = np.zeros(len(x_array))
    cm_array = np.zeros(len(x_array))
    cpTemp = np.zeros_like(cp_array)
    cmTemp = np.zeros_like(cp_array)
    
    front_counter = 0
    profile_counter = 0
    p_counter = 0
    q_counter = 0
    fit_counter = 499
    break_param = False
    
    cp_array[:15] = 1

    start_time = time.time()
    
    for t in np.arange(0,T,dt): 
        
        front_counter = front_counter + 1
        profile_counter = profile_counter + 1
        p_counter = p_counter + 1
        q_counter = q_counter + 1
        
        if t > 0.4*T:
            fit_counter+= 1
        
        if q_counter == q:
            
            cpTemp = np.roll(cp_array, 1)
            
            ### Fixed Boundary Conditions
            cpTemp[0] = cp_array[0] 
            cpTemp[-1] = cp_array[-1]
            
            q_counter = 0
            
        else:
            
            cpTemp = cp_array
        
        if p_counter == p:
            
            cmTemp = np.roll(cm_array, -1)
            
            ### Fixed Boundary Conditions
            
            cmTemp[0] = cm_array[0] 
            cmTemp[-1] = cm_array[-1]

            p_counter = 0
        
        else:
            
            cmTemp = cm_array
            
        
        cpNew = cpTemp - f_p*dt*cpTemp + f_m*dt*cmTemp - d_p*dt*cpTemp
        cpNew = cpNew + r_pp*dt*cpTemp*(1-((cpTemp + cmTemp)/k))*((cpTemp + cmTemp)-c_star) + r_pm*dt*cmTemp*(1-((cpTemp+cmTemp)/k))
        
        cmNew = cmTemp + f_p*dt*cpTemp - f_m*dt*cmTemp - d_m*dt*cmTemp
        cmNew = cmNew + r_mp*dt*cpTemp*(1-((cpTemp + cmTemp)/k)) + r_mm*dt*cmTemp*(1-((cmTemp+cpTemp)/k))
        
        cp_array = cpNew
        cm_array = cmNew
         
        c = cp_array + cm_array
        front_times.append(t)
        front_pos_idx = np.argmin(np.abs(c-10e-4))
        front_pos.append(front_pos_idx)
        front_pos_repeated = []
        front_pos_repeated = [i for i in front_pos if i == front_pos[-1] and i > 0]
         
        if cp_array[-10] > 10e-10:
            cp_array = np.concatenate([cp_array, np.zeros((100),dtype=cp_array.dtype)])
            cm_array = np.concatenate([cm_array, np.zeros((100),dtype=cm_array.dtype)])
        
        if len(front_pos_repeated) > 1000:
            print('Front Velocity Zero, Exiting')
            break_param = True
            break
        
        if fit_counter == 500 and t > 0.4*T:
            n = 50
            vfront = []
            fit_list = partition(n,front_times,front_pos)
        
            for j in range(n):
                fit = np.polyfit(fit_list[0][j], np.multiply(fit_list[1][j],dx),1)
                vfront.append(fit[0])
            
            v2 = np.average(vfront[-5:])
            v1 = np.average(vfront[-10:-5])
            
            if error_calc(v1,v2) < 0.01:
                print('Simulation Converged at T = '+str(t))
                print('Relative Error = '+str(error_calc(v1,v2)))
                break
        
            fit_counter = 0          
           
        if profile_counter == len(np.arange(0,T,dt))//num_saves:
            cp_profile.append(cp_array)
            cm_profile.append(cm_array)
            profile_times.append(t)
            profile_counter = 0
        
    end_time = time.time()
    np.savetxt("Time_loop_"+str(idx), [(np.abs(end_time - start_time))/3600], fmt = '%f') 
    params_sim = {'break_param' : break_param, 'profile_times' : profile_times , 'front_times' : front_times , 'front_pos' : front_pos , 'cp' : cp_profile , 'cm' : cm_profile}
    params_sim.update(params_dict)
    
    outfile = open("Levy_Model_"+figure+"_"+str(idx),'wb')
    pickle.dump(params_sim, outfile)
    outfile.close()
    
simulate(v_p,v_m,f_p,f_m,r_pp,r_pm,r_mp,r_mm,d_p,d_m,k,figure,idx,dt,dx,p,q,T,x,c_star,num_saves)    