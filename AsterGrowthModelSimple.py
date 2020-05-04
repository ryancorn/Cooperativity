#!/share/pkg.7/python3/3.6.5/install/bin/python3.6
#$ -j y
import numpy as np
#import cy_test
import pickle
import argparse
import time
import sys

useSCC = True 

if useSCC == True:
    
    print('Using SCC Methods')
    parser = argparse.ArgumentParser()
    parser.add_argument("IPT", help = "Input dictionary containing params")
    args = parser.parse_args()

    ipt = args.IPT
    infile = open(ipt, 'rb')
    params_dict = pickle.load(infile)
    infile.close()

    T = params_dict['T']
    l = params_dict['l']
    x = params_dict['x']
    vg = params_dict['vg']
    vs = params_dict['vs']
    fcat = params_dict['fcat']
    fres = params_dict['fres']
    r = params_dict['r']
    k = params_dict['k']
    idx = params_dict['idx']
    num_saves = params_dict['num_saves']
    figure = params_dict['figure']
    p = params_dict['p']
    q = params_dict['q']
    dt = params_dict['dt']
    dx = params_dict['dx']
    dl = params_dict['dl']
    c_star = params_dict['c_star']
    nucleation_type = params_dict['nucleation_type']

if useSCC == False: 
    
    print('Not Using SCC Methods')
    os.chdir('C:\\Users\\Ryan\\Documents\\Research\\Microtubule')
    vg = 2 #velocity of growing microtubule ASSUMED TO BE TYPE INT
    vs = 1 #velocity of shrinking microtubule ASSUMED TO BE TYPE INT
    fcat = 3 #switching rate from growing to shrinking : fcat   
    fres = 1  #switching rate from shrinking to growing : fres
    r = 1.8  #nucleation rate for grower from grower
    k = 1 #carrying capacity
    T = 750 #total amount of time
    x = 25
    l = x
    c_star = -0.6
    nucleation_type = 1
    
    print((vs*fcat-vg*fres)**2/(vg*vs*(vg+vs)))
    
    if nucleation_type==1:
        dt = 0.1/max(fcat,fres,r,-r/c_star)
        #dt = 0.1/max(fcat,fres,r)
        dl=vg*dt
        dx=dl
    elif nucleation_type==2 or nucleation_type==3:
        #dt = 0.05/max(fcat,fres,r,-r/c_star)   
        dt = 0.1/max(fcat,fres)
        dl=vg*dt
        if r*dt*dl*k/4 >0.1:
            print("rdt dl was too big")
            print (r*dt*dl*k/4)
            rescale_Factor=0.1/(r*dt*dl*k/4)
            dt=dt*np.sqrt(rescale_Factor)
            dl=vg*dt
            print(r*dt*dl*k/4)
            #exit
        dx=dl
        
    num_saves = len(np.arange(0,T,dt))//10
    figure = 'test'
    ## For vg != vs
    ## q*vg = p*vs
    ## vg/vs = p/q
    p = 2
    q = 1
    idx = 0
    
    input_params = {'vg' : vg,
                    'vs' : vs,
                    'fcat' : fcat,
                    'fres' : fres,
                    'r' : r,
                    'k' : k,
                    'T' : T,
                    'x' : x,
                    'l' : l,
                    'num_saves' : num_saves,
                    'figure' : figure,
                    'p' : p,
                    'q' : q,
                    'idx' : idx}

##define function to calculate growing plus end density
    
CFL = np.round(vg*dt/dl, 5)
r = -1*c_star*r
B = -1/c_star

def calculate_c_g(c_g,rg_array):
    #have to flip the array before doing the sums of the diagonals because trace only does main diagonals, not anti-diagonals             
    rg_flip = np.fliplr(rg_array)
    X = rg_array.shape[0]
    for i in range(1, X+1): #len(x_array) is included therefore the +1 is needed
        c_g[i-1] = np.trace(rg_flip, offset = X-i, axis1=0, axis2=1)
    return c_g

def calculate_c_s(c_s,rs_array):
    #have to flip the array before doing the sums of the diagonals because trace only does main diagonals, not anti-diagonals             
    rs_flip = np.fliplr(rs_array)
    X = len(rs_array)
    for i in range (1, X+1): #len(x_array) is included therefore the +1 is needed
        c_s[i-1] = np.trace(rs_flip, offset = X-i, axis1=0, axis2=1)
    return c_s

#1 function calculates polymer density    
def calculate_D(rg_array):
    X = rg_array.shape[0]
    c=calculate_c_g(np.zeros(X),rg_array)
    rg_x_sum=np.sum(rg_array,axis=1) ## The sum over l values of rg_array
    D = np.zeros(X)
    D[0]=rg_x_sum[0]
    for i in range (1, X):
        D[i]=D[i-1]-c[i-1]+rg_x_sum[i]
    return D
#2

#1 this calculates a "Cg" that is weighted by the length of the microtubules: Initial Conditions will need microtubules of nonzero length for this give you a growing aster!!   
def calculate_W(c,rg_array,l_array):
    l_weighted_array=np.multiply(rg_array,l_array)
    l_weighted_c=calculate_c_g(c,l_weighted_array) 
    return l_weighted_c  

def roll_left(A, Yshape):
    return np.concatenate((A[:,1:], A[:,0].reshape(Yshape)), axis=1)

def roll_right(A, Yshape):
    return np.concatenate((A[:,-1].reshape(Yshape),A[:,:-1]), axis=1)

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

if CFL != 1 and vg!= 0:
    print('CFL != 1, exiting!')
    print('CFL = '+str(CFL))
    sys.exit()

if x < 10:
    x = 0.1*T
    l = x  

if k != 1:
    print('Carrying Capcity != 1, exiting!')
    sys.exit()

print('dt =',dt)
print('dl =',dl)
print('dx =',dx)

#x_array contains the position of negative microtubule lengths from [0,x) with spacing of dx
x_array = np.arange(0, x, dx)

#l_array contains the length from [0,l) with spacing of dl
l_array = np.arange(0, l, dl)

#initialize simulation timer
timer = np.arange(0,T,dt)

X = len(x_array)
L = len(l_array)

leftBoundary = 1
rightBoundary = 0

#create array to hold profiles and times
profile_times = []
front_pos = []
front_times = []
rg_profile_1 = []
rg_profile_2 = []
rg_profile_3 = []
rs_profile_1 = []
rs_profile_2 = []
rs_profile_3 = []

#create array for densities
c_g_array = []
#c_s_array = []

#rg_array contains the values for rho(t,x,l) for growing microtubules 
##rs and rg arrays have form [lattice site, length (multiple of dl)]
rg_array = np.zeros((X, L))
#rs_array contains the values for rho(t,x,l) for shrinking microtubules
rs_array = np.zeros((X, L))

shapey, shapex = rg_array.shape
Yshape = (shapex, 1)

rgNew = np.zeros_like(rg_array)
rsNew = np.zeros_like(rs_array)
rgTemp = np.zeros_like(rg_array)
rsTemp  = np.zeros_like(rs_array)
c_g = np.zeros(X)
D = np.zeros(X) ### polymer density
W = np.zeros(X) ### A new "C_g" where nucleation is proportional to the length of the microtubules.
#c_s = np.zeros(X)
rg_shifted = np.zeros_like(rg_array)
rs_shifted = np.zeros_like(rg_array)
rg_flip = np.zeros_like(rg_array)

#initial condition
if nucleation_type == 3:
    rg_array[:5,0] = 0.0001
    rg_array[0,:5] = 0.0001

elif nucleation_type != 3:
    rg_array[:2,0] = 0.0001


c_g_array.append(calculate_c_g(c_g,rg_array))
    
#c_s_array[0] = calculate_c_s(c_s,rs_array,X,dx)
#c_array[0] = cy_test.function_2(rg_array)
#rg_profile[0] = rg_array
#rs_profile[0] = rs_array
nPoints = 250

profile_counter = 0
q_counter = 0
p_counter = 0
save_counter = 1
rg_profile_counter = 0
extend_count = 0
fit_counter = nPoints - 1
break_param = False
error = []

start_time = time.time()

for t in timer:
    if useSCC == False:
        print(t/T)
    #print(profile_counter)
    profile_counter += 1
    p_counter += 1
    q_counter += 1
    rg_profile_counter += 1
    
    if t > 0.15*np.abs(c_star)*T and nucleation_type == 1:
        fit_counter += 1
    elif t > 0.15*T and nucleation_type == 2:
        fit_counter+= 1
    elif t > 0.1*T and nucleation_type == 3:
        fit_counter += 1
    
    if q_counter == q:  
    
        #rgTemp = roll_right(rg_array, Yshape)
        rgTemp = np.roll(rg_array, 1, axis = 1)
        rgTemp[:,-1] = rg_array[:,-1]
        rgTemp[:,0] = rg_array[:,0]
        q_counter = 0
    
    else:
        
        rgTemp = rg_array
    
    if p_counter == p:
        
        #rsTemp = roll_left(rs_array, Yshape) 
        rsTemp = np.roll(rs_array, -1, axis = 1)
        rsTemp[:,-1] = rs_array[:,-1]
        rsTemp[:,0] = rs_array[:,0]
        p_counter = 0
        
    else: 
        
        rsTemp = rs_array
    
    rgNew = rgTemp - fcat*dt*rgTemp + fres*dt*rsTemp
    rsNew = rsTemp + fcat*dt*rgTemp - fres*dt*rsTemp
    
    rg_array = rgNew
    rs_array = rsNew
    
    rg_array[rg_array < 10e-8] = 0
    rs_array[rg_array < 10e-8] = 0
    
        #1 removed calculation of C in for loop and moved it to function    
    if nucleation_type==1:
        c_g=calculate_c_g(c_g,rg_array) 
    elif nucleation_type==2:
        c_g=calculate_c_g(c_g,rg_array)
        D=calculate_D(rg_array+rs_array) ##D seems to want rg/rs after dynamics
    elif nucleation_type==3:
        c_g=calculate_c_g(c_g,rgTemp)
        W=calculate_W(c_g,rgTemp,l_array)
    else: 
        print("unknown nucleation type specified, exiting")
        break
    
    if nucleation_type == 1:
        rg_array[1:,0] = r*dt*c_g[1:]*(1-(c_g[1:]))*(1+B*c_g[1:])
    
    if nucleation_type == 2:
        rg_array[1:,0] = dx*((r*dt*D*(1-D))*(1+B*D))[1:]
    
    if nucleation_type == 3:
        rg_array[1:,0] = (r*dt*W*(1-(W)))[1:]#*(1+B*W)[1:]
        
    ## Boundary Conditions
    #fixed number of length-zero microtubules at origin 
    rg_array[0,0] = leftBoundary
    #rg_array[0,:] = leftBoundary
    #all shrinking microtubules of length 0 die
    rs_array[:,0] = rightBoundary
    
             
    ## Find front position 
    front_pos_idx = np.argmin(np.abs(np.abs(c_g[1:])-10e-6))
        
    front_pos.append(front_pos_idx)
    front_times.append(t)

    front_pos_repeated = []
    front_pos_repeated = [i for i in front_pos if i == front_pos[-1] and i > 0]
         
    if len(front_pos_repeated) > 1000:
        print('Front Velocity Zero, Exiting')
        break_param = True
        break
    
    if fit_counter == nPoints:
        n = 50
        vfront = []
        fit_list = partition(n,front_times,front_pos)
        
        for j in range(n):
            fit = np.polyfit(fit_list[0][j], np.multiply(fit_list[1][j],dx),1)
            vfront.append(fit[0])
        
        v1 = np.average((vfront[-10:-5]))
        v2 = np.average((vfront[-5:]))
        
        if np.abs((v1 - v2)/v1) < 0.01:
            print('Simulation Converged at T = '+str(t))
            print('Relative Error = '+str(np.abs(v1 - v2)/v1))
            break
        else:
            pass
            #print('Relative Error too large = '+str(np.abs(v1 - v2)/v1))
        
        error.append(np.abs(v1 - v2)/v1)
        fit_counter = 0
    
    if fit_counter > nPoints:
        fit_counter = 0
    
    if np.any(np.isnan(c_g)) == True:
        print('Nan encountered')
        break
    
    if np.any(c_g[-20:] > 10e-8) or np.any(rg_array[:,-20]> 10e-8):
        rg_array = np.concatenate([rg_array, np.zeros((rg_array.shape[0],200),dtype=rg_array.dtype)], axis=1)
        rg_array = np.append(rg_array,np.zeros((200,rg_array.shape[1])), axis = 0)
        
        rs_array = np.concatenate([rs_array, np.zeros((rs_array.shape[0],200),dtype=rs_array.dtype)], axis=1)
        rs_array = np.append(rs_array,np.zeros((200,rs_array.shape[1])), axis = 0)
        
        if nucleation_type == 1:
            c_g = np.zeros(rg_array.shape[0])
        elif nucleation_type == 2:
            D = np.zeros(rg_array.shape[0])
            c_g = np.zeros(rg_array.shape[0])
        elif nucleation_type == 3: 
            W = np.zeros(rg_array.shape[0])
            c_g = np.zeros(rg_array.shape[0])
            l_array = np.linspace(0,x,rg_array.shape[0])
            
    if profile_counter == len(timer)//num_saves:
         #print('C Array '+str(save_counter)+' of '+str(num_saves)+' Saved!')
        
         ## Save plus end density profiles
        c_g_array.append(c_g)

        profile_times.append(t)
         #c_s_array[save_counter] = c_s
         #rg_profile[save_counter] = rg_array
         #rs_profile[save_counter] = rs_array
         
        save_counter = save_counter + 1
        profile_counter = 0
         
#print(front_pos)
end_time = time.time()
rg_profile = [rg_profile_1, rg_profile_2, rg_profile_3]
#rs_profile = [rs_profile_1, rs_profile_2, rs_profile_3]
#print('Simulation time: ',(np.abs(end_time - start_time)))
np.savetxt("Time_loop_"+str(idx)+".txt", [(np.abs(end_time - start_time))/3600], fmt = '%f') 
params_sim = {'break_param' : break_param , 'profile_times' : profile_times , 'front_times' : front_times , 'dt' : dt, 'dx' : dx, 'dl' : dl, 'c_g_array' : c_g_array, 'front_pos' : front_pos, 'rg_profile' : rg_profile, 'error' : error}

if useSCC == True:
    params_sim.update(params_dict) 

if useSCC == False:
    params_sim.update(input_params)

outfile = open("Microtubules_"+figure+'_'+str(idx),'wb')
pickle.dump(params_sim, outfile)
outfile.close()