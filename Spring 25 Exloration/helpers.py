import numpy as np
import time 
import joblib
import random
from numba import jit


def mut_parameter(mut_vector, mut_P, mut_SD, mut_Min, mut_Max,size_Pop):
    # mutation operation for evolving traits number chosen based on poission on mut_P and mutated 
    # based on truncated normal based on mut_SD, mut_Min, mut_Max with mean of unmutated value
    # mut_vector, index_Cheats is list of (size_Pop,) 
    # mut_P, mut_SD, mut_Min, mut_Max are Floats
    # size_Pop is Int
    # returns new mut_vector of size (size_Pop,)
    num_Mut = np.random.poisson(mut_P*size_Pop)
    if num_Mut !=0:
        index_Mut = np.random.choice(size_Pop, size=num_Mut, replace=False)
        for i in index_Mut:
            mutation =  np.random.normal(mut_vector[i], mut_SD)
            if mutation < mut_Min:
                mutation = mut_Min
            if mutation > mut_Max:
                mutation = mut_Max
            mut_vector[i] = mutation
    return mut_vector


def sample_ztp(lam):
    # zero-truncated Poisson distribution based on lambda, returns int
    u = np.random.uniform(np.exp(-lam), 1)
    t = -np.log(u)
    return 1 + np.random.poisson(lam - t)


#@jit tag is for numba
#https://numba.pydata.org/numba-doc/dev/user/jit.html
# basically pre-compiles function to speed up calculations
@jit(cache=True, parallel=True, fastmath=True)
def eval_genotype_Clonal(pro_Rate1, decay_Rate1, induct_Rate1, X_pro_Rate, den_Matrix, gp):
    
    sizepop = int(gp["size_Pop"])
    grid_Size = int(gp["grid_Size"])
    production = np.zeros(sizepop)
    X_star_avg = np.zeros(sizepop)
    Y_star_avg = np.zeros(sizepop)
    #we chunck up m to speed up claculation
    # we usually replace the sums with vector caclculations
    for m in np.linspace(1.5e-7, 1.5e-4, num=100):
        npNPRku = (den_Matrix * pro_Rate1*(1+induct_Rate1)) - gp["k"]*(decay_Rate1+m)
        contribute = 4*gp["k"]*den_Matrix*pro_Rate1*(decay_Rate1+m) + (-1*npNPRku)**2
        total_production = (npNPRku + np.sqrt(contribute)) / (2*(decay_Rate1+m))

        X_star = X_pro_Rate * total_production / (total_production + gp["Ks"]) * den_Matrix / (m+gp["decay_RateX"])
        Y_star = X_star / (X_star + gp["Kx"]) * gp["XY_rate"]/ (gp["Y_consumption"] * den_Matrix + m+gp["decay_RateY"])
        
        production += np.ones(grid_Size).dot(total_production)
        X_star_avg += np.ones(grid_Size).dot(X_star)
        Y_star_avg += np.ones(grid_Size).dot(Y_star)
    signal_cost = production * gp["sig_Cost"]
    cost_sum = X_star_avg * gp["coop_Cost"] 
    benifit_sum = Y_star_avg * gp["coop_Benefit"]
    fitness = gp["baseline"]+benifit_sum-cost_sum-signal_cost
    
    return fitness, benifit_sum, cost_sum, signal_cost


# Baseline function for two signal, however you should double check all the lines to make sure it is right

#  def eval_genotype_Clonal_two_sig(pro_Rate1, pro_Rate2, decay_Rate1, decay_Rate2, induct_Rate1, induct_Rate2,gp):

#     ## this is wrong do not use

#     den_Matrix = np.full((gp["size_Pop"], gp["grid_Size"]), gp["env_CellDen"]).transpose()
#     m_Matrix = np.full((gp["size_Pop"], gp["grid_Size"]), gp["env_CellDen"]).transpose()
#     production_avg = np.zeros(gp["size_Pop"])
#     X_star_avg = np.zeros(gp["size_Pop"])
#     Y_star_avg = np.zeros(gp["size_Pop"])
    
#     for m in np.linspace(1.5e-7, 1.5e-4, num=100):
#         gp["m"] = m
#         denom1 = gp["m"]-induct_Rate1*den_Matrix+decay_Rate1
#         denom2 = gp["m"]-induct_Rate2*den_Matrix+decay_Rate2
#         production1 = (1-np.exp(denom1 * -gp["gen_time"])) *den_Matrix *pro_Rate1 / denom1
#         production2 = (1-np.exp(denom2 * -gp["gen_time"])) *den_Matrix *pro_Rate2 / denom2
#         total_production = production1 + production2
#         X_star = total_production * den_Matrix / (gp["m"]+gp["decay_RateX"])
#         Y_star = X_star * gp["XY_rate"]/ (gp["Y_consumption"] * den_Matrix + gp["m"]+gp["decay_RateY"])
        
#         production_avg += np.ones(gp["grid_Size"]).dot(total_production)
#         X_star_avg += np.ones(gp["grid_Size"]).dot(X_star)
#         Y_star_avg += np.ones(gp["grid_Size"]).dot(Y_star)
#     signal_cost = production_avg * gp["sig_Cost"]
#     cost_sum = X_star_avg * gp["coop_Cost"] 
#     benifit_sum = Y_star_avg * gp["coop_Benefit"]
#     fitness = gp["baseline"]+benifit_sum-cost_sum-signal_cost

#     return fitness, benifit_sum, cost_sum, signal_cost
   



# testing code may or may not be helpful at all
if __name__ == "__main__":
    rng = np.random.default_rng()
    print(rng.integers(0, high=10))
    print(time.strftime("%d-%m %H-%M-%S"))
    # test = np.array([1,2,3])
    # test2 = np.array([4,5,6,7])
    # print(test / [1,2,3])
    # thresholds = np.array([[9,8],[10,0],[11,0]])
    # matrix = np.full((4,3), test).transpose()
    # # print(matrix)
    # matrix2 = np.full((4,4), test2)
    # product = matrix.dot(matrix2) / 4
    # print(product)
    # threshold_matrix = np.full((4,3), thresholds).transpose()
    # print(threshold_matrix)
    # print(product > threshold_matrix)
    # print(((product > threshold_matrix) * test2)> np.full((3,4), 5.5))
    cube = np.full((2,4, 3), [1,2,3])
    cube2 = np.full((3,2, 3), [1,2,3])
    print(cube.sum(np.full((2,4,1),1)))
