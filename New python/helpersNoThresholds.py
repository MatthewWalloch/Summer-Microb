import numpy as np
import time 
import joblib
import random
import fastrand


def mut_parameter(mut_vector, mut_P, mut_SD, mut_Min, mut_Max, index_Cheats,size_Pop):
    # mutation operation for evolving traits number chosen based on poission on mut_P and mutated 
    # based on truncated normal based on mut_SD, mut_Min, mut_Max with mean of unmutated value
    # mut_vector, index_Cheats is list of (size_Pop,) 
    # mut_P, mut_SD, mut_Min, mut_Max are Floats
    # size_Pop is Int
    # returns new mut_vector of size (size_Pop,)
    num_Mut = np.random.poisson(mut_P*size_Pop)
    if num_Mut !=0:
        index_Mut = np.random.choice(size_Pop, size=num_Mut, replace=False)
        index_Diff = [i for i in index_Mut if index_Cheats[i] != 1]
        for i in index_Diff:
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


def eval_genotype_No_Auto_Clonal(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,
    pro_Rate,sensitivity,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
    grid_Size,base_Volume,decay_Rate,median_CellDen):

    contribute = pro_Rate / decay_Rate 
    contribute_Matrix = np.full((grid_Size,size_Pop), contribute).transpose()
    den_Matrix = np.full((grid_Size, grid_Size), env_CellDen)

    H_C_g_j = sensitivity * (env_CellDen * contribute_Matrix).transpose() 
    H_C_g_j = H_C_g_j.transpose()
    cost_sum =  H_C_g_j.dot(np.ones(grid_Size)) * coop_Cost
    H_B_g_j = (H_C_g_j * env_CellDen) > np.full((size_Pop,grid_Size), median_CellDen)
    benifit_sum = H_B_g_j.dot(np.ones(grid_Size)) * coop_Benefit
    signal_cost = pro_Rate * sig_Cost
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost

    return benifit_sum, cost_sum, signal_cost, fitness





# testing code:
if __name__ == "__main__":
    rng = np.random.default_rng()
    print(rng.integers(0, high=10))
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
