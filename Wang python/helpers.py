import numpy as np
import time 
import joblib
import random


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


def eval_genotype_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,auto_pro_Rate,
        pro_Rate,sig_Th,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
        grid_Size,base_Volume,decay_Rate,median_CellDen,K):

    #auto production rate might be fucking up somehwere?? not sure but figure out.
    rng = np.random.default_rng()
    all_index = range(size_Pop)
    mixing_numbers = []
    conbined_rates = []
    max_selection = 0
    total = 0
    den_Matrix = np.full((size_Pop, grid_Size), env_CellDen).transpose()
    npNPRku = (den_Matrix * pro_Rate*(1+auto_R)) - K*decay_Rate
    contribute = 4*K*den_Matrix*pro_Rate*decay_Rate + (-1*npNPRku)**2
    contribute = (npNPRku + contribute ** .5) / (2*decay_Rate)
    thresholds = np.full((size_Pop, 40), np.full((40,), np.inf))
    for i in range(size_Pop):
        mix_Num = sample_ztp(lam)
        if mix_Num > max_selection:
            max_selection = mix_Num
        mixing_numbers.append(1 / mix_Num)
        indexes = np.array(np.append([i], rng.integers(0, high=size_Pop, size=(mix_Num-1,))))
        # bellow is faster but "less random"
        # indexes = np.array(np.append([i], random.choices(all_index, k=mix_Num-1)), dtype=int)
        conbined_rates.append(np.mean(contribute.transpose()[indexes], axis=0))
        thresholds[i][range(mix_Num)] = sig_Th[indexes]
    thresholds = thresholds[:,:max_selection] 

    # work out matricies and broadcasting here, try mathematica notebook?
    threshold_matrix = np.full((grid_Size, size_Pop), thresholds[:,0])
    H_C_g_j = conbined_rates > threshold_matrix.transpose()
    cost_sum =  H_C_g_j.dot(np.ones(grid_Size)) * coop_Cost
    H_B_g_j = (H_C_g_j * env_CellDen)
    for i in range(1, max_selection):
        threshold_matrix = np.full((grid_Size, size_Pop), thresholds[:,i]).transpose()
        H_C_g_j = (conbined_rates ) > threshold_matrix
        H_B_g_j += (H_C_g_j * env_CellDen) 
    H_B_g_j = (mixing_numbers * H_B_g_j.transpose()).transpose() > np.full((size_Pop,grid_Size), median_CellDen)      
    benifit_sum = H_B_g_j.dot(np.ones(grid_Size)) * coop_Benefit
    # test this is the right average bellow, check against paper
    s_star = contribute.transpose().dot(np.ones(grid_Size)) / grid_Size
    auto_pro_Rate = auto_R * (s_star / (K+s_star)) * pro_Rate
    signal_cost = (pro_Rate + auto_pro_Rate) * sig_Cost 
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost - sig_Cost * auto_pro_Rate
    return benifit_sum, cost_sum, auto_pro_Rate, signal_cost, fitness
    

def eval_genotype_No_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,
    pro_Rate,sig_Th,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
    grid_Size,base_Volume,decay_Rate,median_CellDen):
    rng = np.random.default_rng()
    all_index = range(size_Pop)
    mixing_numbers = []
    max_selection = 0
    total = 0
    conbined_rates = np.zeros(size_Pop)
    thresholds = np.full((size_Pop, 40), np.full((40,), np.inf))
    
    for i in range(size_Pop):
        mix_Num = sample_ztp(lam)
        if mix_Num > max_selection:
            max_selection = mix_Num
        mixing_numbers.append(1 / mix_Num)
        indexes = np.array(np.append([i], rng.integers(0, high=size_Pop, size=(mix_Num-1,))))
        # bellow is faster but "less random"
        # indexes = np.array(np.append([i], random.choices(all_index, k=mix_Num-1)), dtype=int)
        conbined_rates[i] = np.average(pro_Rate[indexes])
        thresholds[i][range(mix_Num)] = sig_Th[indexes]
    thresholds = thresholds[:,:max_selection] 
   
    contribute = conbined_rates / decay_Rate 
    contribute_Matrix = np.full((grid_Size,size_Pop), contribute).transpose()
    den_Matrix = np.full((grid_Size, grid_Size), env_CellDen)
    threshold_matrix = np.full((grid_Size, size_Pop), thresholds[:,0]).transpose()
    H_C_g_j = (env_CellDen * contribute_Matrix ) > threshold_matrix
    cost_sum =  H_C_g_j.dot(np.ones(grid_Size)) * coop_Cost
    H_B_g_j = (H_C_g_j * env_CellDen)
    for i in range(1, max_selection):
        threshold_matrix = np.full((grid_Size, size_Pop), thresholds[:,i]).transpose()
        H_C_g_j = (env_CellDen * contribute_Matrix ) > threshold_matrix
        H_B_g_j += (H_C_g_j * env_CellDen) 
    H_B_g_j = (mixing_numbers * H_B_g_j.transpose()).transpose() > np.full((size_Pop,grid_Size), median_CellDen)      
    benifit_sum = H_B_g_j.dot(np.ones(grid_Size)) * coop_Benefit
    signal_cost = pro_Rate * sig_Cost
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost
    return benifit_sum, cost_sum, signal_cost, fitness


def eval_genotype_No_Auto_Clonal(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,
    pro_Rate,sig_Th,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
    grid_Size,base_Volume,decay_Rate,median_CellDen):

    contribute = pro_Rate / decay_Rate 
    contribute_Matrix = np.full((grid_Size,size_Pop), contribute).transpose()
    den_Matrix = np.full((grid_Size, grid_Size), env_CellDen)
    threshold_matrix = np.full((grid_Size, size_Pop), sig_Th).transpose()
    H_C_g_j = (env_CellDen * contribute_Matrix ) > threshold_matrix
    cost_sum =  H_C_g_j.dot(np.ones(grid_Size)) * coop_Cost
    H_B_g_j = (H_C_g_j * env_CellDen) > np.full((size_Pop,grid_Size), median_CellDen)
    benifit_sum = H_B_g_j.dot(np.ones(grid_Size)) * coop_Benefit
    signal_cost = pro_Rate * sig_Cost
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost

    return benifit_sum, cost_sum, signal_cost, fitness


def eval_genotype_Auto_Clonal(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,auto_pro_Rate,
        pro_Rate,sig_Th,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
        grid_Size,base_Volume,decay_Rate,median_CellDen,K):

    den_Matrix = np.full((size_Pop, grid_Size), env_CellDen).transpose()
    npNPRku = (den_Matrix * pro_Rate*(1+auto_R)) - K*decay_Rate
    contribute = 4*K*den_Matrix*pro_Rate*decay_Rate + (-1*npNPRku)**2
    contribute = (npNPRku + contribute ** .5) / (2*decay_Rate)
    
    threshold_matrix = np.full((grid_Size, size_Pop), sig_Th)
    H_C_g_j = (contribute > threshold_matrix).transpose()
    cost_sum =  H_C_g_j.dot(np.ones(grid_Size)) * coop_Cost
    H_B_g_j = (H_C_g_j * env_CellDen) > np.full((size_Pop,grid_Size), median_CellDen)
    benifit_sum = H_B_g_j.dot(np.ones(grid_Size)) * coop_Benefit
    s_star = contribute.transpose().dot(np.ones(grid_Size)) / grid_Size
    auto_pro_Rate = auto_R * (s_star / (K+s_star)) * pro_Rate
    signal_cost = (pro_Rate + auto_pro_Rate) * sig_Cost 
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost

    return benifit_sum, cost_sum, auto_pro_Rate, signal_cost, fitness


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
