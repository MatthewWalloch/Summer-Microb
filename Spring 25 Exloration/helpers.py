import numpy as np
import time 
import joblib
import random



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


def eval_genotype_Clonal_two_sig(pro_Rate1, pro_Rate2, decay_Rate1, decay_Rate2, induct_Rate1, induct_Rate2,gp):
    den_Matrix = np.full((gp["size_Pop"], gp["grid_Size"]), gp["env_CellDen"]).transpose()
    production_avg = np.zeros(gp["size_Pop"])
    X_star_avg = np.zeros(gp["size_Pop"])
    Y_star_avg = np.zeros(gp["size_Pop"])
    for m in np.arange(1.5e-7, 1.5e-4, step=100):
        gp["m"] = m
        denom1 = gp["m"]-induct_Rate1*den_Matrix+decay_Rate1
        denom2 = gp["m"]-induct_Rate2*den_Matrix+decay_Rate2
        production1 = (1-np.exp(denom1 * -gp["gen_time"])) *den_Matrix *pro_Rate1 / denom1
        production2 = (1-np.exp(denom2 * -gp["gen_time"])) *den_Matrix *pro_Rate2 / denom2
        total_production = production1 + production2
        X_star = total_production * den_Matrix / (gp["m"]+gp["decay_RateX"])
        Y_star = X_star * gp["XY_rate"]/ (gp["Y_consumption"] * den_Matrix + gp["m"]+gp["decay_RateY"])
        
        production_avg += np.ones(gp["grid_Size"]).dot(total_production)
        X_star_avg += np.ones(gp["grid_Size"]).dot(X_star)
        Y_star_avg += np.ones(gp["grid_Size"]).dot(Y_star)
    signal_cost = production_avg * gp["sig_Cost"]
    cost_sum = X_star_avg * gp["coop_Cost"] 
    benifit_sum = Y_star_avg * gp["coop_Benefit"]
    fitness = gp["baseline"]+benifit_sum-cost_sum-signal_cost

    return fitness, benifit_sum, cost_sum, signal_cost

def eval_genotype_Clonal(pro_Rate1, decay_Rate1, induct_Rate1,gp):
    den_Matrix = np.full((gp["size_Pop"], gp["grid_Size"]), gp["env_CellDen"]).transpose()
    production_avg = np.zeros(gp["size_Pop"])
    X_star_avg = np.zeros(gp["size_Pop"])
    Y_star_avg = np.zeros(gp["size_Pop"])
    for m in np.arange(1.5e-7, 1.5e-4, step=100):
        npNPRku = (den_Matrix * pro_Rate1*(1+induct_Rate1)) - gp["k"]*(decay_Rate1+m)
        contribute = 4*gp["k"]*den_Matrix*pro_Rate1*(decay_Rate1+m) + (-1*npNPRku)**2
        total_production = (npNPRku + np.sqrt(contribute)) / (2*(decay_Rate1+m))

        X_star = gp["X_prod"] * total_production / (total_production + gp["Ks"]) * den_Matrix / (m+gp["decay_RateX"])
        Y_star = X_star / (X_star + gp["Kx"]) * gp["XY_rate"]/ (gp["Y_consumption"] * den_Matrix + m+gp["decay_RateY"])
        
        production_avg += np.ones(gp["grid_Size"]).dot(total_production)
        X_star_avg += np.ones(gp["grid_Size"]).dot(X_star)
        Y_star_avg += np.ones(gp["grid_Size"]).dot(Y_star)
    signal_cost = production_avg * gp["sig_Cost"]
    cost_sum = X_star_avg * gp["coop_Cost"] 
    benifit_sum = Y_star_avg * gp["coop_Benefit"]
    fitness = gp["baseline"]+benifit_sum-cost_sum-signal_cost

    return fitness, benifit_sum, cost_sum, signal_cost
def eval_genotype_No_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,
    pro_Rate,sensitivity,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
    grid_Size,base_Volume,decay_Rate,median_CellDen):
    rng = np.random.default_rng()    
        
    # conbined_rates = np.array(joblib.Parallel(n_jobs= 4)(joblib.delayed(random_parallel)( i, rng, lam, size_Pop, pro_Rate, sensitivity, env_CellDen, grid_Size, median_CellDen, decay_Rate) for i in range(size_Pop)))

    benifit_sum = np.zeros(size_Pop)
    cost_sum = np.zeros(size_Pop)
    signal_cost = np.zeros(size_Pop)
    for i in range(size_Pop):
        mix_Num = sample_ztp(lam)
        indexes = np.array(np.append([i], rng.integers(0, high=size_Pop, size=(mix_Num-1,))))
        # bellow is faster but "less random"
        # indexes = np.array(np.append([i], random.choices(all_index, k=mix_Num-1)), dtype=int)

        conbined_rates = np.average(pro_Rate[indexes]) / decay_Rate
        den_Matrix = np.full((mix_Num, grid_Size), env_CellDen).transpose()

        contribute_matrix = conbined_rates * den_Matrix
        H_C_g_j = (sensitivity[indexes] * contribute_matrix)
        cost_sum[i] = H_C_g_j.transpose().dot(np.ones(grid_Size))[0] * coop_Cost
        H_B_g_j = (env_CellDen * H_C_g_j.transpose()) > np.full((mix_Num,grid_Size), median_CellDen)
        # H_B_g_j = (env_CellDen * H_C_g_j.transpose() - 10**1.5) / (10**5-10**1.5)
        benifit_sum[i] = np.sum((H_B_g_j.transpose().dot(np.ones(mix_Num)) / mix_Num)) * coop_Benefit
        signal_cost[i] = pro_Rate[i] * sig_Cost
    # signal_cost = pro_Rate * sig_Cost
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost




# testing code:
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
