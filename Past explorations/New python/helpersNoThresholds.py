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


def eval_genotype_No_Auto_All(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,
    pro_Rate,sensitivity,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
    grid_Size,base_Volume,decay_Rate,median_CellDen):

    contribute = pro_Rate / decay_Rate 
    contribute_Matrix = np.full((grid_Size,size_Pop), contribute).transpose()
    den_Matrix = np.full((grid_Size, grid_Size), env_CellDen)

    H_C_g_j = sensitivity * (env_CellDen * contribute_Matrix).transpose() 
    H_C_g_j = H_C_g_j.transpose()
    cost_sum =  H_C_g_j.dot(np.ones(grid_Size)) * coop_Cost
    
    H_B_g_j = np.mean(H_C_g_j * env_CellDen, axis=0) > np.full((grid_Size), median_CellDen)
    benifit_sum = np.full((size_Pop,), H_B_g_j.dot(np.ones(grid_Size)) * coop_Benefit)
    signal_cost = pro_Rate * sig_Cost
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost

    return benifit_sum, cost_sum, signal_cost, fitness


def random_parallel(i, rng, lam, size_Pop, pro_Rate, sensitivity, env_CellDen, grid_Size, median_CellDen, decay_Rate):
    mix_Num = sample_ztp(lam)
    indexes = np.array(np.append([i], rng.integers(0, high=size_Pop, size=(mix_Num-1,))))
    # bellow is faster but "less random"
    # indexes = np.array(np.append([i], random.choices(all_index, k=mix_Num-1)), dtype=int)

    conbined_rates = np.average(pro_Rate[indexes]*sensitivity[indexes])
    contribute = conbined_rates / decay_Rate
    H_C_g_j = (env_CellDen * contribute)
    cost_sum =  H_C_g_j.dot(np.ones(grid_Size)) 
    H_B_g_j = (H_C_g_j * env_CellDen) > np.full((grid_Size,), median_CellDen)      
    benifit_sum = H_B_g_j.dot(np.ones(grid_Size)) 
    return benifit_sum, cost_sum


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


    return benifit_sum, cost_sum, signal_cost, fitness


def eval_genotype_Auto_Clonal(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,auto_pro_Rate,
        pro_Rate,sensitivity,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
        grid_Size,base_Volume,decay_Rate,median_CellDen,K):

    den_Matrix = np.full((size_Pop, grid_Size), env_CellDen).transpose()
    npNPRku = (den_Matrix * pro_Rate*(1+auto_R)) - K*decay_Rate
    contribute = 4*K*den_Matrix*pro_Rate*decay_Rate + (-1*npNPRku)**2
    contribute = (npNPRku + contribute ** .5) / (2*decay_Rate)

    H_C_g_j = (sensitivity * contribute ).transpose()
    cost_sum =  H_C_g_j.dot(np.ones(grid_Size)) * coop_Cost
    H_B_g_j = (H_C_g_j * env_CellDen) > np.full((size_Pop,grid_Size), median_CellDen)
    benifit_sum = H_B_g_j.dot(np.ones(grid_Size)) * coop_Benefit
    s_star = contribute.transpose().dot(np.ones(grid_Size)) / grid_Size
    auto_pro_Rate = auto_R * (s_star / (K+s_star)) * pro_Rate
    signal_cost = (pro_Rate + auto_pro_Rate) * sig_Cost 
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost
    return benifit_sum, cost_sum, auto_pro_Rate, signal_cost, fitness


def eval_genotype_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,auto_pro_Rate,
        pro_Rate,sensitivity,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
        grid_Size,base_Volume,decay_Rate,median_CellDen,K):

    rng = np.random.default_rng()
    benifit_sum = np.zeros(size_Pop)
    cost_sum = np.zeros(size_Pop)
    signal_cost = np.zeros(size_Pop)
    for i in range(size_Pop):
        mix_Num = sample_ztp(lam)
        indexes = np.array(np.append([i], rng.integers(0, high=size_Pop, size=(mix_Num-1,))))
        # bellow is faster but "less random"
        # indexes = np.array(np.append([i], random.choices(all_index, k=mix_Num-1)), dtype=int)
        pro_Rate_env = pro_Rate[indexes]
        auto_R_env = auto_R[indexes]
        sensitivity_env = sensitivity[indexes]
        den_Matrix = np.full((mix_Num, grid_Size), env_CellDen).transpose()
        npNPRku = (den_Matrix * pro_Rate_env*(1+auto_R_env)) - K*decay_Rate
        contribute = 4*K*den_Matrix*pro_Rate_env*decay_Rate + (-1*npNPRku)**2
        contribute = (npNPRku + contribute ** .5) / (2*decay_Rate)
        contribute = contribute.dot(np.ones(mix_Num)) / mix_Num
        contribute_matrix = np.full((mix_Num, grid_Size), contribute).transpose()
        H_C_g_j = (sensitivity_env * contribute_matrix)
        cost_sum[i] = H_C_g_j.transpose().dot(np.ones(grid_Size))[0] * coop_Cost
        H_B_g_j = (env_CellDen * H_C_g_j.transpose()) > np.full((mix_Num,grid_Size), median_CellDen)
        benifit_sum[i] = np.sum((H_B_g_j.transpose().dot(np.ones(mix_Num)) / mix_Num)) * coop_Benefit
        s_star = contribute.transpose().dot(np.ones(grid_Size)) / grid_Size 
        auto_pro_Rate[i] = auto_R[i] * (s_star / (K+s_star)) * pro_Rate[i]
        signal_cost[i] = (pro_Rate[i] + auto_pro_Rate[i]) * sig_Cost 

    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost
   
    return benifit_sum, cost_sum, auto_pro_Rate, signal_cost, fitness

def eval_genotype_No_Auto_Psuedo_Rand(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,
    pro_Rate,sensitivity,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
    grid_Size,base_Volume,decay_Rate,median_CellDen, randlist):   
        
    # conbined_rates = np.array(joblib.Parallel(n_jobs= 4)(joblib.delayed(random_parallel)( i, rng, lam, size_Pop, pro_Rate, sensitivity, env_CellDen, grid_Size, median_CellDen, decay_Rate) for i in range(size_Pop)))

    benifit_sum = np.zeros(size_Pop)
    cost_sum = np.zeros(size_Pop)
    signal_cost = np.zeros(size_Pop)
    for i in range(size_Pop):
        indexes = np.array(np.append([i], randlist[i]))
        # bellow is faster but "less random"
        # indexes = np.array(np.append([i], random.choices(all_index, k=mix_Num-1)), dtype=int)
    
        conbined_rates = np.average(pro_Rate[indexes]) / decay_Rate
        den_Matrix = np.full((lam, grid_Size), env_CellDen).transpose()

        contribute_matrix = conbined_rates * den_Matrix
        H_C_g_j = (sensitivity[indexes] * contribute_matrix)
        cost_sum[i] = H_C_g_j.transpose().dot(np.ones(grid_Size))[0] * coop_Cost
        H_B_g_j = (env_CellDen * H_C_g_j.transpose()) > np.full((lam,grid_Size), median_CellDen)
        # H_B_g_j = (env_CellDen * H_C_g_j.transpose() - 10**1.5) / (10**5-10**1.5)
        benifit_sum[i] = np.sum((H_B_g_j.transpose().dot(np.ones(lam)) / lam)) * coop_Benefit
        signal_cost[i] = pro_Rate[i] * sig_Cost
    # signal_cost = pro_Rate * sig_Cost
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost


    return benifit_sum, cost_sum, signal_cost, fitness


def eval_genotype_No_Auto_Clonal_No_Den(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,
    pro_Rate,sensitivity,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
    grid_Size,base_Volume,decay_Rate,median_CellDen):

    contribute = pro_Rate / decay_Rate 
    contribute_Matrix = np.full((grid_Size,size_Pop), contribute).transpose()
    den_Matrix = np.full((grid_Size, grid_Size), env_CellDen)

    H_C_g_j = sensitivity * (env_CellDen * contribute_Matrix).transpose() 
    H_C_g_j = H_C_g_j.transpose()
    cost_sum =  H_C_g_j.dot(np.ones(grid_Size)) * coop_Cost
    H_B_g_j = .5 * np.log((H_C_g_j * env_CellDen) / median_CellDen + 1)
    benifit_sum = H_B_g_j.dot(np.ones(grid_Size)) * coop_Benefit
    signal_cost = pro_Rate * sig_Cost
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost

    return benifit_sum, cost_sum, signal_cost, fitness

def eval_genotype_No_Auto_Psuedo_Rand_No_Den(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,
    pro_Rate,sensitivity,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
    grid_Size,base_Volume,decay_Rate,median_CellDen, randlist):   
        
    # conbined_rates = np.array(joblib.Parallel(n_jobs= 4)(joblib.delayed(random_parallel)( i, rng, lam, size_Pop, pro_Rate, sensitivity, env_CellDen, grid_Size, median_CellDen, decay_Rate) for i in range(size_Pop)))

    benifit_sum = np.zeros(size_Pop)
    cost_sum = np.zeros(size_Pop)
    signal_cost = np.zeros(size_Pop)
    for i in range(size_Pop):
        indexes = np.array(np.append([i], randlist[i]))
        
        conbined_rates = np.average(pro_Rate[indexes]) / decay_Rate
        den_Matrix = np.full((lam, grid_Size), env_CellDen).transpose()

        contribute_matrix = conbined_rates * den_Matrix
        H_C_g_j = (sensitivity[indexes] * contribute_matrix)
        cost_sum[i] = H_C_g_j.transpose().dot(np.ones(grid_Size))[0] * coop_Cost

        #.5 * np.log((H_C_g_j * env_CellDen) / median_CellDen + 1)
        H_B_g_j =np.log( (env_CellDen * H_C_g_j.transpose())/ median_CellDen + 1)
        # print(1 *H_B_g_j)
        # H_B_g_j = (env_CellDen * H_C_g_j.transpose() - 10**1.5) / (10**5-10**1.5)
        benifit_sum[i] = np.sum((H_B_g_j.transpose().dot(np.ones(lam)) / lam)) * coop_Benefit * sensitivity[i]/ np.average(sensitivity[indexes])
        signal_cost[i] = pro_Rate[i] * sig_Cost
    # signal_cost = pro_Rate * sig_Cost
    fitness = np.full((size_Pop,), baseline) + benifit_sum - cost_sum - signal_cost


    return benifit_sum, cost_sum, signal_cost, fitness


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
