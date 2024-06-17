import numpy as np
import time 
import joblib


def cal_Signal_Concentration(p, N, u, r, K):
    #p, N, u, r: lists of the same size, either usually by genotype
    #K, u float 64,
    # page 4 of wang SI, but that is wrong so fixed here and julia
    # returns caclulated concentration for each population in list
    return_list = []
    for i in range(len(p)):
        a = N[i]**2 * p[i]**2 * r[i]**2
        b = 2 * N[i]**2 * p[i]**2 * r[i]
        c = N[i]**2 * p[i]**2
        d = 2 * N[i] * K * p[i] * r[i] * u
        e = 2 * N[i] * K * p[i] * u
        f = K**2 * u**2
        return_list.append((np.sqrt(a+b+c-d+e+f) + N[i]*p[i] - K*u + N[i]*p[i]*r[i]) / (2*u))
    return np.array(return_list)

# check this one out
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


def fortune_wheel(weights):
    # roulette wheel selection 
    t = time.time_ns()
    accumulation = np.cumsum(weights)
    p = np.random.uniform(0, accumulation[-1])
    for i in range(len(accumulation)):
        if accumulation[i] >= p:
            print((time.time_ns()-t)*10 **-9)
            return i



def eval_genotype_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,auto_pro_Rate,
        pro_Rate,sig_Th,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
        grid_Size,base_Volume,decay_Rate,median_CellDen,K):
    counter = 0
    while counter < size_Pop:
        mix_Num = sample_ztp(lam)

        index_Geno = np.random.choice(size_Pop, mix_Num)
        pool_CellDen = np.full((mix_Num,), env_CellDen[mix_Num])

        p = [pro_Rate[index] for index in index_Geno]
        N = pool_CellDen
        r = [auto_R[index] for index in index_Geno]
        sig_Concentration = cal_Signal_Concentration(p, N, decay_Rate, r, K) / mix_Num

        coop_ON = np.array([int(sig_Concentration[index] > sig_Th[index_Geno[index]]) 
                            for index in range(len(sig_Concentration))])

        for i in index_Geno:
            test_coop_benifit = np.sum(coop_ON) / mix_Num * base_Volume * env_CellDen
            threshold = median_CellDen * base_Volume
            coopPayoff_Pop[i] = coop_Benefit * np.sum([int(indiv > threshold)  for indiv in test_coop_benifit])

            # might throw index error as cell enviorment is low
            high_desncity = [int(env_CellDen[i] > median_CellDen) for i in range(mix_Num)]
            

            coopCost_Pop[i] = coop_Cost * sum(coop_ON)

            mean_sigCon = np.mean(sig_Concentration)
            auto_pro_Rate[i] = pro_Rate[i] * mean_sigCon / (K + mean_sigCon) * auto_R[i]
            sigCost_Pop[i] = sig_Cost*(pro_Rate[i] + auto_pro_Rate[i])

            # the indexing here doesnt really make sense but ok
            fit_Pop[i] = baseline + coopPayoff_Pop[index_Geno[0]]- coopCost_Pop[index_Geno[0]] - sigCost_Pop[index_Geno[0]]

        counter += mix_Num

        return coopPayoff_Pop, coopCost_Pop, auto_pro_Rate, sigCost_Pop, fit_Pop


def parallel_No_auto_faster(size_Pop, mix_Num, pro_Rate, decay_Rate, env_CellDen, sig_Th, median_CellDen, coop_Benefit, coop_Cost, sig_Cost, baseline, index_By_Den):
    new_coopPayoff_Pop = []
    new_coopCost_Pop = []
    new_sigCost_Pop = []
    new_fit_Pop = []
    index_Geno = np.random.choice(size_Pop, mix_Num) # get the indexes
    grid_Size = len(env_CellDen) 
    indexes = index_By_Den[index_Geno].transpose().dot(np.ones(mix_Num)) / mix_Num # H_B_g_j sums from S4 used in S5
    thresholds = sig_Th[index_Geno]
    coop_On = np.full((mix_Num, grid_Size), indexes).transpose() > thresholds # compares the sums to the thresholds to get the full S4
    H_C_sum = coop_On.transpose().dot(np.ones(grid_Size)) # sum used in S3
    H_B_sum = (coop_On.dot(np.ones(mix_Num)) * env_CellDen > np.full((grid_Size,), median_CellDen)).dot(np.ones(grid_Size)) # Sum from S5 and S3 in one
    for i in range(mix_Num):        
        new_coopPayoff_Pop.append(coop_Benefit * H_B_sum)
        new_coopCost_Pop.append(coop_Cost * H_C_sum[i])
        new_sigCost_Pop.append(sig_Cost * pro_Rate[index_Geno[i]])
        new_fit_Pop.append(baseline + coop_Benefit * H_B_sum - coop_Cost * H_C_sum[i] - sig_Cost * pro_Rate[index_Geno[i]])
    return new_coopPayoff_Pop, new_coopCost_Pop, new_sigCost_Pop, new_fit_Pop


def parallel_No_auto(size_Pop, mix_Num, pro_Rate, decay_Rate, env_CellDen, sig_Th, median_CellDen, coop_Benefit, coop_Cost, sig_Cost, baseline):
    new_coopPayoff_Pop = []
    new_coopCost_Pop = []
    new_sigCost_Pop = []
    new_fit_Pop = []
    index_Geno = np.random.choice(size_Pop, mix_Num)
    for geno in index_Geno:
        cost_sum = 0
        benifit_sum = 0

            
        for density in env_CellDen:
            sig_Concentration = sum([pro_Rate[index]/decay_Rate * density / mix_Num for index in index_Geno])
            cost_sum += int(sig_Concentration > sig_Th[geno])
            benifit_sum += int(np.sum([density * int(sig_Concentration > sig_Th[index]) / mix_Num  for index in index_Geno]) > median_CellDen)
    
        new_coopPayoff_Pop.append(coop_Benefit * benifit_sum)
        new_coopCost_Pop.append(coop_Cost * cost_sum)
        new_sigCost_Pop.append(sig_Cost * pro_Rate[geno])
        new_fit_Pop.append(baseline + new_coopPayoff_Pop[-1] - new_coopCost_Pop[-1] - new_sigCost_Pop[-1])
    return new_coopPayoff_Pop, new_coopCost_Pop, new_sigCost_Pop, new_fit_Pop
     

def eval_genotype_No_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,
    pro_Rate,sig_Th,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,
    grid_Size,base_Volume,decay_Rate,median_CellDen):
    counter = 0
    new_coopPayoff_Pop = []
    new_coopCost_Pop = []
    new_sigCost_Pop = []
    new_fit_Pop = []
    mixing_numbers = []
    contribute = pro_Rate / decay_Rate 
    matrix = np.full((grid_Size,size_Pop), contribute)
    index_By_Den = env_CellDen * matrix.transpose()
    while np.sum(mixing_numbers) < size_Pop:
        mixing_numbers.append(sample_ztp(lam))
    # try looping through and doing everythign with numpy maybe? will probably work a bit better.
    everything = joblib.Parallel(n_jobs=8)(joblib.delayed(parallel_No_auto_faster)(
        size_Pop, mix_Num, pro_Rate, decay_Rate, env_CellDen, sig_Th, median_CellDen, coop_Benefit, coop_Cost, sig_Cost, baseline, index_By_Den) for mix_Num in mixing_numbers)

    for item in everything:
        new_coopPayoff_Pop += item[0]
        new_coopCost_Pop += item[1]
        new_sigCost_Pop += item[2]
        new_fit_Pop += item[3]

    return np.array(new_coopPayoff_Pop[0:size_Pop]), np.array(new_coopCost_Pop[0:size_Pop]), np.array(new_sigCost_Pop[0:size_Pop]), np.array(new_fit_Pop[0:size_Pop])

def eval_genotype_No_Auto_No_Probability(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,
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


# testing code:
if __name__ == "__main__":
    test = np.array([1,2,3])
    test2 = np.array([4,5,6,7])
    thresholds = np.array([9,10,11])
    matrix = np.full((4,3), test).transpose()
    # print(matrix)
    matrix2 = np.full((4,4), test2)
    product = matrix.dot(matrix2) / 4
    print(product)
    threshold_matrix = np.full((4,3), thresholds).transpose()
    print(threshold_matrix)
    print(product > threshold_matrix)
    print(((product > threshold_matrix) * test2)> np.full((3,4), 5.5))
