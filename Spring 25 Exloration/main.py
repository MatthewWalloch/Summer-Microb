import numpy as np
import time
from helpers import *
from graph import *
import matplotlib.pyplot as plt
import json
import joblib



def main_QS(lam, testing, Auto=False, max_G=5000, clonal=True):
    ################################################################################
    # set parameters
    ################################################################################
    # population size
    # baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen
    # maximum cellular density (cells per microliter)
    max_CellDen = 10.0 ** 5
    # minimum cellular density
    min_CellDen = 10.0 ** 1.5
    size_Pop = 5000
    grid_Size= 100
    # general parameters
    
    gp = {
        "baseline": 100,
        "coop_Benefit": .75 / 500,
        "coop_Cost":2.5 / 500,
        "sig_Cost": 16 / 50000,
        "size_Pop": size_Pop,
        "lam": lam,
        "env_CellDen": np.array(list(np.linspace(min_CellDen, max_CellDen,num=grid_Size))),
        "grid_Size": grid_Size,
        "k": 50,
        "Ks": 600,
        "X_prod": 5e-9,
        "decay_RateX": 4e-6,
        "decay_RateY": 4e-6,
        "Kx": 500,
        "Y_consumption": 1e-8,
        "XY_rate": 2e0,
        "median_CellDen": (max_CellDen + min_CellDen) * .5,
        "m": "10.0 ** -7.0 to 1e-4 100 steps"
    }
    gp_no_np= {k: gp[k] for k in set(list(gp.keys())) - set(["env_CellDen"])}
    # mutation rate
    mu_Production = 0.01
    mu_DecayRate = 0.01
    mu_induct_Rate = 0.01
    mu_X_pro_Rate = 0.01
    # print("mutation rates increased for testing")
    # mu_Production = 0.1
    # mu_sensitiity = 0.1
    # mu_R = 0.1

    
 
   
    # maximum cellular production rate
    max_ProRate = 2e-08
    # minimum cellular production rate
    min_ProRate = 0
    # initial cellular production rate
    init_pro_Rate = testing
    # SD for mutation of cellular production rate
    mu_SD_ProRate = 0.1e-09
    gp_no_np["Production mutation"] = [max_ProRate, min_ProRate, init_pro_Rate, mu_SD_ProRate]
    # maximum cellular production rate
    max_DecayRate = 10
    # minimum cellular production rate
    min_DecayRate = 10e-11
    # initial cellular production rate
    init_DecayRate = 4e-6
    # SD for mutation of cellular production rate
    mu_SD_DecayRate = 1e-7
    gp_no_np["Decay mutation"] = [max_DecayRate, min_DecayRate, init_DecayRate, mu_SD_DecayRate]
    # maximum ratio of autoinduction production
    max_induct_Rate = 20
    # minimum ratio of autoinduction production
    min_induct_Rate = 0
    # initial ratio of autoinduction production
    init_induct_Rate = 5
    # SD for mutation of ratio of autoinduction production
    mu_SD_induct_Rate = 1e-09
    gp_no_np["Induction mutation"] = [max_induct_Rate, min_induct_Rate, init_induct_Rate, mu_SD_induct_Rate]


    # maximum cellular production rate
    max_X_ProRate = 2e-08
    # minimum cellular production rate
    min_X_ProRate = 0
    # initial cellular production rate
    init_X_pro_Rate = 5e-9
    # SD for mutation of cellular production rate
    mu_SD_X_ProRate = 0.1e-09
    gp_no_np["X Production Mutation"] = [max_X_ProRate, min_X_ProRate, init_X_pro_Rate, mu_SD_X_ProRate]
    ################################################################################
    # initialization
    ################################################################################
    # initialize production rate
    pro_Rate1 = np.full((size_Pop,), init_pro_Rate)
    pro_Rate2 = np.full((size_Pop,), init_pro_Rate)
    # initialize signal threshold
    decay_Rate1 = np.full((size_Pop,), init_DecayRate)
    decay_Rate2 = np.full((size_Pop,), init_DecayRate)
    # initialize ratio of autoinduction production
    induct_Rate1 = np.full((size_Pop,), init_induct_Rate)
    induct_Rate2 = np.full((size_Pop,), init_induct_Rate)

    X_pro_Rate = np.full((size_Pop,), init_X_pro_Rate)

    # initialize genotype fitness
    fit_Pop = np.zeros((size_Pop,))
    # initialize cooperation payoff
    coopPayoff_Pop = np.zeros(size_Pop)
    # initialize cost for signaling
    sigCost_Pop = np.zeros(size_Pop)
    # initialize cost for cooperation
    coopCost_Pop = np.zeros(size_Pop)
    # autoinduction production rate
    auto_pro_Rate = np.zeros(size_Pop)
    # efficiency for cooperation


    # initialize result saving space
    # all lineages
    fit_Evo =  np.zeros(max_G)
    pro_Rate_Evo1 =  np.zeros(max_G)
    pro_Rate_Evo2 =  np.zeros(max_G)
    induct_Rate_Evo1 =  np.zeros(max_G)
    induct_Rate_Evo2 =  np.zeros(max_G)
    decay_Rate_Evo1 =  np.zeros(max_G)
    decay_Rate_Evo2 =  np.zeros(max_G)
    coopPayoff_Evo =  np.zeros(max_G)
    sigCost_Evo =  np.zeros(max_G)
    coopCost_Evo =  np.zeros(max_G)
    coopEff_Evo =  np.zeros(max_G)
    auto_pro_Rate_Evo =  np.zeros(max_G)
    X_pro_Rate_Evo = np.zeros(max_G)

    ################################################################################
    # Evolution
    ################################################################################
    # initial evaluation

    # genotype eval function
    fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop = eval_genotype_Clonal(pro_Rate1, decay_Rate1, induct_Rate1, X_pro_Rate, gp)
    
    
    g = 0
    fit_Evo[g] = np.mean(fit_Pop)
    pro_Rate_Evo1[g] = np.mean(pro_Rate1)
    pro_Rate_Evo2[g] = np.mean(pro_Rate2)
    decay_Rate_Evo1[g] = np.mean(decay_Rate1)
    decay_Rate_Evo2[g] = np.mean(decay_Rate2)
    induct_Rate_Evo1[g] =  np.mean(induct_Rate1)
    induct_Rate_Evo2[g] =  np.mean(induct_Rate2)
    coopPayoff_Evo[g] = np.mean(coopPayoff_Pop)
    sigCost_Evo[g] = np.mean(sigCost_Pop)
    coopCost_Evo[g] = np.mean(coopCost_Pop)
    X_pro_Rate_Evo[g] = np.mean(X_pro_Rate)
    rng = np.random.default_rng()
    t= time.time_ns()
    for g in range(1,max_G):
        gp["m"] = rng.uniform(low=1.5e-7, high=1.5e-4)
        temp_pro_Rate1 = np.zeros(size_Pop)
        temp_pro_Rate2 = np.zeros(size_Pop)
        temp_decay_Rate1 = np.zeros(size_Pop)
        temp_decay_Rate2 = np.zeros(size_Pop)
        temp_induct_Rate1 = np.zeros(size_Pop)
        temp_induct_Rate2 = np.zeros(size_Pop)
        temp_X_pro_Rate = np.zeros(size_Pop)

        temp_fit_Pop = np.zeros(size_Pop)
        temp_coopPayoff_Pop = np.zeros(size_Pop)
        temp_coopCost_Pop = np.zeros(size_Pop)
        temp_sigCost_Pop = np.zeros(size_Pop)
        temp_auto_pro_Rate = np.zeros(size_Pop)
        # print(fit_Evo)
        for i in range(len(fit_Pop)):
            if fit_Pop[i] < 0:
                fit_Pop[i] = 0.0
        fit_pop_probability = fit_Pop
        if sum(fit_pop_probability) == 0:
            fit_pop_probability = np.ones(size_Pop)
        fit_pop_probability = fit_pop_probability / sum(fit_pop_probability)
        index_Select = np.random.choice(size_Pop, size=size_Pop, p=fit_pop_probability)

        for n in range(len(index_Select)):
            choice = int(index_Select[n])
            temp_pro_Rate1[n] = pro_Rate1[choice]
            temp_pro_Rate2[n] = pro_Rate2[choice]
            temp_decay_Rate1[n] = decay_Rate1[choice]
            temp_decay_Rate2[n] = decay_Rate2[choice]
            temp_induct_Rate1[n] = induct_Rate1[choice]
            temp_induct_Rate2[n] = induct_Rate2[choice]
            temp_X_pro_Rate[n] = X_pro_Rate[choice]

            temp_fit_Pop[n] = fit_Pop[choice]
            temp_coopPayoff_Pop[n] = coopPayoff_Pop[choice]
            temp_coopCost_Pop[n] = coopCost_Pop[choice]

            temp_sigCost_Pop[n] = sigCost_Pop[choice]
        
        pro_Rate1 = temp_pro_Rate1
        pro_Rate2 = temp_pro_Rate2
        decay_Rate1 = temp_decay_Rate1
        decay_Rate2 = temp_decay_Rate2
        induct_Rate1 = temp_induct_Rate1
        induct_Rate2 = temp_induct_Rate2
        fit_Pop = temp_fit_Pop
        coopPayoff_Pop = temp_coopPayoff_Pop 
        coopCost_Pop = temp_coopCost_Pop
        sigCost_Pop = temp_sigCost_Pop
        X_pro_Rate = temp_X_pro_Rate

        
        # mutate production rate
        pro_Rate1 = mut_parameter(pro_Rate1,mu_Production,mu_SD_ProRate,min_ProRate,max_ProRate,size_Pop)
        pro_Rate2 = mut_parameter(pro_Rate2,mu_Production,mu_SD_ProRate,min_ProRate,max_ProRate,size_Pop)
        # mutate signal threshold
        # decay_Rate1 = mut_parameter(decay_Rate1,mu_DecayRate,mu_SD_DecayRate,min_DecayRate,max_DecayRate,size_Pop)
        # decay_Rate2 = mut_parameter(decay_Rate2,mu_DecayRate,mu_SD_DecayRate,min_DecayRate,max_DecayRate,size_Pop)
        # # mutate ratio of autoinduction production
        # induct_Rate1 = mut_parameter(induct_Rate1,mu_induct_Rate,mu_SD_induct_Rate,min_induct_Rate,max_induct_Rate,size_Pop)
        # induct_Rate2 = mut_parameter(induct_Rate2,mu_induct_Rate,mu_SD_induct_Rate,min_induct_Rate,max_induct_Rate,size_Pop)
        
        X_pro_Rate = mut_parameter(X_pro_Rate, mu_X_pro_Rate, mu_SD_X_ProRate, min_X_ProRate, max_X_ProRate, size_Pop)
        
        # save results  0 sec ish
        fit_Evo[g] = np.mean(fit_Pop)
        pro_Rate_Evo1[g] = np.mean(pro_Rate1)
        pro_Rate_Evo2[g] = np.mean(pro_Rate2)
        decay_Rate_Evo1[g] = np.mean(decay_Rate1)
        decay_Rate_Evo2[g] = np.mean(decay_Rate2)
        induct_Rate_Evo1[g] =  np.mean(induct_Rate1)
        induct_Rate_Evo2[g] =  np.mean(induct_Rate2)
        coopPayoff_Evo[g] = np.mean(coopPayoff_Pop)
        sigCost_Evo[g] = np.mean(sigCost_Pop)
        coopCost_Evo[g] = np.mean(coopCost_Pop)
        X_pro_Rate_Evo[g] = np.mean(X_pro_Rate)

        # genotype eval function
        fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop = eval_genotype_Clonal(pro_Rate1, decay_Rate1, induct_Rate1, X_pro_Rate, gp)
        
        
        if g % 500 == 499:
            print(f"{testing}: {g+1}    {(time.time_ns()-t)* 10 **-9}")
            t = time.time_ns()
            data = {"fit_Pop": fit_Pop.tolist(),
                    "pro_Rate1": pro_Rate1.tolist(),
                    "pro_Rate2": pro_Rate2.tolist(),
                    "decay_Rate1": decay_Rate1.tolist(),
                    "decay_Rate2": decay_Rate2.tolist(),
                    "induct_Rate1": induct_Rate1.tolist(),
                    "induct_Rate2": induct_Rate2.tolist(),
                    "X_pro_Rate": X_pro_Rate.tolist(),
                    "coopPayoff_Pop": coopPayoff_Pop.tolist(),
                    "sigCost_Pop": sigCost_Pop.tolist(),
                    "coopCost_Pop": coopCost_Pop.tolist(),
                    "params": gp_no_np}
            file = f"Spring 25 Exloration\json\\testing\\gen {g+1}.json"
            with open(file, "w") as f:
                json.dump(data, f,  ensure_ascii=False, indent=4)
    timestr = time.strftime("%m-%d %H-%M-%S")
    file = f"Spring 25 Exloration\json\Production rate testing\\{timestr}.json"
    
    data = {"fit_Pop_Evo": fit_Evo.tolist(),
                "pro_Rate_Evo1": pro_Rate_Evo1.tolist(),
                "pro_Rate_Evo2": pro_Rate_Evo2.tolist(),
                "decay_Rate_Evo1": decay_Rate_Evo1.tolist(),
                "decay_Rate_Evo2": decay_Rate_Evo2.tolist(),
                "induct_Rate_Evo1": induct_Rate_Evo1.tolist(),
                "induct_Rate_Evo2": induct_Rate_Evo2.tolist(),
                "X_pro_Rate_Evo": X_pro_Rate_Evo.tolist(),
                "coopPayoff_Pop_Evo": coopPayoff_Evo.tolist(),
                "sigCost_Pop_Evo": sigCost_Evo.tolist(),
                "coopCost_Pop_Evo": coopCost_Evo.tolist(),
                "params": gp_no_np}
    with open(file, "w") as f:
        json.dump(data, f,  ensure_ascii=False, indent=4)
    return file



if __name__ == "__main__":

    # joblib.Parallel(n_jobs=6)(joblib.delayed(vary_signal)(sig_Cost * 10**8) for sig_Cost in range(5,105,5))
    # joblib.Parallel(n_jobs=6)(joblib.delayed(vary_genotype)(np.round(lam, decimals=1), True) for lam in np.arange(0,10,step=.1))
    # joblib.Parallel(n_jobs=5)(joblib.delayed(main_QS)(1, p) for p in np.linspace(0, 2e-08, 100))
    clonal = False
    # # t = time.time_ns()
    max_G = 5000
    file = main_QS(1, 5e-9, max_G=max_G, clonal=clonal)
    graph(file)
    # # graphNoThresholds.graph(file)
