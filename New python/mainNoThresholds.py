import numpy as np
import time
from helpersNoThresholds import *
import matplotlib.pyplot as plt
import json
import joblib
import graphNoThresholds

def main_QS(init_SN, sig_Cost ,lam , K, mu_Cheats, Auto=False, max_G=500, clonal=True):
    ################################################################################
    # QS with autoregulation
    ################################################################################
    # coop_Cost: constant for cost of cooperation, Float
    # sig_Cost: constant for cost of signaling, Float
    # lambda: parameter used in the zero-truncated Poisson distribution for
    #         generating number of mixing genotypes, Float
    # K: half concentration value, Float
    # mu_Cheats: parameter used in the zero-truncated Poisson distribution for
    #            generating number of cheaters, FLoat
    # rep_Num: index number of replications, Int




    ################################################################################
    # set parameters
    ################################################################################
    # population size
    size_Pop = 5000
    # environments
    grid_Size = 100
    # baseline fitness
    baseline = 100
    # payoff for gene turned 'ON' & cooperation
    coop_Benefit = 1.5
    coop_Cost = .5

    # mutation rate

    # mu_Production = 0.01
    # mu_sensitiity = 0.01
    # mu_R = 0.01
    # print("mutation rates increased for testing")
    mu_Production = 0.1
    mu_sensitiity = 0.1
    mu_R = 0.1

    # maximum cellular density (cells per microliter)
    max_CellDen = 10.0 ** 5
    # minimum cellular density
    min_CellDen = 10.0 ** 1.5
    # median cellular density
    median_CellDen = (max_CellDen + min_CellDen) * .5
    # baseline volume
    base_Volume = 10.0
    # signal decay rate
    decay_Rate = 10.0 ** -4.0

    # initial testing environments
    env_CellDen = np.array(list(np.linspace(min_CellDen, max_CellDen,num=grid_Size)))
    # maximum cellular production rate
    max_ProRate = 2e-08
    # minimum cellular production rate
    min_ProRate = 0
    # initial cellular production rate
    init_pro_Rate = 0.5e-08 
    # SD for mutation of cellular production rate
    mu_SD_ProRate = 0.01e-08

    # maximum signal concentration
    max_sensitivity = 20
    # minimum signal concentration
    min_sensitivity = 0.00000001
    # initial signal threshold
    # init_SN = .01
    # SD for mutation of signal concentration
    mu_SD_sensitivity = .01

    # maximum ratio of autoinduction production
    max_R = 20.
    # minimum ratio of autoinduction production
    min_R = 0.01
    # initial ratio of autoinduction production
    init_R = 5.
    # SD for mutation of ratio of autoinduction production
    mu_SD_R = 1.0

    # define cheats parameters
    cheats_proRate = 0.
    cheats_sigTh = 2e1

    ################################################################################
    # initialization
    ################################################################################
    # initialize production rate
    pro_Rate = np.full((size_Pop,), init_pro_Rate)
    # initialize signal threshold
    sensitivity = np.full((size_Pop,), init_SN)
    # initialize ratio of autoinduction production
    auto_R = np.full((size_Pop,), init_R)
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
    # cheats index
    index_Cheats = np.zeros(size_Pop)
    # selection index
    index_Select = np.zeros(size_Pop)

    # initialize result saving space
    # all lineages
    fit_Evo =  np.zeros(max_G)
    pro_Rate_Evo =  np.zeros(max_G)
    sensitivity_Evo =  np.zeros(max_G)
    auto_R_Evo =  np.zeros(max_G)
    coopPayoff_Evo =  np.zeros(max_G)
    sigCost_Evo =  np.zeros(max_G)
    coopCost_Evo =  np.zeros(max_G)
    coopEff_Evo =  np.zeros(max_G)
    auto_pro_Rate_Evo =  np.zeros(max_G)

    # record cheats
    numCheats_Evo =  np.zeros(max_G)
    ################################################################################
    # Evolution
    ################################################################################
    # initial evaluation
    if Auto:
        if clonal:
            coopPayoff_Pop, coopCost_Pop, auto_pro_Rate, sigCost_Pop, fit_Pop = eval_genotype_Auto_Clonal(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,auto_pro_Rate,pro_Rate,sensitivity,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen,K)
        else:
            coopPayoff_Pop, coopCost_Pop, auto_pro_Rate, sigCost_Pop, fit_Pop = eval_genotype_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,auto_pro_Rate,pro_Rate,sensitivity,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen,K)  
    else:
        if clonal:
            coopPayoff_Pop, coopCost_Pop, sigCost_Pop, fit_Pop = eval_genotype_No_Auto_Clonal(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,pro_Rate,sensitivity,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen)
        else:
            coopPayoff_Pop, coopCost_Pop, sigCost_Pop, fit_Pop = eval_genotype_No_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,pro_Rate,sensitivity,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen)  
   
    g = 0
    numCheats_Evo[g] = np.sum(index_Cheats)
    fit_Evo[g] = np.mean(fit_Pop)
    pro_Rate_Evo[g] = np.mean(pro_Rate)
    sensitivity_Evo[g] = np.mean(sensitivity)
    auto_R_Evo[g] = np.mean(auto_R)
    coopPayoff_Evo[g] = np.mean(coopPayoff_Pop)
    sigCost_Evo[g] = np.mean(sigCost_Pop)
    coopCost_Evo[g] = np.mean(coopCost_Pop)
    auto_pro_Rate_Evo[g] = np.mean(auto_pro_Rate)

    t= time.time_ns()
    for g in range(1,max_G):

        temp_pro_Rate = np.zeros(size_Pop)
        temp_sensitivity = np.zeros(size_Pop)
        temp_fit_Pop = np.zeros(size_Pop)
        temp_coopPayoff_Pop = np.zeros(size_Pop)
        temp_coopCost_Pop = np.zeros(size_Pop)
        temp_sigCost_Pop = np.zeros(size_Pop)
        temp_index_Cheats = np.zeros(size_Pop)
        temp_auto_R = np.zeros(size_Pop)
        temp_auto_pro_Rate = np.zeros(size_Pop)
        # print(fit_Evo)
        for i in range(len(fit_Pop)):
            if fit_Pop[i] < 0:
                fit_Pop[i] = 0
        fit_pop_probability = fit_Pop
        fit_pop_probability = fit_pop_probability / sum(fit_pop_probability)
        index_Select = np.random.choice(size_Pop, size=size_Pop, p=fit_pop_probability)

        for n in range(len(index_Select)):
            temp_pro_Rate[n] = pro_Rate[int(index_Select[n])]
            temp_sensitivity[n] = sensitivity[int(index_Select[n])]
            temp_fit_Pop[n] = fit_Pop[int(index_Select[n])]
            temp_coopPayoff_Pop[n] = coopPayoff_Pop[int(index_Select[n])]
            temp_coopCost_Pop[n] = coopCost_Pop[int(index_Select[n])]

            temp_sigCost_Pop[n] = sigCost_Pop[int(index_Select[n])]
            temp_index_Cheats[n] = index_Cheats[int(index_Select[n])]
            temp_auto_R[n] = auto_R[int(index_Select[n])]
            temp_auto_pro_Rate[n] = auto_pro_Rate[int(index_Select[n])]
        
        pro_Rate = temp_pro_Rate
        sensitivity = temp_sensitivity
        fit_Pop = temp_fit_Pop
        coopPayoff_Pop = temp_coopPayoff_Pop 
        coopCost_Pop = temp_coopCost_Pop
        igCost_Pop = temp_sigCost_Pop
        index_Cheats = temp_index_Cheats
        auto_pro_Rate = temp_auto_pro_Rate
        auto_R = temp_auto_R
        
        # # generate cheats -- 0 sec
        # if numCheats_Evo[g] < size_Pop:
        #     num_Cheats = np.random.poisson(lam=mu_Cheats)
        #     if num_Cheats!=0:
        #         temp_Index = np.random.choice(size_Pop,size=num_Cheats, replace=False)
        #         for n in temp_Index:
        #             index_Cheats[n] = 1
        #             pro_Rate[n] = cheats_proRate
        #             sig_Th[n] = cheats_sigTh
        
        
        # while section <<.005
        # mutate production rate
        
        pro_Rate = mut_parameter(pro_Rate,mu_Production,mu_SD_ProRate,min_ProRate,max_ProRate,index_Cheats,size_Pop)

        # mutate signal threshold
        sensitivity = mut_parameter(sensitivity,mu_sensitiity,mu_SD_sensitivity,min_sensitivity,max_sensitivity,index_Cheats,size_Pop)

        # mutate ratio of autoinduction production
        auto_R = mut_parameter(auto_R,mu_R,mu_SD_R,min_R,max_R,index_Cheats,size_Pop)
        # record cheats
        numCheats_Evo[g] = np.sum(index_Cheats)
        
        # save results  0 sec ish
        fit_Evo[g] = np.mean(fit_Pop)
        pro_Rate_Evo[g] = np.mean(pro_Rate)
        sensitivity_Evo[g] = np.mean(sensitivity)
        auto_R_Evo[g] = np.mean(auto_R)
        coopPayoff_Evo[g] = np.mean(coopPayoff_Pop)
        sigCost_Evo[g] = np.mean(sigCost_Pop)
        coopCost_Evo[g] = np.mean(coopCost_Pop)
        auto_pro_Rate_Evo[g] = np.mean(auto_pro_Rate)

        
        
        # genotype evaluation
        
        if Auto:
            if clonal:
                coopPayoff_Pop, coopCost_Pop, auto_pro_Rate, sigCost_Pop, fit_Pop = eval_genotype_Auto_Clonal(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,auto_pro_Rate,pro_Rate,sensitivity,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen,K)
            else:
                coopPayoff_Pop, coopCost_Pop, auto_pro_Rate, sigCost_Pop, fit_Pop = eval_genotype_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,auto_pro_Rate,pro_Rate,sensitivity,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen,K)  
        else:
            if clonal:
                coopPayoff_Pop, coopCost_Pop, sigCost_Pop, fit_Pop = eval_genotype_No_Auto_Clonal(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,pro_Rate,sensitivity,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen)
            else:
                coopPayoff_Pop, coopCost_Pop, sigCost_Pop, fit_Pop = eval_genotype_No_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,sigCost_Pop,pro_Rate,sensitivity,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen)  

        if g % 500 == 499:
            print(f"{lam}: {g+1}    {(time.time_ns()-t)* 10 **-9}")
            t = time.time_ns()
            data = {"fit_Pop": fit_Pop.tolist(),
                "pro_Rate": pro_Rate.tolist(),
                "sensitivity": sensitivity.tolist(),
                "auto_R": auto_R.tolist(),
                "coopPayoff_Pop": coopPayoff_Pop.tolist(),
                "sigCost_Pop": sigCost_Pop.tolist(),
                "coopCost_Pop": coopCost_Pop.tolist(),
                "auto_pro_Rate": auto_pro_Rate.tolist()}
            tgen = time.asctime().replace(":", "-" )
            if Auto:
                file = f"New python\\auto json\\generations\gen {g+1} {lam} {tgen} {init_SN} {sig_Cost} {max_G} {clonal}.json"
            else:
                file = f"New python\\json\\50 genotypes\gen {g+1} {lam} {tgen} {init_SN} {sig_Cost} {max_G} {clonal}.json"
            with open(file, "w") as f:
                json.dump(data, f,  ensure_ascii=False, indent=4)
            # graphNoThresholds.graph_last_gen(file)


        # if g == 1000:
        #     data = {"fit_Pop": fit_Pop.tolist(),
        #         "pro_Rate": pro_Rate.tolist(),
        #         "sensitivity": sensitivity.tolist(),
        #         "auto_R": auto_R.tolist(),
        #         "coopPayoff_Pop": coopPayoff_Pop.tolist(),
        #         "sigCost_Pop": sigCost_Pop.tolist(),
        #         "coopCost_Pop": coopCost_Pop.tolist(),
        #         "auto_pro_Rate": auto_pro_Rate.tolist()}
        #     tgen = time.asctime().replace(":", "-" )
        #     with open(f"New python\json\\testing results\gen {g} {tgen} {init_SN} {sig_Cost} {max_G} {clonal}.json", "w") as f:
        #         json.dump(data, f,  ensure_ascii=False, indent=4)
    # data = {"fit_Pop": fit_Pop.tolist(),
    #         "pro_Rate": pro_Rate.tolist(),
    #         "sensitivity": sensitivity.tolist(),
    #         "auto_R": auto_R.tolist(),
    #         "coopPayoff_Pop": coopPayoff_Pop.tolist(),
    #         "sigCost_Pop": sigCost_Pop.tolist(),
    #         "coopCost_Pop": coopCost_Pop.tolist(),
    #         "auto_pro_Rate": auto_pro_Rate.tolist()}
    # t = time.asctime().replace(":", "-" )
    # if Auto:
    #     file = f"New python\\auto json\\50 genotypes\gen {g+1} {lam} {tgen} {init_SN} {sig_Cost} {max_G} {clonal}.json"
    # else:
    #     file = f"New python\json\\50 genotypes\gen {g+1} {lam} {tgen} {init_SN} {sig_Cost} {max_G} {clonal}.json"
            
    # with open(file , "w") as f:
    #     json.dump(data, f,  ensure_ascii=False, indent=4)
    timestr = time.strftime("%d-%m %H-%M-%S")
    if Auto:
        file = f"New python\\auto json\\{timestr} {lam} {init_SN} {sig_Cost} {max_G} {clonal}.json"
    else:
        file = f"New python\\json\\{timestr} {lam} {init_SN} {sig_Cost} {max_G} {clonal}.json"
    
    data = {"fit_Evo": fit_Evo.tolist(),
            "pro_Rate_Evo": pro_Rate_Evo.tolist(),
            "sensitivity_Evo": sensitivity_Evo.tolist(),
            "auto_R_Evo": auto_R_Evo.tolist(),
            "coopPayoff_Evo": coopPayoff_Evo.tolist(),
            "sigCost_Evo": sigCost_Evo.tolist(),
            "coopCost_Evo": coopCost_Evo.tolist(),
            "auto_pro_Rate_Evo": auto_pro_Rate_Evo.tolist()}
    with open(file, "w") as f:
        json.dump(data, f,  ensure_ascii=False, indent=4)
    
    
    return file


def vary_signal(sig_Cost, auto, clonal):

    # t = time.time_ns()
    coop_Cost ,lam , K, mu_Cheats, max_G = 0.5, 0, 50.0, 10.0 ** -4, 5000
    returnValues = main_QS(coop_Cost, sig_Cost ,lam , K, mu_Cheats, max_G=max_G, Auto=auto, clonal=clonal)
    fit_Evo, pro_Rate_Evo, sensitivity_Evo, auto_R_Evo, coopPayoff_Evo, sigCost_Evo, coopCost_Evo, auto_pro_Rate_Evo = returnValues
    # print((time.time_ns()-t)* 10 **-9)
    print(f'{sig_Cost // 10 ** 8} done') 
    data = {"fit_Evo": fit_Evo.tolist(),
            "pro_Rate_Evo": pro_Rate_Evo.tolist(),
            "sensitivity_Evo": sensitivity_Evo.tolist(),
            "auto_R_Evo": auto_R_Evo.tolist(),
            "coopPayoff_Evo": coopPayoff_Evo.tolist(),
            "sigCost_Evo": sigCost_Evo.tolist(),
            "coopCost_Evo": coopCost_Evo.tolist(),
            "auto_pro_Rate_Evo": auto_pro_Rate_Evo.tolist()}
    t = time.asctime().replace(":", "-" )
    with open(f"Wang python\json\\autoreg\\{t} {coop_Cost} {sig_Cost // 10**8} {lam} {mu_Cheats} {auto} {max_G}.json", "w") as f:
        json.dump(data, f,  ensure_ascii=False, indent=4)

def vary_genotype(lam, Auto):
    clonal = True
    # t = time.time_ns()
    init_SN, sig_Cost, K, mu_Cheats, max_G = .2, 10**9, 50.0, 10.0 ** -4, 5000
    returnValues = main_QS(init_SN, sig_Cost ,lam , K, mu_Cheats, max_G=max_G, Auto=Auto, clonal=clonal)
    # fit_Evo, pro_Rate_Evo, sensitivity_Evo, auto_R_Evo, coopPayoff_Evo, sigCost_Evo, coopCost_Evo, auto_pro_Rate_Evo = returnValues
    # print((time.time_ns()-t)* 10 **-9)
    # print(f'{sig_Cost // 10 ** 8} done') 

if __name__ == "__main__":
    # joblib.Parallel(n_jobs=6)(joblib.delayed(vary_signal)(sig_Cost * 10**8) for sig_Cost in range(5,105,5))
    joblib.Parallel(n_jobs=6)(joblib.delayed(vary_genotype)(np.round(lam, decimals=1), True) for lam in np.arange(.5,50,step=1))
    joblib.Parallel(n_jobs=6)(joblib.delayed(vary_genotype)(np.round(lam, decimals=1), False) for lam in np.arange(.5,50,step=1))
    # Auto = True
    # clonal = False
    # # t = time.time_ns()
    # init_SN, sig_Cost, lam, K, mu_Cheats, max_G = .2, 10**9, 4, 50.0, 10.0 ** -4, 5000
    # file = main_QS(init_SN, sig_Cost ,lam , K, mu_Cheats, max_G=max_G, Auto=Auto, clonal=clonal)
    # graphNoThresholds.graph(file)