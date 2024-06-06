import numpy as np
import time
from helpers import *
import matplotlib.pyplot as plt

def main_QS(coop_Cost, sig_Cost ,lam , K, mu_Cheats, Auto=False):
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
    # maxeimum generation
    max_G = 1000
    # population size
    size_Pop = 5000
    # environments
    grid_Size = 100
    # baseline fitness
    baseline = 100.0
    # payoff for gene turned 'ON' & cooperation
    coop_Benefit = 1.5

    # mutation rate
    mu_Production = 0.01
    mu_Th_Signal = 0.01
    mu_R = 0.01

    # maximum cellular density (cells per microliter)
    max_CellDen = 10.0 ** 5
    # minimum cellular density
    min_CellDen = 10.0 ** 1.5
    # median cellular density
    median_CellDen = (max_CellDen + min_CellDen) * .5
    # baseline volume
    base_Volume = 10.0
    # signal decay rate
    decay_Rate = 10.0 ** -4.

    # initial testing environments
    env_CellDen = np.array(list(np.linspace(min_CellDen, max_CellDen,num=grid_Size)))
    # maximum cellular production rate
    max_ProRate = 2e-08
    # minimum cellular production rate
    min_ProRate = 0.
    # initial cellular production rate
    init_pro_Rate = 0.5e-08 
    # SD for mutation of cellular production rate
    mu_SD_ProRate = 0.01e-08

    # maximum signal concentration
    max_sigTh = 2e1
    # minimum signal concentration
    min_sigTh = 0.0001e1
    # initial signal threshold
    init_sig_Th = 0.5e1
    # SD for mutation of signal concentration
    mu_SD_sigTh = 0.01e1

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
    sig_Th = np.full((size_Pop,), init_sig_Th)
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
    coopEff_Pop = np.zeros(size_Pop)
    # cheats index
    index_Cheats = np.zeros(size_Pop)
    # selection index
    index_Select = np.zeros(size_Pop)

    # initialize result saving space
    # all lineages
    fit_Evo =  np.zeros(max_G)
    pro_Rate_Evo =  np.zeros(max_G)
    sig_Th_Evo =  np.zeros(max_G)
    auto_R_Evo =  np.zeros(max_G)
    coopPayoff_Evo =  np.zeros(max_G)
    sigCost_Evo =  np.zeros(max_G)
    coopCost_Evo =  np.zeros(max_G)
    coopEff_Evo =  np.zeros(max_G)
    auto_pro_Rate_Evo =  np.zeros(max_G)
    # non-cheats
    fit_nonCheats_Evo =  np.zeros(max_G)
    pro_Rate_nonCheats_Evo =  np.zeros(max_G)
    sig_Th_nonCheats_Evo =  np.zeros(max_G)
    auto_R_nonCheats_Evo =  np.zeros(max_G)
    coopPayoff_nonCheats_Evo =  np.zeros(max_G)
    sigCost_nonCheats_Evo =  np.zeros(max_G)
    coopCost_nonCheats_Evo =  np.zeros(max_G)
    coopEff_nonCheats_Evo =  np.zeros(max_G)
    auto_pro_Rate_nonCheats_Evo =  np.zeros(max_G)

    # record cheats
    numCheats_Evo =  np.zeros(max_G)
    ################################################################################
    # Evolution
    ################################################################################
    # initial evaluation
    if Auto:   
        coopPayoff_Pop, coopEff_Pop, coopCost_Pop, auto_pro_Rate, sigCost_Pop, fit_Pop = eval_genotype_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,coopEff_Pop,sigCost_Pop,auto_pro_Rate,pro_Rate,sig_Th,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen,K)
    else:
        coopPayoff_Pop, coopEff_Pop, coopCost_Pop, sigCost_Pop, fit_Pop = eval_genotype_No_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,coopEff_Pop,sigCost_Pop,pro_Rate,sig_Th,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen)
    
    for g in range(1,max_G):
        t = time.time_ns()
        # roulette-wheel selection
        # for n in range(0,size_Pop):
        #     index_Select[n] = int(fortune_wheel(fit_Pop))
        fit_pop_probability = fit_Pop / sum(fit_Pop)
        index_Select = np.random.choice(size_Pop, size_Pop, p=fit_pop_probability)
        
        # update
        #0.01... sec
        for n in range(len(index_Select)):
            pro_Rate[n] = pro_Rate[int(index_Select[n])]
            sig_Th[n] = sig_Th[int(index_Select[n])]
            fit_Pop[n] = fit_Pop[int(index_Select[n])]
            coopPayoff_Pop[n] = coopPayoff_Pop[int(index_Select[n])]
            coopCost_Pop[n] = coopCost_Pop[int(index_Select[n])]
            coopEff_Pop[n] = coopEff_Pop[int(index_Select[n])]
            sigCost_Pop[n] = sigCost_Pop[int(index_Select[n])]
            index_Cheats[n] = index_Cheats[int(index_Select[n])]
            if Auto: 
                auto_R[n] = auto_R[int(index_Select[n])]
                auto_pro_Rate[n] = auto_pro_Rate[int(index_Select[n])]
            else:
                pass
        
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
        sig_Th = mut_parameter(sig_Th,mu_Th_Signal,mu_SD_sigTh,min_sigTh,max_sigTh,index_Cheats,size_Pop)

        # mutate ratio of autoinduction production
        if Auto:
            auto_R = mut_parameter(auto_R,mu_R,mu_SD_R,min_R,max_R,index_Cheats,size_Pop)
        else: 
            pass
        # record cheats
        numCheats_Evo[g] = np.sum(index_Cheats)
        
        # save results  0 sec ish
        fit_Evo[g] = np.mean(fit_Pop)
        pro_Rate_Evo[g] = np.mean(pro_Rate)
        sig_Th_Evo[g] = np.mean(sig_Th)
        auto_R_Evo[g] = np.mean(auto_R)
        coopPayoff_Evo[g] = np.mean(coopPayoff_Pop)
        sigCost_Evo[g] = np.mean(sigCost_Pop)
        coopCost_Evo[g] = np.mean(coopCost_Pop)
        auto_pro_Rate_Evo[g] = np.mean(auto_pro_Rate)
        coopEff_Evo[g] = np.mean(coopEff_Pop)
        
        if numCheats_Evo[g] < size_Pop:

            index_nonCheats = np.ones(size_Pop)-index_Cheats
            fit_nonCheats_Evo[g] = np.dot(fit_Pop,index_nonCheats) / size_Pop
            pro_Rate_nonCheats_Evo[g] = np.dot(pro_Rate,index_nonCheats) / size_Pop
            sig_Th_nonCheats_Evo[g] = np.dot(sig_Th,index_nonCheats) / size_Pop
            auto_R_nonCheats_Evo[g] = np.dot(auto_R,index_nonCheats) / size_Pop
            coopPayoff_nonCheats_Evo[g] = np.dot(coopPayoff_Pop,index_nonCheats) / size_Pop
            sigCost_nonCheats_Evo[g] = np.dot(sigCost_Pop,index_nonCheats) / size_Pop
            coopCost_nonCheats_Evo[g] = np.dot(coopCost_Pop,index_nonCheats) / size_Pop
            auto_pro_Rate_nonCheats_Evo[g] = np.dot(auto_pro_Rate,index_nonCheats) / size_Pop
            coopEff_nonCheats_Evo[g] = np.dot(coopEff_Pop,index_nonCheats) / size_Pop
        else:
            fit_nonCheats_Evo[g] = fit_Evo[g]
            pro_Rate_nonCheats_Evo[g] = pro_Rate_Evo[g]
            sig_Th_nonCheats_Evo[g] = sig_Th_Evo[g]
            auto_R_nonCheats_Evo[g] = auto_R_Evo[g]
            coopPayoff_nonCheats_Evo[g] = coopPayoff_Evo[g]
            sigCost_nonCheats_Evo[g] = sigCost_Evo[g]
            coopCost_nonCheats_Evo[g] = coopCost_Evo[g]
            auto_pro_Rate_nonCheats_Evo[g] = auto_pro_Rate_Evo[g]
            coopEff_nonCheats_Evo[g] = coopEff_Evo[g]
        
        # genotype evaluation
        
        if Auto:   
            coopPayoff_Pop, coopEff_Pop, coopCost_Pop, auto_pro_Rate, sigCost_Pop, fit_Pop = eval_genotype_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,coopEff_Pop,sigCost_Pop,auto_pro_Rate,pro_Rate,sig_Th,auto_R,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen,K)
        else:
            coopPayoff_Pop, coopEff_Pop, coopCost_Pop, sigCost_Pop, fit_Pop = eval_genotype_No_Auto(fit_Pop,coopPayoff_Pop,coopCost_Pop,coopEff_Pop,sigCost_Pop,pro_Rate,sig_Th,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lam,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen)
        for i in range(len(fit_Pop)):
            if fit_Pop[i] < 0:
                fit_Pop[i] = 0
        
    return pro_Rate_Evo, sig_Th_Evo




pro_Rate_Evo, sig_Th_Evo = main_QS(0.5, 10.0 ** 9, 2, 50.0, 10.0 ** -4)
plt.plot(range(1000), pro_Rate_Evo)
plt.plot(range(1000), sig_Th_Evo)
plt.show()