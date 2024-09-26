import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import font_manager
import matplotlib
import json
import os
import numpy as np
import testing


# 
def graph(filename):
    
    with open(filename, "r") as f:
        data = json.load(f)
    fig, ax = plt.subplots(4,2)
    max_G= len(data["fit_Evo"])
    
    # maximum cellular density (cells per microliter)
    max_CellDen = 10.0 ** 5
    # minimum cellular density
    min_CellDen = 10.0 ** 1.5
    grid_Size = 100
    env_CellDen = np.array(list(np.linspace(min_CellDen, max_CellDen,num=grid_Size)))
    cm = plt.get_cmap("plasma")
    cNorm = colors.Normalize(vmin=min_CellDen, vmax= max_CellDen)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    for j in env_CellDen:
        ax[0,0].plot(range(max_G), np.array(data["sensitivity_Evo"])*np.array(data["pro_Rate_Evo"]) /  10.0 ** -4 * j ** 2, color=scalarMap.to_rgba(j)) # (10.0 ** 1.5)**2
    ax[0,0].set_title("testing")
    ax[0,0].plot(range(max_G), np.full((max_G,), 50016 ), color="black")

    
    ax[3,1].plot(range(max_G), data["fit_Evo"])
    ax[3,1].set_title("fit_Evo")

    ax[1,0].plot(range(max_G), data["pro_Rate_Evo"])
    ax[1,0].set_title("pro_Rate_Evo")

    ax[2,0].plot(range(max_G), data["sensitivity_Evo"])
    ax[2,0].set_title("sensitivity_Evo")

    ax[0,1].plot(range(max_G), data["coopPayoff_Evo"])
    ax[0,1].set_title("coopPayoff_Evo")

    ax[1,1].plot(range(max_G), data["sigCost_Evo"])
    ax[1,1].set_title("sigCost_Evo")

    ax[2,1].plot(range(max_G), data["coopCost_Evo"])
    ax[2,1].set_title("coopCost_Evo")

    # ax[3,1].plot(range(max_G), np.array(data["auto_pro_Rate_Evo"])+np.array(data["pro_Rate_Evo"]))
    # ax[3,1].set_title("total production")

    ax[3,0].plot(range(max_G), np.array(data["sensitivity_Evo"])*np.array(data["pro_Rate_Evo"])/  10.0 ** -4 )
    ax[3,0].set_title("S*p/u")
    ax[3,0].plot(range(max_G), np.full((max_G,), 1.509*10**-5 ), color="#3fe61e", linewidth=2.0, linestyle="dashdot")
    ax[3,0].plot(range(max_G), np.full((max_G,), .498*10**-5 ), color="black", linewidth=2.0, linestyle="dashdot")
    
    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.09, hspace=.524)
    plt.show()


def graph_auto(filename):
    
    with open(filename, "r") as f:
        data = json.load(f)
    fig, ax = plt.subplots(4,2)
    max_G= len(data["fit_Evo"])
    
    # maximum cellular density (cells per microliter)
    max_CellDen = 10.0 ** 5
    # minimum cellular density
    min_CellDen = 10.0 ** 1.5
    grid_Size = 100
    env_CellDen = np.array(list(np.linspace(min_CellDen, max_CellDen,num=grid_Size)))
    cm = plt.get_cmap("plasma")
    cNorm = colors.Normalize(vmin=min_CellDen, vmax= max_CellDen)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    for j in env_CellDen:
        ax[0,0].plot(range(max_G), np.array(data["sensitivity_Evo"])*np.array(data["pro_Rate_Evo"]) /  10.0 ** -4 * j ** 2, color=scalarMap.to_rgba(j)) # (10.0 ** 1.5)**2
    ax[0,0].set_title("testing")
    ax[0,0].plot(range(max_G), np.full((max_G,), 50016 ), color="black")

    
    ax[3,1].plot(range(max_G), data["fit_Evo"])
    ax[3,1].set_title("fit_Evo")

    ax[1,0].plot(range(max_G), data["pro_Rate_Evo"], color="blue")
    ax[1,0].plot(range(max_G), np.array(data["auto_pro_Rate_Evo"])+np.array(data["pro_Rate_Evo"]) , color="red")
    ax[1,0].set_title("Basil production(B) and total(R)")

    ax[2,0].plot(range(max_G), data["sensitivity_Evo"])
    ax[2,0].set_title("sensitivity_Evo")

    ax[0,1].plot(range(max_G), data["coopPayoff_Evo"], color="blue")
    ax[0,1].plot(range(max_G), data["coopCost_Evo"], color="red")
    ax[0,1].set_title("coopPayoff_Evo(B)/coopCost_Evo(R)")

    ax[1,1].plot(range(max_G), data["sigCost_Evo"])
    ax[1,1].set_title("sigCost_Evo")

    ax[2,1].plot(range(max_G), data["auto_R_Evo"])
    ax[2,1].set_title("auto_R_Evo")

    # ax[3,1].plot(range(max_G), np.array(data["auto_pro_Rate_Evo"])+np.array(data["pro_Rate_Evo"]))
    # ax[3,1].set_title("total production")
    ax[3,0].plot(range(max_G), np.full((max_G,), 1.509*10**-5 ), color="#3fe61e", linewidth=2.0, linestyle="dashdot")
    ax[3,0].plot(range(max_G), np.full((max_G,), .498*10**-5 ), color="black", linewidth=2.0, linestyle="dashdot")
    
    ax[3,0].plot(range(max_G), np.array(data["sensitivity_Evo"])*np.array(data["pro_Rate_Evo"])/  10.0 ** -4 )
    ax[3,0].set_title("S*p/u")

    
    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.09, hspace=.524)
    plt.show()


def graph_multiple(split_value, folder_name, minimum, maximum, step):
    cm = plt.get_cmap("plasma")
    cNorm = colors.Normalize(vmin=minimum, vmax= maximum)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig, ax = plt.subplots(4,2)
    max_G= 0
    prod = []
    threshold  = []
    costs = []
    
    for path, directories, files in os.walk(folder_name):
        for file in files:
            try:
                # cost = int(file.split(" ")[split_value])
                cost = np.round(float(file.split(" ")[split_value]), decimals=1)
                filename = path+"\\"+file
                with open(filename, "r") as f:
                    data = json.load(f)
                max_G = len(data["fit_Evo"])
                ax[0,0].plot(range(max_G), np.array(data["sensitivity_Evo"])*np.array(data["pro_Rate_Evo"]) /  10.0 ** -4 , color=scalarMap.to_rgba(cost))
                ax[0,0].set_title("testing")

                ax[3,1].plot(range(max_G), data["fit_Evo"], color=scalarMap.to_rgba(cost))
                ax[3,1].set_title("fit_Evo vs generation")

                ax[2,0].plot(range(max_G), data["pro_Rate_Evo"], color=scalarMap.to_rgba(cost))
                ax[2,0].set_title("pro_Rate_Evo vs generation")

                ax[3,0].plot(range(max_G), data["sensitivity_Evo"], color=scalarMap.to_rgba(cost))
                ax[3,0].set_title("sensitivity_Evo vs generation")

                ax[0,1].plot(range(max_G), data["coopPayoff_Evo"], color=scalarMap.to_rgba(cost))
                ax[0,1].set_title("coopPayoff_Evo vs generation")
                
                ax[1,1].plot(range(max_G), data["sigCost_Evo"], color=scalarMap.to_rgba(cost))
                ax[1,1].set_title("sigCost_Evo vs generation")

                ax[2,1].plot(range(max_G), data["coopCost_Evo"], color=scalarMap.to_rgba(cost))
                ax[2,1].set_title("coopCost_Evo vs generation")

                prod.append(np.mean(data["pro_Rate_Evo"][-50:]))
                threshold.append(np.mean(data["sensitivity_Evo"][-50:]))
                costs.append(cost)
            except:
                pass
    # color="#3fe61e"
    ax[0,0].plot(range(max_G), np.full((max_G,), 1.509*10**-5 ), color="#3fe61e", linewidth=2.0, linestyle="dashdot")
    ax[0,0].plot(range(max_G), np.full((max_G,), .498*10**-5 ), color="black", linewidth=2.0, linestyle="dashdot")
    
    # ax[0,0].scatter(costs, prod,color="blue")

    # another = ax[0,0].twinx()
    # another.scatter(costs, threshold,color="red")
    # ax[0,0].set_xlim([0,50])
    # another.spines['right'].set_color('red')
    # another.spines['left'].set_color('blue')
    # ax[0,0].tick_params(axis="y", colors="blue")
    # another.tick_params(axis="y", colors="red")
    # ax[0,0].set_xlabel("Average genotype per group")
    # ax[0,0].set_ylabel("Evolved production rate")
    # another.set_ylabel("Evolved signal Threshold")
    # ax[0,0].set_title("Evolved production rate (red) and Signal Threshold(red) vs genotype")

    ax[1,0].scatter(prod, threshold, color= scalarMap.to_rgba(costs))
    ax[1,0].set_xlabel("evolved production rate")
    ax[1,0].set_ylabel("evolved siginaling threshold")

    plt.colorbar(scalarMap, ax=ax[1,0])
    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.09, hspace=.524)
    plt.show()


def graph_last_gen(file):
    
    with open(file, "r") as f:
        data = json.load(f)
    fig, ax = plt.subplots(4,2)
    binNo = 100
    fit_Pop = data["fit_Pop"]
    cm = plt.get_cmap("plasma")
    q1, mead, q3 = np.quantile(fit_Pop, [.25,.5,.75])
    cNorm = colors.CenteredNorm(vcenter=mead, halfrange=(q3-q1))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

    name = "fit_Pop"
    bins = np.digitize(data[name], bins=np.histogram_bin_edges(data[name],bins=binNo))
    bincolor= []
    x = np.histogram_bin_edges(data[name],bins=binNo)[:-1]
    counts =[]
    for i in range(binNo):
        counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
        bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    ax[0,0].bar(x, counts, color= bincolor, width=(np.max(data[name])-np.min(data[name]))/110)
    ax[0,0].set_title("fit_Pop historgram")

    name = "coopPayoff_Pop"
    bins = np.digitize(data[name], bins=np.histogram_bin_edges(data[name],bins=binNo))
    bincolor= []
    x = np.histogram_bin_edges(data[name],bins=binNo)[:-1]
    counts =[]
    for i in range(binNo):
        counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
        bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    ax[1,0].bar(x, counts, color= bincolor, width=(np.max(data[name])-np.min(data[name]))/110)
    ax[1,0].set_title("coopPayoff_Pop historgram")

    
    name = "sigCost_Pop"
    bins = np.digitize(data[name], bins=np.histogram_bin_edges(data[name],bins=binNo))
    bincolor= []
    x = np.histogram_bin_edges(data[name],bins=binNo)[:-1]
    counts =[]
    for i in range(binNo):
        counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
        bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    ax[2,0].bar(x, counts, color=bincolor, width=(np.max(data[name])-np.min(data[name]))/110)
    ax[2,0].set_title("sigCost_Pop historgram")
 
    name = "coopCost_Pop"
    bins = np.digitize(data[name], bins=np.histogram_bin_edges(data[name],bins=binNo))
    bincolor= []
    x = np.histogram_bin_edges(data[name],bins=binNo)[:-1]
    counts =[]
    for i in range(binNo):
        counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
        bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    ax[3,0].bar(x, counts, color= bincolor, width=(np.max(data[name])-np.min(data[name]))/110)
    ax[3,0].set_title("coopCost_Pop historgram")

    name = "pro_Rate"
    bins = np.digitize(data[name], bins=np.histogram_bin_edges(data[name],bins=binNo))
    bincolor= []
    x = np.histogram_bin_edges(data[name],bins=binNo)[:-1]
    counts =[]
    for i in range(binNo):
        counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
        bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    ax[0,1].bar(x, counts, color= bincolor, width=(np.max(data[name])-np.min(data[name]))/110)
    ax[0,1].set_title("pro_Rate historgram")

    name = "sensitivity"
    bins = np.digitize(data[name], bins=np.histogram_bin_edges(data[name],bins=binNo))
    bincolor= []
    x = np.histogram_bin_edges(data[name],bins=binNo)[:-1]
    counts =[]
    for i in range(binNo):
        counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
        bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    ax[1,1].bar(x, counts, color= bincolor, width=(np.max(data[name])-np.min(data[name]))/110)
    ax[1,1].set_title("sensitivity historgram")


    ax[2,1].scatter(data["pro_Rate"], data["sensitivity"],color=scalarMap.to_rgba(fit_Pop), s=3)
    ax[2,1].set_title("production rate vs sensitivity")

    data["ratio"] = np.array(data["sensitivity"])* np.array(data["pro_Rate"]) / 10 ** -4
    name = "ratio"
    bins = np.digitize(data[name], bins=np.histogram_bin_edges(data[name],bins=binNo))
    bincolor= []
    x = np.histogram_bin_edges(data[name],bins=binNo)[:-1]
    counts =[]
    for i in range(binNo):
        counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
        bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    ax[3,1].bar(x, counts, color= bincolor, width=(np.max(data[name])-np.min(data[name]))/110)
    ax[3,1].set_title("ratio historgram")

    

    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.09, hspace=.524)
    plt.show()


def graph_multiple_auto(split_value, folder_name, minimum, maximum, step):
    cm = plt.get_cmap("plasma")
    cNorm = colors.Normalize(vmin=minimum, vmax= maximum)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig, ax = plt.subplots(4,2)
    max_G= 0
    prod = []
    threshold  = []
    costs = []
    
    for path, directories, files in os.walk(folder_name):
        for file in files:
            try:
                # cost = int(file.split(" ")[split_value])
                cost = np.round(float(file.split(" ")[split_value]), decimals=1)
                filename = path+"\\"+file
                with open(filename, "r") as f:
                    data = json.load(f)
                max_G = len(data["fit_Evo"])
                ax[0,0].plot(range(max_G), np.array(data["sensitivity_Evo"])*(np.array(data["pro_Rate_Evo"])+ np.array(data["auto_pro_Rate_Evo"]))/  10.0 ** -4 , color=scalarMap.to_rgba(cost))
                ax[0,0].set_title("testing")

                ax[3,1].plot(range(max_G), data["fit_Evo"], color=scalarMap.to_rgba(cost))
                ax[3,1].set_title("fit_Evo vs generation")

                ax[2,0].plot(range(max_G), data["pro_Rate_Evo"], color=scalarMap.to_rgba(cost))
                ax[2,0].set_title("pro_Rate_Evo vs generation")

                ax[3,0].plot(range(max_G), data["sensitivity_Evo"], color=scalarMap.to_rgba(cost))
                ax[3,0].set_title("sensitivity_Evo vs generation")

                ax[0,1].plot(range(max_G), data["coopPayoff_Evo"], color=scalarMap.to_rgba(cost))
                ax[0,1].set_title("coopPayoff_Evo vs generation")
                
                ax[1,1].plot(range(max_G), data["sigCost_Evo"], color=scalarMap.to_rgba(cost))
                ax[1,1].set_title("sigCost_Evo vs generation")

                ax[2,1].plot(range(max_G), data["coopCost_Evo"], color=scalarMap.to_rgba(cost))
                ax[2,1].set_title("coopCost_Evo vs generation")

                prod.append(np.mean(data["pro_Rate_Evo"][-50:]))
                threshold.append(np.mean(data["sensitivity_Evo"][-50:]))
                costs.append(cost)
            except:
                pass
    # color="#3fe61e"
    ax[0,0].plot(range(max_G), np.full((max_G,), 1.509*10**-5 ), color="#3fe61e", linewidth=2.0, linestyle="dashdot")
    ax[0,0].plot(range(max_G), np.full((max_G,), .498*10**-5 ), color="black", linewidth=2.0, linestyle="dashdot")
    
    # ax[0,0].scatter(costs, prod,color="blue")

    # another = ax[0,0].twinx()
    # another.scatter(costs, threshold,color="red")
    # ax[0,0].set_xlim([0,50])
    # another.spines['right'].set_color('red')
    # another.spines['left'].set_color('blue')
    # ax[0,0].tick_params(axis="y", colors="blue")
    # another.tick_params(axis="y", colors="red")
    # ax[0,0].set_xlabel("Average genotype per group")
    # ax[0,0].set_ylabel("Evolved production rate")
    # another.set_ylabel("Evolved signal Threshold")
    # ax[0,0].set_title("Evolved production rate (red) and Signal Threshold(red) vs genotype")

    ax[1,0].scatter(prod, threshold, color= scalarMap.to_rgba(costs))
    ax[1,0].set_xlabel("evolved production rate")
    ax[1,0].set_ylabel("evolved siginaling threshold")

    plt.colorbar(scalarMap, ax=ax[1,0])
    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.09, hspace=.524)
    plt.show()


def optimum(p,r,s):
    a =  1.00450376e+03 * np.exp(-5.04190644e-02 / s) - 1.00444712e+03
    b =  2.0118419e-09 * s + 9.9727154e-08
    return 100 - 10** 8*np.abs(a *p +b - (p*(1+r)))



def graph_last_gen_standard(file, lam, gen):
    # print(file.split("\\")[1].split(" ")[0])
    if file.split("\\")[1].split(" ")[0] == "auto":
        auto = True
    else:
        auto = False
    generation = gen
    with open(file, "r") as f:
        data = json.load(f)
    fig, ax = plt.subplots(1,2, figsize=(18, 6))
    binNo = 50
    fit_Pop = data["fit_Pop"]
    cm = plt.get_cmap("winter")
    q1, mead, q3 = np.quantile(fit_Pop, [.25,.5,.75])
    cNorm = colors.CenteredNorm(vcenter=mead, halfrange=(q3-q1))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    name = "fit_Pop"
    bins = np.digitize(data[name], bins=np.histogram_bin_edges(data[name],bins=binNo))
    bincolor= []
    x = np.histogram_bin_edges(data[name],bins=binNo)[:-1]
    fs = 35
    counts =[]
    for i in range(binNo):
        counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
        bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    ax[0].bar(x, counts, color= bincolor, width=(np.max(data[name])-np.min(data[name]))/55)
    ax[0].set_title("Fitness Histogram", fontsize=fs)
    ax[0].set_xlabel("Fitness (FU)", fontsize=fs)

    # difference =  np.array(data["coopPayoff_Pop"]) /  np.array(data["sigCost_Pop"])
    # bins = np.digitize(difference, bins=np.histogram_bin_edges(difference,bins=binNo))
    # bincolor= []
    # x = np.histogram_bin_edges(difference,bins=binNo)[:-1]
    # counts =[]
    # for i in range(binNo):
    #     counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
    #     bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    # if (np.max(difference)-np.min(difference)) ==0:
    #     width = .05
    # else: 
    #     width=(np.max(difference)-np.min(difference))/50
    # ax[1].bar(x, counts, color= bincolor, width=width)
    # ax[1].set_title("Cooperation Payoff vs. Cost Difference")
    # ax[1].set_xlabel("Cooperation Payoff vs. Cost Difference (FU)")


    # ax[1].scatter(np.array(data["coopPayoff_Pop"]),  np.array(data["coopCost_Pop"])+np.array(data["sigCost_Pop"]), color=scalarMap.to_rgba(fit_Pop), s=10)
    # ax[1].set_title("Cooperation Payoff vs. Cost Difference", fontsize=fs)
    # ax[1].set_xlabel("Cooperation Payoff (FU)", fontsize=fs)
    # ax[1].set_ylabel("Cooperation Cost (FU)", fontsize=fs)

    if auto:
        production = np.array(data["pro_Rate"]) + np.array(data["auto_pro_Rate"])
    else:
        production = np.array(data["pro_Rate"]) 
    extreme = np.array([np.min(production), np.max(production)])
    # ax[2].plot(extreme, extreme * 50016 / 10 ** -4, color="black", linewidth = 2,  linestyle="dashed", dashes=(5, 10))


    sortestlist = sorted(zip(data["fit_Pop"], production, data["sensitivity"]))
    production = [p for f, p, s in sortestlist]
    sensitivity =[s for f, p, s in sortestlist]
    fitness = [f for f, p, s in sortestlist]
    ax[1].scatter(production, sensitivity,color=scalarMap.to_rgba(fitness), s=2)
    ax[1].set_title("Production Rate vs. Sensitivty", fontsize=fs)
    ax[1].set_xlabel("Production Rate $\left(\\frac{\mu M}{s}\\right)$", fontsize=fs)
    ax[1].set_ylabel("Sensitivity $\left(\mu M\\right)$", fontsize=fs)

    x = np.linspace(0, 25*10**-9, num = 100)
    # x = np.linspace(np.min(production),np.max(production), num = 100)
    # ax[2].plot(x, , color="black", linewidth = 2,  linestyle="dashed", )
    ax[1].plot(x, 10**-4 * 1.62*10**-5 / x , color="black", linewidth=2.0, linestyle="dashdot")
    ax[1].plot(x, 10**-4 * .4992*10**-5 / x, color="black", linewidth=2.0, linestyle="dashed")
    if auto:
        ax[1].set_xlim([0, 23 * 10 **-9])
        ax[1].set_ylim([0, .3])
    else:
        ax[1].set_xlim([0, 9 * 10 **-9])
        ax[1].set_ylim([0, .75])
    fig.suptitle(f"Generation {generation}", fontsize = fs)
    # 
   

    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.162, hspace=.524)
    # plt.show()
    if auto:
        plt.savefig(f"Pictures\\auto {lam} {generation}.png")
    else:
        plt.savefig(f"Pictures\\no auto {lam} {generation}.png")


def graph_multiple_standard1(split_value, folder_name, minimum, maximum, step):

    fig, ax = plt.subplots(1,3, figsize=(25, 6))
    cm = plt.get_cmap("plasma")
    cNorm = colors.Normalize(vmin=minimum, vmax= maximum)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    max_G= 5000
    prod = []
    totalprod =[]
    threshold  = []
    splitname = ""
    splitname = "Number of Genotypes per Testing Environment"
    # fig.suptitle("Genetic Mixing", fontsize= 20)
    if split_value == 6:
        splitname = "Cost of Signaling $(10^8 FU)$"
        fig.suptitle("Variation of Cost of Singaling", fontsize= 20)
    if split_value == 7:
        splitname = "Number of Genotypes per Testing Environment"
        fig.suptitle("Genetic Mixing", fontsize= 20)
    costs = []
    fs = 30
    for path, directories, files in os.walk(folder_name):
        for file in files:
            # cost = int(file.split(" ")[split_value])
            cost = np.round(float(file.split(" ")[split_value]), decimals=1) 
            filename = path+"\\"+file
            with open(filename, "r") as f:
                data = json.load(f)
            max_G = len(data["fit_Evo"])
           
            p = np.array(data["pro_Rate_Evo"])
            r = np.array(data["auto_pro_Rate_Evo"])
            s = np.array(data["sensitivity_Evo"])

            ax[0].plot(range(max_G), np.array(data["sensitivity_Evo"])*(p) /  10.0 ** -4, color=scalarMap.to_rgba(cost))

            # ax[0].plot(range(max_G), s * 10 **-4 / (p+r), color=scalarMap.to_rgba(cost))
            ax[1].plot(range(max_G), np.array(data["pro_Rate_Evo"]), color=scalarMap.to_rgba(cost))
            ax[2].plot(range(max_G), data["sensitivity_Evo"], color=scalarMap.to_rgba(cost))
    plt.colorbar(scalarMap, ax=ax[0], label= splitname)
    ax[0].plot(range(max_G), np.full((max_G,), 1.62*10**-5 ), color="black", linewidth=2.0, linestyle="dashdot")
    ax[0].plot(range(max_G), np.full((max_G,), .4992*10**-5 ), color="black", linewidth=2.0, linestyle="dashed")
    ax[0].set_title("$\\frac{Sn_i p_{total}}{u}$ Evolution", fontsize=fs)
    ax[0].set_xlabel("Generation", fontsize=fs)
    ax[0].set_ylabel('$\\frac{Sn_i p_{total}}{u}$  $\left(\mu L\\right)$', fontsize=fs)

    ax[1].set_title("Total Production Rate Evolution", fontsize=fs)
    ax[1].set_xlabel("Generation", fontsize=fs)
    ax[1].set_ylabel('Production Rate $\left(\\frac{\mu M}{s}\\right)$', fontsize=fs)

    ax[2].set_title("Sensitivity Evolution", fontsize=fs)
    ax[2].set_xlabel("Generation", fontsize=fs)
    ax[2].set_ylabel('Sensitivity', fontsize=fs)

    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.162, hspace=.524)
    plt.savefig(f"Pictures\\No Auto genotype 1.png")

def graph_multiple_standard2(split_value, folder_name, minimum, maximum, step):
    fig, ax = plt.subplots(1,3, figsize=(25, 6))
    cm = plt.get_cmap("plasma")
    cNorm = colors.Normalize(vmin=minimum, vmax= maximum)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    max_G= 5000
    prod = []
    totalprod =[]
    threshold  = []
    splitname = ""
    splitname = "Number of Genotypes per Testing Environment"
    # fig.suptitle("Genetic Mixing", fontsize= 20)
    if split_value == 6:
        splitname = "Cost of Signaling $(10^8 FU)$"
        fig.suptitle("Variation of Cost of Singaling", fontsize= 20)
    if split_value == 7:
        splitname = "Number of Genotypes per Testing Environment"
        fig.suptitle("Genetic Mixing", fontsize= 20)
    costs = []
    fs = 30
    for path, directories, files in os.walk(folder_name):
        for file in files:
            # cost = int(file.split(" ")[split_value])
            cost = np.round(float(file.split(" ")[split_value]), decimals=1) 
            filename = path+"\\"+file
            with open(filename, "r") as f:
                data = json.load(f)
            max_G = len(data["fit_Evo"])
           
            p = np.array(data["pro_Rate_Evo"])
            r = np.array(data["auto_pro_Rate_Evo"])
            s = np.array(data["sensitivity_Evo"])

            ax[0].plot(range(max_G), np.array(data["sensitivity_Evo"])*(p+ r) /  10.0 ** -4, color=scalarMap.to_rgba(cost))

            # ax[0].plot(range(max_G), s * 10 **-4 / (p+r), color=scalarMap.to_rgba(cost))
            ax[1].plot(range(max_G), p+r, color=scalarMap.to_rgba(cost))
            ax[2].plot(range(max_G), data["sensitivity_Evo"], color=scalarMap.to_rgba(cost))
    plt.colorbar(scalarMap, ax=ax[0], label= splitname)
    ax[0].plot(range(max_G), np.full((max_G,), 1.62*10**-5 ), color="black", linewidth=2.0, linestyle="dashdot")
    ax[0].plot(range(max_G), np.full((max_G,), .4992*10**-5 ), color="black", linewidth=2.0, linestyle="dashed")
    ax[0].set_title("$\\frac{Sn_i p_{total}}{u}$ Evolution", fontsize=fs)
    ax[0].set_xlabel("Generation", fontsize=fs)
    ax[0].set_ylabel('$\\frac{Sn_i p_{total}}{u}$  $\left(\mu L\\right)$', fontsize=fs)

    ax[1].set_title("Total Production Rate Evolution", fontsize=fs)
    ax[1].set_xlabel("Generation", fontsize=fs)
    ax[1].set_ylabel('Production Rate $\left(\\frac{\mu M}{s}\\right)$', fontsize=fs)

    ax[2].set_title("Sensitivity Evolution", fontsize=fs)
    ax[2].set_xlabel("Generation", fontsize=fs)
    ax[2].set_ylabel('Sensitivity', fontsize=fs)

    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.162, hspace=.524)
    plt.savefig(f"Pictures\\Auto genotype 1.png")


if __name__ == "__main__":
    file = "New python\json\\25-09 14-27-54 3 0.2 1000000000 5000 True.json"
    plt.rcParams["savefig.directory"]  = "Pictures"
    font_path = 'Inter-Regular.otf'  # Your font path goes here
    font_manager.fontManager.addfont("C:\\Users\\matth\\Downloads\\Lato\\Lato-Regular.ttf")
    prop = font_manager.FontProperties(fname=font_path)
    font = {'family' : 'Lato',
        # 'weight' : 'bold',
        'size'   : 15}
    # print(matplotlib.font_manager.get_font_names())
    matplotlib.rc('font', **font)
    graph(file)
    # graph_last_gen(file)
    # graph_multiple(2, "New python\json\\actuall 50 genotype evo", 0, 10, .5)
    # fig, ax = plt.subplots(1, figsize=(8, 6))
    # steps = np.linspace(0, 6e-5, num=10000)
    # y = [testing.fitness_sum(c) for c in steps]
    # minmax = [np.min(y), np.max(y) + 10]


    # plt.plot([4.992e-6, 4.992e-6], minmax,color ="black", linestyle="dashed")
    # plt.plot([1.620e-5, 1.62e-5], minmax,color ="black", linestyle="dashdot")

    # plt.plot(steps, [testing.fitness_sum(c) for c in steps])
    # plt.xlabel("$Sn_i S_{N_j}$",fontsize=20)
    # plt.ylabel("Net Cooperation Impact",fontsize=20)
    
    
    # plt.show()
    # lam=""
    # gen=5000
    # graph_last_gen_standard(f"New python\\json\\clonal gen\\gen {gen}.json", lam, gen)
    # graph_multiple_standard2(2,"New python\\auto json\\10 geno", 0, 10, .1)
    # graph_multiple_standard1(2,"New python\json\\10 genotypes", 0, 10, .1)
    
    # for gen in [250,750,1000,5000]:
    #     graph_last_gen_standard(f"New python\\json\\clonal gen\\gen {gen}.json", "0", gen)
    #     graph_last_gen_standard(f"New python\\json\\generations for 10\\gen {gen} 6.0.json", "6.0", gen)
    
    # for gen in [250,1250,1500,5000]:
        
    #     graph_last_gen_standard(f"New python\\auto json\\generations for 10\\gen {gen} 3.0.json", "3.0", gen)