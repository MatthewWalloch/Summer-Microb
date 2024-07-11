import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import json
import os
import numpy as np


# 
def graph(filename): 
    with open(filename, "r") as f:
        data = json.load(f)
    fig, ax = plt.subplots(4,2)
    max_G= len(data["fit_Evo"])
    ax[0,0].plot(range(max_G), data["fit_Evo"])
    ax[0,0].set_title("fit_Evo")

    ax[1,0].plot(range(max_G), data["pro_Rate_Evo"])
    ax[1,0].set_title("pro_Rate_Evo")

    ax[2,0].plot(range(max_G), data["sig_Th_Evo"])
    ax[2,0].set_title("sig_Th_Evo")

    ax[0,1].plot(range(max_G), data["coopPayoff_Evo"])
    ax[0,1].set_title("coopPayoff_Evo")

    ax[1,1].plot(range(max_G), data["sigCost_Evo"])
    ax[1,1].set_title("sigCost_Evo")

    ax[2,1].plot(range(max_G), data["coopCost_Evo"])
    ax[2,1].set_title("coopCost_Evo")

    ax[3,1].plot(range(max_G), np.array(data["auto_pro_Rate_Evo"])+np.array(data["pro_Rate_Evo"]))
    ax[3,1].set_title("total production")

    ax[3,0].plot(range(max_G), data["auto_R_Evo"])
    ax[3,0].set_title("auto_R_Evo")
    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.09, hspace=.524)
    plt.show()

def graph_multiple(split_value, folder_name, minimum, maximum, step):
    cm = plt.get_cmap("plasma")
    cNorm = colors.Normalize(vmin=minimum, vmax= maximum)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig, ax = plt.subplots(4,2)

    max_G= 5000
    prod = []
    totalprod =[]
    threshold  = []
    costs = []
    ax[0,0].plot(range(max_G), np.full((max_G,), 50016 ), color="black")
    for path, directories, files in os.walk(folder_name):
        for file in files:
            # cost = int(file.split(" ")[split_value])
            cost = np.round(float(file.split(" ")[split_value]), decimals=1) // 10**8
            filename = path+"\\"+file
            with open(filename, "r") as f:
                data = json.load(f)
            max_G = len(data["fit_Evo"])
            ax[0,0].plot(range(max_G), np.array(data["sig_Th_Evo"]) /np.array(data["pro_Rate_Evo"]) *  10.0 ** -4, color=scalarMap.to_rgba(cost))
            ax[0,0].set_title("testing")

            ax[3,1].plot(range(max_G), data["fit_Evo"], color=scalarMap.to_rgba(cost))
            ax[3,1].set_title("fit_Evo vs generation")

            ax[2,0].plot(range(max_G), np.array(data["pro_Rate_Evo"])+np.array(data["auto_pro_Rate_Evo"]), color=scalarMap.to_rgba(cost))
            ax[2,0].set_title("total production vs generation")

            ax[3,0].plot(range(max_G), data["sig_Th_Evo"], color=scalarMap.to_rgba(cost))
            ax[3,0].set_title("sig_Th_Evo vs generation")

            ax[0,1].plot(range(max_G), np.array(data["coopPayoff_Evo"]), color=scalarMap.to_rgba(cost))
            ax[0,1].set_title("coopPayoff_Evo vs generation")
            
            ax[1,1].plot(range(max_G), data["sigCost_Evo"], color=scalarMap.to_rgba(cost))
            ax[1,1].set_title("sigCost_Evo vs generation")

            ax[2,1].plot(range(max_G), data["coopCost_Evo"], color=scalarMap.to_rgba(cost))
            ax[2,1].set_title("coopCost_Evo vs generation")

            prod.append(np.mean(data["pro_Rate_Evo"][-50:]))
            totalprod.append(np.mean(data["pro_Rate_Evo"][-50:])+ np.mean(data["auto_pro_Rate_Evo"][-50:]))
            threshold.append(np.mean(data["sig_Th_Evo"][-50:]))
            costs.append(cost)

    
    # ax[0,0].scatter(costs, totalprod,color="blue", marker="^")
    # ax[0,0].scatter(costs, prod, color="blue")
    # ax[0,0].set_ylim([0,1.5*10**-8])
    # another = ax[0,0].twinx()
    # another.scatter(costs, threshold,color="red")
    # another.set_ylim([0,20])
    # ax[0,0].set_xlim([0,10])
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


def graph_multiple(split_value, folder_name, minimum, maximum, step):
    cm = plt.get_cmap("plasma")
    cNorm = colors.Normalize(vmin=minimum, vmax= maximum)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

    max_G= 5000
    prod = []
    totalprod =[]
    threshold  = []
    costs = []
    
    for path, directories, files in os.walk(folder_name):
        for file in files:
            # cost = int(file.split(" ")[split_value])
            cost = np.round(float(file.split(" ")[split_value]), decimals=1) // 10**8
            filename = path+"\\"+file
            with open(filename, "r") as f:
                data = json.load(f)
            max_G = len(data["fit_Evo"])
            plt.plot(range(max_G), np.array(data["sig_Th_Evo"]) /np.array(data["pro_Rate_Evo"]) *  10.0 ** -4, color=scalarMap.to_rgba(cost))
            plt.title("testing")

            # ax.plot(range(max_G), data["fit_Evo"], color=scalarMap.to_rgba(cost))
            # ax.set_title("fit_Evo vs generation")

            # ax.plot(range(max_G), np.array(data["pro_Rate_Evo"])+np.array(data["auto_pro_Rate_Evo"]), color=scalarMap.to_rgba(cost))
            # ax.set_title("total production vs generation")

            # ax.plot(range(max_G), data["sig_Th_Evo"], color=scalarMap.to_rgba(cost))
            # ax.set_title("sig_Th_Evo vs generation")

            # ax.plot(range(max_G), np.array(data["coopPayoff_Evo"]), color=scalarMap.to_rgba(cost))
            # ax.set_title("coopPayoff_Evo vs generation")
            
            # ax.plot(range(max_G), data["sigCost_Evo"], color=scalarMap.to_rgba(cost))
            # ax.set_title("sigCost_Evo vs generation")

            # ax.plot(range(max_G), data["coopCost_Evo"], color=scalarMap.to_rgba(cost))
            # ax.set_title("coopCost_Evo vs generation")

            # prod.append(np.mean(data["pro_Rate_Evo"][-50:]))
            # totalprod.append(np.mean(data["pro_Rate_Evo"][-50:])+ np.mean(data["auto_pro_Rate_Evo"][-50:]))
            # threshold.append(np.mean(data["sig_Th_Evo"][-50:]))
            # costs.append(cost)

    
    # ax[0,0].scatter(costs, totalprod,color="blue", marker="^")
    # ax[0,0].scatter(costs, prod, color="blue")
    # ax[0,0].set_ylim([0,1.5*10**-8])
    # another = ax[0,0].twinx()
    # another.scatter(costs, threshold,color="red")
    # another.set_ylim([0,20])
    # ax[0,0].set_xlim([0,10])
    # another.spines['right'].set_color('red')
    # another.spines['left'].set_color('blue')
    # ax[0,0].tick_params(axis="y", colors="blue")
    # another.tick_params(axis="y", colors="red")
    # ax[0,0].set_xlabel("Average genotype per group")
    # ax[0,0].set_ylabel("Evolved production rate")
    # another.set_ylabel("Evolved signal Threshold")
    # ax[0,0].set_title("Evolved production rate (red) and Signal Threshold(red) vs genotype")

    # ax[1,0].scatter(prod, threshold, color= scalarMap.to_rgba(costs))
    # ax[1,0].set_xlabel("evolved production rate")
    # ax[1,0].set_ylabel("evolved siginaling threshold")
    plt.plot(range(max_G), np.full((max_G,), 50016 ), color="black", linewidth=3)
    plt.colorbar(scalarMap, ax= plt.gca(),label="Cost of Signaling $(10^8 FU)$")
    plt.ylim([0,110000])
    plt.xlabel("Generations", fontsize= 20)
    plt.ylabel('$\\frac{S_{Th} u}{p}$', fontsize= 25)
    plt.title("$\\frac{S_{Th} u}{p}$ Compared to $N_Th$ vs. Generations", fontsize= 25)
    # plt.tight_layout()
    # plt.subplots_adjust(left= .05, wspace=0.09, hspace=.524)
    plt.show()

def plotSupovernth(split_value, folder_name, minimum, maximum, step):
    cm = plt.get_cmap("plasma")
    cNorm = colors.Normalize(vmin=minimum, vmax= maximum)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

    max_G= 5000
    prod = []
    totalprod =[]
    threshold  = []
    costs = []
    
    for path, directories, files in os.walk(folder_name):
        for file in files:
            # cost = int(file.split(" ")[split_value])
            cost = np.round(float(file.split(" ")[split_value]), decimals=1) // 10**8
            filename = path+"\\"+file
            with open(filename, "r") as f:
                data = json.load(f)
            max_G = len(data["fit_Evo"])
            plt.plot(range(max_G), np.array(data["sig_Th_Evo"]) /np.array(data["pro_Rate_Evo"]) *  10.0 ** -4, color=scalarMap.to_rgba(cost))
    plt.plot(range(max_G), np.full((max_G,), 50016 ), color="black", linewidth=3)
    plt.colorbar(scalarMap, ax= plt.gca(),label="Cost of Signaling $(10^8 FU)$")
    plt.ylim([0,110000])
    plt.xlabel("Generations", fontsize= 20)
    plt.ylabel('$\\frac{S_{Th} u}{p}$', fontsize= 25)
    plt.title("$\\frac{S_{Th} u}{p}$ Compared to $N_Th$ vs. Generations", fontsize= 25)
    # plt.tight_layout()
    # plt.subplots_adjust(left= .05, wspace=0.09, hspace=.524)
    plt.show()


def optimum(p,r,s):
    a =  1.00450376e+03 * np.exp(-5.04190644e-02 / s) - 1.00444712e+03
    b =  2.0118419e-09 * s + 9.9727154e-08
    return 100 - 10** 8*np.abs(a *p +b - (p*(1+r)))

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

    name = "sig_Th"
    bins = np.digitize(data[name], bins=np.histogram_bin_edges(data[name],bins=binNo))
    bincolor= []
    x = np.histogram_bin_edges(data[name],bins=binNo)[:-1]
    counts =[]
    for i in range(binNo):
        counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
        bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    ax[1,1].bar(x, counts, color= bincolor, width=(np.max(data[name])-np.min(data[name]))/110)
    ax[1,1].set_title("sig_Th historgram")


    ax[2,1].scatter(data["pro_Rate"], data["sig_Th"],color=scalarMap.to_rgba(fit_Pop), s=3)
    ax[2,1].set_title("production rate vs sensitivity")

    data["ratio"] = np.array(data["sig_Th"])/ np.array(data["pro_Rate"]) * 10 ** -4
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


def graph_last_gen_standard(file):
    generation = file.split("\\gen ")[1].split(" ")[0]
 
    with open(file, "r") as f:
        data = json.load(f)
    fig, ax = plt.subplots(1,3, figsize=(18, 6))
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
    ax[0].bar(x, counts, color= bincolor, width=(np.max(data[name])-np.min(data[name]))/110)
    ax[0].set_title("Fitness Histogram")
    ax[0].set_xlabel("Fitness (FU)")

    difference = np.array(data["coopPayoff_Pop"]) - np.array(data["coopCost_Pop"])
    bins = np.digitize(difference, bins=np.histogram_bin_edges(difference,bins=binNo))
    bincolor= []
    x = np.histogram_bin_edges(difference,bins=binNo)[:-1]
    counts =[]
    for i in range(binNo):
        counts.append(np.count_nonzero([fit_Pop[j] for j in range(5000) if i == bins[j]]))
        bincolor.append(scalarMap.to_rgba(np.mean([fit_Pop[j] for j in range(5000) if i == bins[j]])))
    if (np.max(difference)-np.min(difference))/100 ==0:
        width = .05
    else: 
        width=(np.max(difference)-np.min(difference))/100
    ax[1].bar(x, counts, color= bincolor, width=width)
    ax[1].set_title("Cooperation Payoff vs. Cost Difference")
    ax[1].set_xlabel("Cooperation Payoff vs. Cost Difference (FU)")


    production = np.array(data["pro_Rate"]) + np.array(data["auto_pro_Rate"])
    extreme = np.array([np.min(production), np.max(production)])
    ax[2].plot(extreme, extreme * 50016 / 10 ** -4, color="black", linewidth = 2,  linestyle="dashed", dashes=(5, 10))


        
    ax[2].scatter(production, data["sig_Th"],color=scalarMap.to_rgba(fit_Pop), s=10)
    ax[2].set_title("Production Rate vs. Signaling Threshold")
    ax[2].set_xlabel("Production Rate $\left(\\frac{\mu M}{s}\\right)$")
    ax[2].set_ylabel("Signaling Threshold $\left(\mu M\\right)$")
    
    fig.suptitle(f"Generation {generation}", fontsize = 20)
    # ax[2],plot()
   

    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.162, hspace=.524)
    plt.show()


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
    if split_value == 6:
        splitname = "Cost of Signaling $(10^8 FU)$"
        fig.suptitle("Variation of Cost of Singaling", fontsize= 20)
    if split_value == 7:
        splitname = "Number of Genotypes per Testing Environment"
        fig.suptitle("Genetic Mixing", fontsize= 20)
    costs = []
    
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
            s = np.array(data["sig_Th_Evo"])

            ax[0].plot(range(max_G), s * 10 **-4 / (p+r), color=scalarMap.to_rgba(cost))
            ax[1].plot(range(max_G), np.array(data["pro_Rate_Evo"]) + np.array(data["auto_pro_Rate_Evo"]), color=scalarMap.to_rgba(cost))
            ax[2].plot(range(max_G), data["sig_Th_Evo"], color=scalarMap.to_rgba(cost))
    ax[0].plot(range(max_G), np.full((max_G,), 50016 ), color="black", linewidth=3)
    plt.colorbar(scalarMap, ax=ax[0], label= splitname)

    ax[0].set_title("$\\frac{S_{Th} u}{p_{total}}$ Compared to $N_Th$ Evolution")
    ax[0].set_xlabel("Generation")
    ax[0].set_ylabel('$\\frac{S_{Th} u}{p}$  $\left(\mu L\\right)$')

    ax[1].set_title("Total Production Rate Evolution")
    ax[1].set_xlabel("Generation")
    ax[1].set_ylabel('Production Rate $\left(\\frac{\mu M}{s}\\right)$')

    ax[2].set_title("Signaling Threshold Evolution")
    ax[2].set_xlabel("Generation")
    ax[2].set_ylabel('Signaling Threshold $\left(\mu M\\right)$')

    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.162, hspace=.524)
    plt.show()

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
    if split_value == 6:
        splitname = "Cost of Signaling $(10^8 FU)$"
        fig.suptitle("Variation of Cost of Singaling", fontsize= 20)
    if split_value == 7:
        splitname = "Number of Genotypes per Testing Environment"
        fig.suptitle("Genetic Mixing", fontsize= 20)
    costs = []
    
    for path, directories, files in os.walk(folder_name):
        for file in files:
            # cost = int(file.split(" ")[split_value])
            cost = np.round(float(file.split(" ")[split_value]), decimals=1) 
            filename = path+"\\"+file
            with open(filename, "r") as f:
                data = json.load(f)
            max_G = len(data["fit_Evo"])
           
            

            ax[0].plot(range(max_G), data["fit_Evo"], color=scalarMap.to_rgba(cost))
            ax[1].plot(range(max_G), data["coopPayoff_Evo"], color=scalarMap.to_rgba(cost))
            ax[2].plot(range(max_G), data["coopCost_Evo"], color=scalarMap.to_rgba(cost))

    plt.colorbar(scalarMap, ax=ax[0], label= splitname)

    ax[0].set_title("Fitness Evolution")
    ax[0].set_xlabel("Generation")
    ax[0].set_ylabel('Fitness (FU)')

    ax[1].set_title("Cooperation Payoff Evolution")
    ax[1].set_xlabel("Generation")
    ax[1].set_ylabel('Cooperation Payoff (FU)')

    ax[2].set_title("Cooperation Cost Evolution")
    ax[2].set_xlabel("Generation")
    ax[2].set_ylabel('Cooperation Cost (FU)')

    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.162, hspace=.524)
    plt.show()


if __name__ == "__main__":
    plt.rcParams["savefig.directory"]  = "Pictures"
    # graph("Wang python\json\Mon Jun 24 11-53-59 2024 0.5 1000000000 4 0.0001 True 5000.json")
    # graph_multiple_standard1(7, "Wang python\json\genotype with auto",0, 10, .1)
    for path, directories, files in os.walk("Wang python\json\\7 evo"):
        for file in files:
            graph_last_gen_standard(path+"\\"+file)
    # graph_last_gen_standard("Wang python\json\gen 1500 7 09-07 15-14-41 1000000000 5000 False False.json")
    # graph_last_gen_standard("Wang python\json\gen 2000 7 09-07 15-17-27 1000000000 5000 False False.json")
    # graph_last_gen_standard("Wang python\json\gen 2500 7 09-07 15-20-03 1000000000 5000 False False.json")
    # graph_last_gen_standard("Wang python\json\gen 3000 7 09-07 15-22-40 1000000000 5000 False False.json")