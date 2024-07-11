import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
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
if __name__ == "__main__":
    file = "New python\\auto json\\1719860836745544100 50 0.2 1000000000 5000 True.json"
    graph_auto(file)
    # graph_last_gen(file)
    # graph_multiple(2, "New python\json\\actuall 50 genotype evo", 0, 10, .5)
    


    # steps = np.linspace(0, 6e-5, num=10000)
    # plt.plot(steps, [testing.fitness_sum(c) for c in steps])
    # plt.xlabel("S*p/u")
    # plt.ylabel("Benifit minus cooperation")
    # plt.show()