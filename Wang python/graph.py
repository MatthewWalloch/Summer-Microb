import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import json
import os
import numpy as np


# 

filename = "Wang python\json\Mon Jun 17 16-22-08 2024 0.5 1000000000 3 0.0001 False 5000.json"
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

ax[3,1].plot(range(max_G), data["auto_pro_Rate_Evo"])
ax[3,1].set_title("auto_pro_Rate_Evo")
plt.tight_layout()
plt.show()

def graph_multiple(split_value):
    cm = plt.get_cmap("plasma")
    cNorm = colors.Normalize(vmin=5, vmax=110)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig, ax = plt.subplots(4,2)
    max_G= 0


    prod = []
    threshold  = []
    costs = []
    for path, directories, files in os.walk("Wang python\json\Fig1"):
        for file in files:
            cost = int(file.split(" ")[6]) // 10** 8
            filename = path+"\\"+file
            with open(filename, "r") as f:
                data = json.load(f)
            max_G = len(data["fit_Evo"])
            ax[3,1].plot(range(max_G), data["fit_Evo"], color=scalarMap.to_rgba(cost))
            ax[3,1].set_title("fit_Evo vs generation")

            ax[2,0].plot(range(max_G), data["pro_Rate_Evo"], color=scalarMap.to_rgba(cost))
            ax[2,0].set_title("pro_Rate_Evo vs generation")

            ax[3,0].plot(range(max_G), data["sig_Th_Evo"], color=scalarMap.to_rgba(cost))
            ax[3,0].set_title("sig_Th_Evo vs generation")

            ax[0,1].plot(range(max_G), data["coopPayoff_Evo"], color=scalarMap.to_rgba(cost))
            ax[0,1].set_title("coopPayoff_Evo vs generation")

            ax[1,1].plot(range(max_G), data["sigCost_Evo"], color=scalarMap.to_rgba(cost))
            ax[1,1].set_title("sigCost_Evo vs generation")

            ax[2,1].plot(range(max_G), data["coopCost_Evo"], color=scalarMap.to_rgba(cost))
            ax[2,1].set_title("coopCost_Evo vs generation")

            prod.append(np.mean(data["pro_Rate_Evo"][-50:]))
            threshold.append(np.mean(data["sig_Th_Evo"][-50:]))
            costs.append(cost)

        

        print(path)

    ax[0,0].bar([i for i in range(5,105, 5)], 
            [1 for i in range(5,105,5)], 
            width=5,
            tick_label=[i for i in range(5,105, 5)],
            color=[scalarMap.to_rgba(i) for i in range(5,105, 5)])
    ax[0,0].set_xlabel("Cost of signaling (10^8 FU)")

    ax[1,0].scatter(prod, threshold, color= scalarMap.to_rgba(costs))
    ax[1,0].set_xlabel("evolved production rate")
    ax[1,0].set_ylabel("evolved siginaling threshold")
    plt.tight_layout()
    plt.show()