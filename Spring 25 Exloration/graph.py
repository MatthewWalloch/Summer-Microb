import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import font_manager
import matplotlib
import json
import os
import numpy as np


# what I used to graph things
# the layout may or may not make sense to you
def graph(filename):
    
    with open(filename, "r") as f:
        data = json.load(f)
    fig, ax = plt.subplots(4,2)
    max_G= len(data["fit_Pop_Evo"])
    fig.suptitle(f"$c=${filename.split('.j')[0].split(' ')[-1]}")

    
    ax[0,1].plot(range(max_G), data["fit_Pop_Evo"])
    ax[0,1].set_title("Fitness per generation")


    ax[1,1].plot(range(max_G), data["coopPayoff_Pop_Evo"])
    ax[1,1].set_title("Y sum per generation")

    ax[2,1].plot(range(max_G), data["coopCost_Pop_Evo"])
    ax[2,1].set_title("X sum per generation")

    ax[3,1].plot(range(max_G), data["sigCost_Pop_Evo"])
    ax[3,1].set_title("S sum per generation")

    ax[0,0].plot(range(max_G), data["X_pro_Rate_Evo"])
    ax[0,0].set_title("X Production Rate")

    ax[1,0].plot(range(max_G), data["pro_Rate_Evo1"])
    ax[1,0].set_title("pro_Rate_Evo")

    ax[2,0].plot(range(max_G), data["decay_Rate_Evo1"])
    ax[2,0].set_title("decay_Rate_Evo")

    ax[3,0].plot(range(max_G), data["induct_Rate_Evo1"])
    ax[3,0].set_title("induct_Rate_Evo")

    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.09, hspace=.524)
    plt.show()


def graph_last_gen(file):
    
    with open(file, "r") as f:
        data = json.load(f)
    fig, ax = plt.subplots(4,2, figsize=(18, 6))
    binNo = 100
    fit_Pop = data["fit_Pop"]
    cm = plt.get_cmap("plasma")
    q1, mead, q3 = np.quantile(fit_Pop, [.25,.5,.75])
    cNorm = colors.CenteredNorm(vcenter=mead, halfrange=(q3-q1))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig.suptitle(file.split("\\")[3].split(".")[0])
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


    ax[2,1].scatter(data["pro_Rate"], data["sensitivity"],color=scalarMap.to_rgba(fit_Pop), s=5)
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


if __name__ == "__main__":
    # if you want to graph a whole bunch of stuff
    for path, directories, files in os.walk("Spring 25 Exloration\json\Production rate testing"):
        for file in files:
            if file.split(" ")[0] in ["04-06", "04-07"]:
                print(file)
                graph(f"Spring 25 Exloration\json\\Production rate testing\\{file}")