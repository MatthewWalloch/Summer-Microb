import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.widgets import Slider
import matplotlib.cm as cmx
import json
import os
import numpy as np



def graph_last_gen(file):
    for i in range(4):
        for j in range(2):
            ax[i,j].clear()

    file = f"New python\json\\50 genotypes\gen {file*100} 49.0.json"
    with open(file, "r") as f:
        data = json.load(f)
    
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
    if np.max(data[name])-np.min(data[name]) == 0:
        ax[1,0].bar(x, counts, color= bincolor, width=.01)
        ax[1,0].set_title("coopPayoff_Pop historgram")
    else:
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

    
    return fig, ax
    


if __name__ == "__main__":

    # graph(file)
    fig, ax = plt.subplots(4,2)
    axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    freq_slider = Slider(
        ax=axfreq,
        label='generation in hundres',
        valmin=5,
        valmax=50,
        valstep=5,
        valinit=5
    )

    
    plt.tight_layout()
    plt.subplots_adjust(left= .05, wspace=0.09, hspace=.524, bottom=0.185)
    freq_slider.on_changed(graph_last_gen)
    plt.show()
    
    # # graph_multiple(7, "Wang python\json\Evolve genotype", 0, 10, .1)
    


    # steps = np.linspace(0, 6e-5, num=10000)
    # plt.plot(steps, [testing.fitness_sum(c) for c in steps])
    # plt.xlabel("S*p/u")
    # plt.ylabel("Benifit minus cooperation")
    # plt.show()