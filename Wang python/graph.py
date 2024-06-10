import matplotlib.pyplot as plt
import json

filename = "Wang python\json\Mon Jun 10 15-52-34 2024 0.5 1000000000.0 2 0.0001 False 1000.json"
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

ax[3,0].plot(range(max_G), data["auto_R_Evo"])
ax[3,0].set_title("auto_R_Evo")

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