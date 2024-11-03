import scipy
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import collections
import json
def fitness_sum_threshold(c):
    grid_Size = 100
    max_CellDen = 10.0 ** 5
    # minimum cellular density
    min_CellDen = 10.0 ** 1.5
    # median cellular density
    median_CellDen = (max_CellDen + min_CellDen) * .5
    # baseline volume
    base_Volume = 10.0
    # signal decay rate
    decay_Rate = 10.0 ** -4.
    benifit = 1.5
    # initial testing environments
    env_CellDen = np.array(list(np.linspace(min_CellDen, max_CellDen,num=grid_Size)))

    cost_sum = np.sum([c*j for j in env_CellDen])
    benifit_sum = np.sum([int(c*j**2>median_CellDen) for j in env_CellDen])
    return (benifit*benifit_sum - (2-benifit)*cost_sum)

def fitness_sum(c):
    grid_Size = 100
    max_CellDen = 10.0 ** 5
    # minimum cellular density
    min_CellDen = 10.0 ** 1.5
    # median cellular density
    median_CellDen = (max_CellDen + min_CellDen) * .5
    # baseline volume
    base_Volume = 10.0
    # signal decay rate
    decay_Rate = 10.0 ** -4.
    benifit = 1.1
    cost = .5
    # initial testing environments
    env_CellDen = np.array(list(np.linspace(min_CellDen, max_CellDen,num=grid_Size)))

    cost_sum = np.sum([c*j for j in env_CellDen])
    benifit_sum = np.sum([np.log(c*j**2/median_CellDen + 1) for j in env_CellDen])
    return (benifit*benifit_sum - cost*cost_sum)


def s_star(N, p, r):
    u = 10.0 ** -4
    K=50
    square = np.sqrt((K*u-N*p*r-N*p)**2+4*K*N*p*u)
    return (square - K*u+N*p*r+N*p) / (2* u)

def fitness_auto(p,r, s):
    grid_Size = 100
    max_CellDen = 10.0 ** 5
    # minimum cellular density
    min_CellDen = 10.0 ** 1.5
    # median cellular density
    median_CellDen = (max_CellDen + min_CellDen) * .5
    # baseline volume
    base_Volume = 10.0
    # signal decay rate
    u = 10.0 ** -4.
    K=50
    benifit = 1.5
    # initial testing environments
    env_CellDen = np.array(list(np.linspace(min_CellDen, max_CellDen,num=grid_Size)))

    cost_sum = np.sum([int(s_star(j,p,r) > s) for j in env_CellDen])
    benifit_sum = np.sum([int(j * int(s_star(j,p,r) > s) >median_CellDen) for j in env_CellDen])
    return (benifit*benifit_sum - (2-benifit)*cost_sum)

def fitness_sum_all(Sn, p):
    grid_Size = 100
    max_CellDen = 10.0 ** 5
    # minimum cellular density
    min_CellDen = 10.0 ** 1.5
    # median cellular density
    median_CellDen = (max_CellDen + min_CellDen) * .5
    # baseline volume
    base_Volume = 10.0
    # signal decay rate
    decay_Rate = 10.0 ** -4.
    benifit = 1.5
    cost = .5
    # initial testing environments
    env_CellDen = np.array(list(np.linspace(min_CellDen, max_CellDen,num=grid_Size)))
    sig_cost =  10**9
    cost_sum = np.sum([Sn*j*p / decay_Rate for j in env_CellDen])
    benifit_sum = Sn*np.sum([np.log(.5 * Sn*j**2*p / decay_Rate /median_CellDen + 1) for j in env_CellDen])
    return 100+benifit*benifit_sum - cost*cost_sum-sig_cost*p


# print(fitness_auto(5*10**-9,5,5))
# res = scipy.optimize.minimize_scalar(fitness_sum, bounds=[0,5])
# print(res.x)
# print(res.fun)
# print(res.success)
# print(4.998*10**-6/(0.5e-08/10.0 ** -4))

# for path, directories, files in os.walk("New python\\auto json\generations for 10"):
#         for file in files:
#             # print(file.split("Sat")[0])

#             if file.split(" ")[3] == "11-07":
#                 os.rename("New python\\auto json\generations for 10\\"+file, "New python\\auto json\generations for 10\\"+file.split(" 11-07")[0]+".json")
#             if file.split(" ")[3] == "12-07":
#                 os.rename("New python\\auto json\generations for 10\\"+file, "New python\\auto json\generations for 10\\"+file.split(" 12-07")[0]+".json")

def optimum(p,r,s):
    a =  1.00450376e+03 * np.exp(-5.04190644e-02 / s) - 1.00444712e+03
    b =  2.0118419e-09 * s + 9.9727154e-08
    return np.abs(a *p +b - (p*(1+r)))


def plot_no_change():
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    x = []
    y = []
    z = []
    z2 = []
    for p in np.linspace(0.0,10*10**-9, 100):
        for r in np.linspace(0, 20, 100):
            x.append(p)
            y.append(r)
            # x.append(p)
            # y.append(r)
            # x.append(p)
            # y.append(r)
            z.append(fitness_auto(p,r,10))
            # z.append(fitness_auto(p,r,2))
            # z.append(fitness_auto(p,r,10))
            z2.append(optimum(p,r,10))
            # z2.append(optimum(p,r,2))
            # z2.append(optimum(p,r,10))
    cm = plt.get_cmap("plasma_r")
    cNorm = colors.Normalize(vmin=np.min(z2), vmax= np.max(z2))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    ax.scatter(x,y,z, color=scalarMap.to_rgba(z2))
    # ax.scatter(x,y,z2, color=scalarMap.to_rgba(z2))
    plt.show()


# plot_no_change()
# s= "8.88888888888889"
# s= "4.040404040404041"
# s= "1.8181818181818181"
# cm = plt.get_cmap("plasma")
# cNorm = colors.Normalize(vmin=0, vmax= 20)
# scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
# with open("s value.json", "r") as f:
#         data = json.load(f)
# def func(x, a, b):
#     return a*x+b

# def exp(x, a, b, c):
#     return  a * np.exp(b / x) -c

# s_list = []
# slopelist = []
# for s in data.keys():
#     if float(s) > .5:
#         x = np.array(data[s]["p"])
#         y = np.array(data[s]["r"])
#         if len(y) > 0:
#             s_list.append(float(s))
#             # plt.plot(x, y,color="blue")
#             popt, pcov = scipy.optimize.curve_fit(func, x, x*(1+y))
#             # # print(x/y)
#             slopelist.append(popt[1])
#             # plt.plot(x, func(x, *popt), color=scalarMap.to_rgba(float(s)))
#             # plt.plot(x, , color=scalarMap.to_rgba(float(s)))
#             # print(popt)

# s_list= np.array(s_list)
# slopelist = np.array(slopelist)
# popt, pcov = scipy.optimize.curve_fit(func, s_list, slopelist)
# print(popt)
# plt.plot(s_list, func(s_list, *popt),color="blue")
# # plt.plot(x, , color=scalarMap.to_rgba(float(s)))
# plt.plot(s_list, slopelist, color="red")
# plt.show()

# 2.0118419e-09 * s + 9.9727154e-08
# 1.00450376e+03 * np.exp(-5.04190644e-02 * p *(1+r) / s) - 1.00444712e+03

# fig, ax = plt.subplots(1, figsize=(8, 6))
# steps = np.linspace(0, 6e-5, num=10000)
# y = [fitness_sum(c) for c in steps]
# minmax = [np.min(y), np.max(y) + 10]


# # plt.plot([4.992e-6, 4.992e-6], minmax,color ="black", linestyle="dashed")
# # plt.plot([1.620e-5, 1.62e-5], minmax,color ="black", linestyle="dashdot")

# plt.plot(steps, [fitness_sum(c) for c in steps])
# plt.plot(steps, [fitness_sum_threshold(c) for c in steps])
# plt.xlabel("$Sn_i S_{N_j}$",fontsize=20)
# plt.ylabel("Net Cooperation Impact",fontsize=20)


# plt.show() 

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# S = np.arange(0,2,.01)
# P = np.arange(0,10**-8, 10**-10)
# X , Y = np.meshgrid(S,P)
# zs= np.array([fitness_sum_all(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
# Z = zs.reshape(X.shape)
# # x = []
# # y = []
# # z = []
# # for s in S:
# #     for p in P:
# #         x.append(s)
# #         y.append(p)
# #         z.append(fitness_sum_all(s,p))

# ax.plot_surface(X,Y,Z, cmap=cmx.coolwarm, linewidth=0)
# ax.set_xlabel("Sensitivity")
# ax.set_ylabel("Production")
# ax.set_zlabel("fitness")
# plt.show()
