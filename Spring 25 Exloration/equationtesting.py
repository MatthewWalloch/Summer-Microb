import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import scipy.optimize as optimize
import json

def s_star(N,m,p,r,u, K):
    square = np.sqrt((K*(u+m)-N*p*r-N*p)**2+4*K*N*p*(u+m))
    return (square - K*(u+m)+N*p*r+N*p) / (2* (u+m))

def x_star(N, m,S,f, P, kx= 500):
    s_adj = S / (S+kx)
    return (N*P*s_adj)/(m+f)

def y_star(N, m, q, X, c, e, ky= 500):
    x_ad = X / (X+ky)
    return (q*x_ad)/(c*N+m+e)

def fitness(p,r,u,K,f,q,c,e, P,kx=500,ky=500):
    s_sum = 0
    x_sum = 0
    y_sum = 0
    for m in np.linspace(1.5e-7, 1.5e-4, num=100):
        for N in np.linspace(10**1.5, 1e5, num=100):
            S = s_star(N,m,p,r,u, K)
            X = x_star(N, m,S,f, P, kx=kx)
            Y = y_star(N, m, q, X, c, e,ky=ky)
            s_sum += S
            y_sum += Y
            x_sum += X
    return -( 100 + coop_Benefit * y_sum - coop_cost * x_sum - sig_Cost * s_sum)


coop_Benefit = 0.0015
coop_cost = 0.005
sig_Cost = 0.00032
K = 50
p = 1e-08
r = 5
u = 1e-4
f = 4e-6
e = 4e-6
q = 2e0
P = 1e6
c =  1e-8
N = 50015.81138830084
m =  1E-05
p= 3.686e-08
P= 1.780e-07

p_val = [0, 2e-08]
r_val = [0,20]
u_val = [10e-11,10]
N_val = [0,5e4]  #[10**1.5, 1e5]
m_val = [0, 3e-4]  # [1.5e-7, 1.5e-4]

N_space = np.linspace(N_val[0],N_val[1], 100)
m_space = np.linspace(m_val[0],m_val[1], 100)
p_space = np.linspace(p_val[0],p_val[1], 100)
r_space = np.linspace(r_val[0],r_val[1], 100)
u_space = np.linspace(u_val[0],u_val[1], 100)

kx=400
ky=600

fig, ax = plt.subplots(1,2, figsize=(12,6))
ax[0].plot(N_space, sig_Cost * s_star(N_space,m,p,r,u,K), label="S Star")
ax[0].plot(N_space, coop_cost * x_star(N_space, m,s_star(N_space,m,p,r,u,K),f, P,kx=kx) ,  label="X Star")
ax[0].plot(N_space, coop_Benefit * y_star(N_space, m, q, x_star(N_space, m,s_star(N_space,m,p,r,u,K),f, P,kx=kx), c, e,ky=ky),  label="Y Star")

ax[0].legend()
ax[0].set_xlabel("Cellular density")
ax[1].plot(m_space, sig_Cost * s_star(N,m_space,p,r,u,K), label="S Star")
ax[1].plot(m_space, coop_cost * x_star(N, m_space,s_star(N,m_space,p,r,u,K),f, P,kx=kx),  label="X Star")
ax[1].plot(m_space, coop_Benefit * y_star(N, m_space, q, x_star(N, m_space,s_star(N,m_space,p,r,u,K),f, P,kx=kx), c, e,ky=ky),  label="Y Star")
ax[1].legend()
ax[1].set_xlabel("Mass transfer")
ax[0].set_ylim(0,1.1)
ax[1].set_ylim(0,2.9)
plt.show() 

# fig, ax = plt.subplots(1,2)
# ax[0].plot(N_space, coop_Benefit * y_star(N_space, m, q, x_star(N_space, m,s_star(N_space,m,p,r,u,K),f, P), c, e)- coop_cost * x_star(N_space, m,s_star(N_space,m,p,r,u,K),f, P)-sig_Cost * s_star(N_space,m,p,r,u,K),  label="fitness")
# ax[1].plot(m_space, coop_Benefit * y_star(N, m_space, q, x_star(N, m_space,s_star(N,m_space,p,r,u,K),f, P), c, e)-coop_cost * x_star(N, m_space,s_star(N,m_space,p,r,u,K),f, P) - sig_Cost * s_star(N,m_space,p,r,u,K),  label="fitness")
# ax[0].legend()
# ax[0].set_xlabel("Cellular density")
# ax[1].legend()
# ax[1].set_xlabel("Mass transfer")
# plt.show() 


# s_space = np.linspace(0,1600, 1000)
# plt.plot(s_space,  x_star(N, m,s_space,f,P,kx=kx), label="X_star")
# plt.plot(s_space, y_star(N, m, q, x_star(N, m,s_space,f,P,kx=kx), c, e,ky=ky),  label="Y_star")
# plt.legend()
# plt.ylim(0,2500)
# plt.xlabel("S_star")
# plt.show() 

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# x = []
# y = []
# z = []
# z2 = []

# for p in np.linspace(0.0,5e-7, 50):
#     for P in np.linspace(0.0,5e-7, 50):
#         fit = -fitness(p,r,u,K,f,q,c,e, P,ks=400, kx=800)
#         if fit > 0:
#             x.append(p)
#             y.append(P)
#             z.append(fit)
# # # for p in np.linspace(1e-7,5e-7, 100):
# # #     for P in np.linspace(0.0,1e-7, 50):
# # #         fit = fitness(p,r,u,K,f,q,c,e, P)
# # #         if fit > 0:
# # #             x.append(p)
# # #             y.append(P)
# # #             z.append(fit)
# cm = plt.get_cmap("plasma")
# z = np.array(z).clip(0)
# cNorm = colors.Normalize(vmin=0, vmax= np.max(z))
# scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
# ax.scatter(x,y,z, color=scalarMap.to_rgba(z))
# ax.set_xlabel("S production")
# ax.set_ylabel("X production")
# ax.set_zlabel("Fiteness")
# plt.show()
# data = {"p":list(x), "P":list(y), "fit":list(z)}
# with open("allPp500kx", "w") as f:
#     json.dump(data, f,  ensure_ascii=False, indent=4)

# PpDict = {"P":[], "p":[], "fit":[]}
# for P in np.linspace(0, 1e-5, 400):
#     print(P)
#     opt = optimize.minimize(fitness, x0=.5e-7, args=(r,u,K,f,q,c,e, P), method="Nelder-Mead")
#     PpDict["P"].append(P)
#     PpDict["p"].append(list(opt["x"]))
#     PpDict["fit"].append(-opt["fun"])

# with open("Ppdata_extended", "w") as f:
#     json.dump(PpDict, f,  ensure_ascii=False, indent=4)

# with open("PpData", "r") as f:
#     data = json.load(f)
# with open("PpData_extended", "r") as f:
#     data = json.load(f)
#     # data["p"] += data2["p"]
#     # data["P"] += data2["P"] 
#     # data["fit"] += data2["fit"] 
# cm = plt.get_cmap("winter")
# cNorm = colors.Normalize(vmin=200, vmax= 400)
# scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
# plt.scatter(np.log(data["p"][0:]), np.log(data["P"][0:]), color=scalarMap.to_rgba(data["fit"][0:]))
# plt.xlabel("Signal Production (Log)")

# plt.ylabel("X Production (Log)")
# plt.title("Optimal Signal vs X Production")
# plt.show()

def fitness2d(Product,r,u,K,ks,kx,f,q,c,e):
    p = Product[0]
    P= Product[1]
    s_sum = 0
    x_sum = 0
    y_sum = 0
    for m in np.linspace(1.5e-7, 1.5e-4, num=100):
        for N in np.linspace(10**1.5, 1e5, num=100):
            S = s_star(N,m,p,r,u, K)
            X = x_star(N, m,S,f, P, ks=ks)
            Y = y_star(N, m, q, X, c, e,kx=kx)
            s_sum += S
            y_sum += Y
            x_sum += X
    return -(100 + coop_Benefit * y_sum - coop_cost * x_sum - sig_Cost * s_sum)

# halfkdict = {"ks":[], "kx":[], "p, eta":[],  "fit":[]}
# for ks in np.linspace(50, 1000, num=10):
#     print(ks)
#     for kx in np.linspace(50, 1000, num=10):
#         optimium = optimize.minimize(fitness2d, (1e-7, .5e-7), args=(r,u,K,ks, kx,f,q,c,e), method="Nelder-Mead")
#         halfkdict["ks"].append(ks)
#         halfkdict["kx"].append(kx)
#         halfkdict["p, eta"].append(list(optimium["x"]))
#         halfkdict["fit"].append(optimium["fun"])
# with open("Half concentration maxes", "w") as f:
    json.dump(halfkdict, f,  ensure_ascii=False, indent=4)

# with open("Half concentration maxes", "r") as f:
#     data = json.load(f)

# ks = np.array(data["ks"])
# kx = np.array(data["kx"])
# p = []
# eta = []
# for g in data["p, eta"]:
#     p.append(g[0])
#     eta.append(g[1])
# p = np.array(p)
# eta = np.array(eta)
# fit = -np.array(data["fit"])
# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(projection='3d')
# cm = plt.get_cmap("coolwarm")
# cNorm = colors.Normalize(vmin=np.min(fit), vmax= np.max(fit) )
# scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
# ax.scatter(ks,kx, p, color=scalarMap.to_rgba(fit))
# ax.set_xlabel("$K_x$")
# ax.set_ylabel("$K_y$")
# ax.set_zlabel("$p$")
# ax.set_zlim(0,.75e-6)
# plt.show()



# print(optimize.minimize(fitness2d, (1e-7, .5e-7), args=(r,u,K,f,q,c,e), method="Nelder-Mead"))

# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(projection='3d')
# x = []
# y = []
# z = []
# p= 3.686e-08
# P= 1.780e-07
# for N in np.linspace(10**1.5, 1e5, 100):
#     for m in np.linspace(1.5e-7, 1.5e-4, 100):
#         x.append(N)
#         y.append(m)
#         z.append(coop_Benefit * y_star(N, m, q, x_star(N, m,s_star(N,m,p,r,u,K),f, P, ks=ks), c, e,kx=kx)-coop_cost * x_star(N, m,s_star(N,m,p,r,u,K),f, P,ks=ks) - sig_Cost * s_star(N,m,p,r,u,K))
 
# # data = {"N":list(x), "m":list(y), "fit":list(z)}
# # with open("CeldenMassTranData", "w") as f:
# #     json.dump(data, f,  ensure_ascii=False, indent=4)

# with open("CeldenMassTranData", "r") as f:
#     data = json.load(f)
# x = data["N"]
# y = data["m"]
# z = data["fit"]
# cm = plt.get_cmap("coolwarm")
# # z = np.array(z).clip(0)
# cNorm = colors.Normalize(vmin=-.1, vmax= .1 )
# scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
# ax.scatter(x,y,z, color=scalarMap.to_rgba(z))
# ax.set_xlabel("Celular Density")
# ax.set_ylabel("Mass Transfer")
# ax.set_zlabel("Fiteness")
# plt.show()


# fig, ax = plt.subplots(1,2)
# s_space = np.linspace(0,10000,100)
# x_space = np.linspace(0,10000,100)
# ax[0].plot(s_space, coop_cost * x_star(N, m,s_space,f, P) ,  label="X Star", color="orange")
# ax[0].plot(s_space, coop_Benefit * y_star(N, m, q, x_star(N, m,s_space,f, P), c, e),  label="Y Star", color="green")

# ax[0].legend()
# ax[0].set_xlabel("S Star")
# ax[1].plot(x_space, coop_Benefit * y_star(N, m, q, x_space, c, e),  label="Y Star", color="green")
# ax[1].legend()
# ax[1].set_xlabel("X Star")
# plt.show() 