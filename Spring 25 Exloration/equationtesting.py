import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import scipy.optimize as optimize
import json

def s_star(N,m,p,r,u, K):
    square = np.sqrt((K*(u+m)-N*p*r-N*p)**2+4*K*N*p*(u+m))
    return (square - K*(u+m)+N*p*r+N*p) / (2* (u+m))

def x_star(N, m,S,f, P):
    s_adj = S / (S+500)
    return (N*P*s_adj)/(m+f)

def y_star(N, m, q, X, c, e):
    x_ad = X / (X+500)
    return (q*x_ad)/(c*N+m+e)

def fitness(p,r,u,K,f,q,c,e, P):
    s_sum = 0
    x_sum = 0
    y_sum = 0
    for m in np.arange(1.5e-7, 1.5e-4, step=100):
        for N in np.arange(10**1.5, 1e5, step=100):
            S = s_star(N,m,p,r,u, K)
            X = x_star(N, m,S,f, P)
            Y = y_star(N, m, q, X, c, e)
            s_sum += S
            y_sum += Y
            x_sum += X
    return ( 100 + coop_Benefit * y_sum - coop_cost * x_sum - sig_Cost * s_sum)


coop_Benefit = 0.0015
coop_cost = 0.005
sig_Cost = 0.00032
K = 50
p = 1.027e-07
r = 5
u = 1e-4
f = 4e-6
e = 4e-6
q = 2e0
P = 2.997e-08
c = 1e-8
N = 50015.81138830084
m =  1E-05

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
# fig, ax = plt.subplots(1,2)
# ax[0].plot(N_space, sig_Cost * s_star(N_space,m,p,r,u,K), label="S Star")
# ax[0].plot(N_space, coop_cost * x_star(N_space, m,s_star(N_space,m,p,r,u,K),f, P) ,  label="X Star")
# ax[0].plot(N_space, coop_Benefit * y_star(N_space, m, q, x_star(N_space, m,s_star(N_space,m,p,r,u,K),f, P), c, e),  label="Y Star")

# ax[0].legend()
# ax[0].set_xlabel("Cellular density")
# ax[1].plot(m_space, sig_Cost * s_star(N,m_space,p,r,u,K), label="S Star")
# ax[1].plot(m_space, coop_cost * x_star(N, m_space,s_star(N,m_space,p,r,u,K),f, P),  label="X Star")
# ax[1].plot(m_space, coop_Benefit * y_star(N, m_space, q, x_star(N, m_space,s_star(N,m_space,p,r,u,K),f, P), c, e),  label="Y Star")
# ax[1].legend()
# ax[1].set_xlabel("Mass transfer")
# plt.show() 


# ax[0].plot(N_space, coop_Benefit * y_star(N_space, m, q, x_star(N_space, m,s_star(N_space,m,p,r,u,K),f, P), c, e)- coop_cost * x_star(N_space, m,s_star(N_space,m,p,r,u,K),f, P)-sig_Cost * s_star(N_space,m,p,r,u,K),  label="fitness")
# ax[1].plot(m_space, coop_Benefit * y_star(N, m_space, q, x_star(N, m_space,s_star(N,m_space,p,r,u,K),f, P), c, e)-coop_cost * x_star(N, m_space,s_star(N,m_space,p,r,u,K),f, P) - sig_Cost * s_star(N,m_space,p,r,u,K),  label="fitness")
# ax[0].legend()
# ax[0].set_xlabel("Cellular density")
# ax[1].legend()
# ax[1].set_xlabel("Mass transfer")
# plt.show() 

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# x = []
# y = []
# z = []
# z2 = []

# for p in np.linspace(0.0,1e-7, 50):
#     for P in np.linspace(0.0,5e-7, 500):
#         fit = fitness(p,r,u,K,f,q,c,e, P)
#         if fit > 0:
#             x.append(p)
#             y.append(P)
#             z.append(fit)
# for p in np.linspace(1e-7,5e-7, 500):
#     for P in np.linspace(0.0,1e-7, 50):
#         fit = fitness(p,r,u,K,f,q,c,e, P)
#         if fit > 0:
#             x.append(p)
#             y.append(P)
#             z.append(fit)
# cm = plt.get_cmap("plasma")
# z = np.array(z).clip(0)
# cNorm = colors.Normalize(vmin=0, vmax= np.max(z))
# scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
# ax.scatter(x,y,z, color=scalarMap.to_rgba(z))
# ax.set_xlabel("S production")
# ax.set_ylabel("X production")
# ax.set_zlabel("Fiteness")
# plt.show()
# data = {"p":x, "P":y, "fit":z}
# with open("allPp500kx", "w") as f:
#     json.dump(data, f,  ensure_ascii=False, indent=4)

# PpDict = {"P":[], "p":[], "fit":[]}
# for P in np.logspace(0, 1e-5, 100):
#     opt = optimize.minimize(fitness, x0=.5e-7, args=(r,u,K,f,q,c,e, P), method="Nelder-Mead")
#     PpDict["P"].append(P)
#     PpDict["p"].append(list(opt["x"]))
#     PpDict["fit"].append(-opt["fun"])

# with open("Ppdata_extended", "w") as f:
#     json.dump(PpDict, f,  ensure_ascii=False, indent=4)

# with open("PpData", "r") as f:
#     data = json.load(f)
# with open("PpData_extended", "r") as f:
#     data2 = json.load(f)
#     data["p"] += data2["p"]
#     data["P"] += data2["P"] 
#     data["fit"] += data2["fit"] 
# cm = plt.get_cmap("winter")
# cNorm = colors.Normalize(vmin=100, vmax= 150)
# scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
# plt.scatter(np.log(data["p"][1:]), np.log(data["P"][1:]), color=scalarMap.to_rgba(data["fit"][1:]))
# plt.xlabel("Signal Production (Log)")
# plt.ylabel("X Production (Log)")
# plt.title("Optimal Signal vs X Production")
# plt.show()

# def fitness2d(Product,r,u,K,f,q,c,e):
#     p = Product[0]
#     P= Product[1]
#     s_sum = 0
#     x_sum = 0
#     y_sum = 0
#     for m in np.arange(1.5e-7, 1.5e-4, step=100):
#         for N in np.arange(10**1.5, 1e5, step=100):
#             S = s_star(N,m,p,r,u, K)
#             X = x_star(N, m,S,f, P)
#             Y = y_star(N, m, q, X, c, e)
#             s_sum += S
#             y_sum += Y
#             x_sum += X
#     return -(coop_Benefit * y_sum - coop_cost * x_sum - sig_Cost * s_sum)

# print(optimize.minimize(fitness2d, (1e-7, .5e-7), args=(r,u,K,f,q,c,e), method="Nelder-Mead"))

# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(projection='3d')
x = []
y = []
z = []

# for N in np.linspace(10**1.5, 1e5, 100):
#     for m in np.linspace(1.5e-7, 1.5e-4, 100):
#         fit = fitness(p,r,u,K,f,q,c,e, P)
#         if fit > 0:
#             x.append(N)
#             y.append(m)
#             z.append(coop_Benefit * y_star(N, m, q, x_star(N, m,s_star(N,m,p,r,u,K),f, P), c, e)-coop_cost * x_star(N, m,s_star(N,m,p,r,u,K),f, P) - sig_Cost * s_star(N,m,p,r,u,K))

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
# data = {"N":list(x), "m":list(y), "fit":list(z)}


fig, ax = plt.subplots(1,2)
s_space = np.linspace(0,10000,100)
x_space = np.linspace(0,10000,100)
ax[0].plot(s_space, coop_cost * x_star(N, m,s_space,f, P) ,  label="X Star", color="orange")
ax[0].plot(s_space, coop_Benefit * y_star(N, m, q, x_star(N, m,s_space,f, P), c, e),  label="Y Star", color="green")

ax[0].legend()
ax[0].set_xlabel("S Star")
ax[1].plot(x_space, coop_Benefit * y_star(N, m, q, x_space, c, e),  label="Y Star", color="green")
ax[1].legend()
ax[1].set_xlabel("X Star")
plt.show() 