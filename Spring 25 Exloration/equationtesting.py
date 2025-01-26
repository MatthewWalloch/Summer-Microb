import numpy as np
import matplotlib.pyplot as plt

def s_star(N,m,p,r,u, K):
    square = np.sqrt((K*(u+m)-N*p*r-N*p)**2+4*K*N*p*(u+m))
    return (square - K*(u+m)+N*p*r+N*p) / (2* (u+m))

def x_star(N, m,S,f):
    s_adj = S / (S+1.63)
    return (N*5e-9*s_adj)/(m+f)

def y_star(N, m, q, X, c, e):
    x_ad = X / (X+500)
    return (q*x_ad)/(c*N+m+e)

def fitness(p,r,u,K,f,q,c,e):
    s_sum = 0
    x_sum = 0
    y_sum = 0
    for m in np.arange(1.5e-7, 1.5e-4, step=100):
        for N in np.arange(10**1.5, 1e5, step=100):
            S = s_star(N,m,p,r,u, K)
            X = x_star(N, m,S,f)
            Y = y_star(N, m, q, X, c, e)
            s_sum += S
            y_sum += Y
            x_sum += X
    s_avg = s_sum / 50015.81138830084
    return 100 + 1e-2 * y_sum - 1e-5 * x_sum - 1e-10* (1+r*s_avg / (K +s_avg))*p

K = 50
p = 5e-9
r = 5
u = 1e-4
f = 4e-6
e = 4e-6
q = 1e0
c = 1e-8
N = 50015.81138830084
m =  7.5e-5

p_val = [0, 2e-08]
r_val = [0,20]
u_val = [10e-11,10]
N_val = [10**1.5, 1e5]
m_val = [1.5e-7, 1.5e-4]

N_space = np.linspace(N_val[0],N_val[1], 100)
m_space = np.linspace(m_val[0],m_val[1], 100)
p_space = np.linspace(p_val[0],p_val[1], 100)
r_space = np.linspace(r_val[0],r_val[1], 100)
u_space = np.linspace(u_val[0],u_val[1], 100)
fig, ax = plt.subplots(1,2)
# ax[0].plot(N_space, s_star(N_space,m,p,r,u,K), label="S Star")
ax[0].plot(N_space, x_star(N_space, m,s_star(N_space,m,p,r,u,K),f) +s_star(N_space,m,p,r,u,K) ,  label="X Star")
ax[0].plot(N_space, y_star(N_space, m, q, x_star(N_space, m,s_star(N_space,m,p,r,u,K),f), c, e),  label="Y Star")
ax[0].legend()
ax[0].set_xlabel("Cellular density")
# ax[1].plot(m_space, s_star(N,m_space,p,r,u,K), label="S Star")
ax[1].plot(m_space, s_star(N,m_space,p,r,u,K)+x_star(N, m_space,s_star(N,m_space,p,r,u,K),f),  label="X Star")
ax[1].plot(m_space, y_star(N, m_space, q, x_star(N, m_space,s_star(N,m_space,p,r,u,K),f), c, e),  label="Y Star")
ax[1].legend()
ax[1].set_xlabel("Mass transfer")
plt.show() 

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# x = []
# y = []
# z = []
# z2 = []
# for p in np.linspace(0.0,2e-08, 50):
#     for u in np.linspace(10e-11,10e-7, 50):
#         x.append(p)
#         y.append(u)
#         z.append(fitness(p,r,u,K,f,q,c,e))
# ax.scatter(x,y,z)
# plt.show()