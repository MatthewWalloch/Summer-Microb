import scipy
import numpy as np

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
    benifit = 1.5
    # initial testing environments
    env_CellDen = np.array(list(np.linspace(min_CellDen, max_CellDen,num=grid_Size)))

    cost_sum = np.sum([c*j for j in env_CellDen])
    benifit_sum = np.sum([int(c*j**2>median_CellDen) for j in env_CellDen])
    return (benifit*benifit_sum - (2-benifit)*cost_sum)

# res = scipy.optimize.minimize_scalar(fitness_sum, bounds=[0,5])
# print(res.x)
# print(res.fun)
# print(res.success)
# print(4.998*10**-6/(0.5e-08/10.0 ** -4))