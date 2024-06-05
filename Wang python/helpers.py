import numpy as np
index_Cheats = [0,1,0,0,0]
rand = [1,5,4,3,2]
index_non = np.ones(5)-index_Cheats 
print(np.dot(rand,index_non) / 5)
