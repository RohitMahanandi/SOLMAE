import numpy as np
#from rng import ChaCha20
from os import urandom
import matplotlib.pyplot as plt

def Unifcrown(R_min, R_max, randombytes=urandom):
    u_rho = int.from_bytes(randombytes(8), 'little')
    u_theta = int.from_bytes(randombytes(8), 'little')
    u_rho = (u_rho & 0x1fffffffffffff) * 2**(-53)
    u_theta = (u_theta & 0x1fffffffffffff) * 2**(-53)
    
    rho = np.sqrt(R_min**2+u_rho*(R_max**2-R_min**2))
    #It differs a bit from the written spec as far as I understand: instead of sampling u_theta in [0,1] and use cos(2pi*u_theta), it samples instead in [0, 1/4].
    #This way, one does not need to take absolute values in the next step of the algorithm. by Alex 2023.2.1.
    x = rho*np.cos(np.pi/2*u_theta)
    y = rho*np.sin(np.pi/2*u_theta) 
    # This is equivalent to Algorithm 9
    return x, y

if __name__ == '__main__':
    x_list = []
    y_list = []
    for _ in range(2):
        x, y = Unifcrown(3,1)
        x_list.append(x)
        y_list.append(y)
    print(x_list);
    print(y_list);
    plt.plot(x_list,y_list, 'o', markersize=3)

    plt.show()