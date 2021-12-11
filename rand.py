import math
import random
import numpy as np
import matplotlib.pyplot as plt


def get_array(mu=0, sigma=1, number_el=int(10e5)):
    return np.random.normal(mu, sigma, number_el)


def momentum(x, mu, sigma):
    tmp1 = (x - mu)
    tmp2 = tmp1**2
    tmp3 = tmp1**3
    tmp4 = tmp1**4
    print(0, " : ", tmp1.mean())
    print(sigma**2, " : ", tmp2.mean())
    print(0, " : ", tmp3.mean())
    print(3 * sigma**4, " : ", tmp4.mean(), "\n")


x1 = get_array()
momentum(x1, 0, 1)

mu1 = 1
sigma1 = 3
mu2 = 2
sigma2 = 2

y1 = get_array(mu1, sigma1)
y2 = get_array(mu2, sigma2)

y = y1 + y2
momentum(y, mu1 + mu2, math.sqrt(sigma1**2 + sigma2**2))
