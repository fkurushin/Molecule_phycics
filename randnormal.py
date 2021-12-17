import math
import random
import numpy as np


def get_array(mu=0, sigma=1, number_el=int(10e5)):
    return np.random.normal(mu, sigma, number_el)


# Здесь возможно лежит баг
def get_number(mu=0, sigma=1):
    return get_array(mu, sigma, number_el=1)[0]


def momentum(x, mu, sigma):
    tmp1 = (x - mu)
    tmp2 = tmp1**2
    tmp3 = tmp1**3
    tmp4 = tmp1**4
    print(0, " : ", tmp1.mean())
    print(sigma**2, " : ", tmp2.mean())
    print(0, " : ", tmp3.mean())
    print(3 * sigma**4, " : ", tmp4.mean(), "\n")
