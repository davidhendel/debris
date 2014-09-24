# mcint/integrate.py

import itertools
import math
import numpy as np
from streamteam.util import get_pool

def integrate(integrand, sampler, args=(), measure=1.0, n=100):
    # Sum elements and elements squared
    total = 0.0
    total_sq = 0.0
    for x in itertools.islice(sampler, n):
        f = integrand(x, M=args)
        total += f
        total_sq += (f**2)
    # Return answers
    sample_mean = total/n
    sample_var = (total_sq - ((total/n)**2)/n)/(n-1.0)
    return (measure*sample_mean, measure*math.sqrt(sample_var/n))


def integrate_2(integrand, sampler, args=(), measure=1.0, n=100):
    # Sum elements and elements squared
    total_1_arr    = np.zeros(n)
    total_2_arr    = np.zeros(n)
    total_sq_1_arr = np.zeros(n)
    total_sq_2_arr = np.zeros(n)
    i = 0.
    for x in itertools.islice(sampler, n):
        total_1_arr[i], total_2_arr[i] = integrand(x, args)
        total_sq_1_arr[i], total_sq_2_arr[i] = total_1_arr[i]**2, total_2_arr[i]**2
        i = i+1
    # Return answer
    total_1 = np.sum(total_1_arr)
    total_2 = np.sum(total_2_arr)
    total_sq_1 = np.sum(total_sq_1_arr)
    total_sq_2 = np.sum(total_sq_2_arr)

    sample_mean_1 = total_1/n
    sample_mean_2 = total_2/n
    sample_var_1 = (total_sq_1 - ((total_1/n)**2)/n)/(n-1.0)
    sample_var_2 = (total_sq_2 - ((total_2/n)**2)/n)/(n-1.0)
    return (measure*sample_mean_1, measure*math.sqrt(sample_var_1/n),measure*sample_mean_2, measure*math.sqrt(sample_var_2/n))
