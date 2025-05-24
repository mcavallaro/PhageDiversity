# Code from
# Johanna Elena Schmitz, Sven Rahmann, A comprehensive review and evaluation of species richness estimation, Briefings in Bioinformatics, Volume 26, Issue 2, March 2025, bbaf158, https://doi.org/10.1093/bib/bbaf158


import argparse
import math
import numpy as np
from scipy.stats import poisson, binom


def GT_BI(freq, n, t):
    """
    EFRON-THISTED
    """
    r = len(freq)

    if t > 1:
        p = 1 / (1 + t)
        k = int(0.5 * np.log2(n * t**2 / (t - 1)))
        cum = binom.pmf(k=0, n=k, p=p)

    U = 0
    for i in range(1, r + 1):
        if t <= 1:
            U += math.pow(-t, i) * freq[i - 1]
        else:
            U += math.pow(-t, i) * (1 - cum) * freq[i - 1]
            if i <= k:
                cum += binom.pmf(k=i, n=k, p=p)
    return - U


def GT_PO(freq, n, t):
    """
    Smoothed Poisson Good Toulmin
    """
    r = len(freq)

    if t > 1:
        lam = 0.5 / t * np.log(n * math.pow(t + 1, 2) / (t - 1))
        cum = poisson.pmf(k=0, mu=lam)

    U = 0
    for i in range(1, r + 1):
        if t <= 1:
            U += math.pow(-t, i) * freq[i - 1]
        else:
            U += math.pow(-t, i) * (1 - cum) * freq[i - 1]
            cum += poisson.pmf(k=i, mu=lam)
    return - U


import matplotlib.pyplot as plt


freq = [118, 74, 44, 24, 29, 22, 20, 19, 20, 15, 12, 14, 6, 1, 6]
n = 10
gt_bi = [GT_BI(freq,n,t) for t in np.linspace(0.01,1.2,100)]
gt_po = [GT_PO(freq,n,t) for t in np.linspace(0.01,1.2,100)]

plt.plot(np.linspace(0.01,1.2,100)* sum(freq), gt_bi, label='bi 10')
plt.plot(np.linspace(0.01,1.2,100)* sum(freq), gt_po, label='po 10')


n = 100
gt_bi = [GT_BI(freq,n,t) for t in np.linspace(0.01,1.2,100)]
gt_po = [GT_PO(freq,n,t) for t in np.linspace(0.01,1.2,100)]

plt.plot(np.linspace(0.01,1.2,100) * sum(freq), gt_bi, label='bi 100')
plt.plot(np.linspace(0.01,1.2,100) * sum(freq), gt_po, label='po 100')

plt.legend()

plt.savefig('GT_pycode.pdf')


# parser = argparse.ArgumentParser(description='Calculate species richness estimates.')
# parser.add_argument('-d', '--data', metavar='D', type=str, nargs=1, 
#             help='input data with frequency counts (value in row k corresponds to the number of clones with exactly k individuals in the sample)')
# parser.add_argument('-o', '--output', metavar='O', type=str, nargs=1, 
#             help='output file with estimate of number of unseen species')
# parser.add_argument('-p', '--population', metavar='P', type=float, nargs=1, 
#             help='Population size')
# parser.add_argument('-s', '--seen', metavar='S', type=int, nargs=1, 
#             help='Percentage of seen population.')

# args = parser.parse_args()
# input = args.data[0]
# output = args.output[0]
# p = int(args.population[0])
# s = args.seen[0]

# with open(input, "r") as file:
#     freq = np.array([int(x) for x in file.readlines()])

# n = p * s / 100
# m = p - n
# t = m / n

# try:
#     UEF = GT_BI(freq, n, t)
#     if UEF < 0:
#         UEF = np.nan
# except:
#     UEF = np.nan

# try:
#     ULP = GT_PO(freq, n, t)
#     if ULP < 0:
#         ULP = np.nan
# except:
#     ULP = np.nan

# with open(output, "w") as file:
#     file.write("Efron-Thisted Good Toulmin\t{est:.3f}\n".format(est=UEF))
#     file.write("Smoothed Poisson Good Toulmin\t{est:.3f}".format(est=ULP))
