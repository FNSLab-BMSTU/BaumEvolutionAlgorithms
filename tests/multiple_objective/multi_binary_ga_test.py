from baumeva import FFGA, VEGA
from baumeva.ga import DynamicPenalty
import math
import time


class timex:
    def __enter__(self):
        # Фиксация времени старта процесса
        self.t = time.time()
        return self

    def __exit__(self, type, value, traceback):
        # Вывод времени работы
        print('Processing time: {:.2f} sec.'.format(time.time() - self.t))


def griewank_func(value_list):
    res_sum = 0
    res_mpl = 1
    for i, value in enumerate(value_list):
        res_sum += (value*value) / 4000
        res_mpl *= math.cos(value/math.sqrt(i+1))
    return res_sum - res_mpl + 1


def sum_of_squares(value_list):
    result = 0
    for value in value_list:
        result += value*value
    return result

# ---------------------------- https://github.com/P-N-Suganthan/2021-RW-MOP/blob/main/RWMOP.pdf -----------------------
# 2.1.3 Two Bar Truss Design (page 5)

# VEGA - don't work: num_generations=100 - 2.62s, num_generations=1000 - 27.18s
# FFGA - don't work: num_generations=10 - 3.85s, num_generations=20 - 25.65s
def rcm03(x: list) -> tuple:
    """
    conditions=['optimize']*2 + ['<=']*3
    gens=((0.00001, 100, 0.01), (0.00001, 100, 0.01), (1, 3, 0.01))
    """
    f1 = x[0] * math.sqrt(16 + x[2] ** 2) + x[1] * math.sqrt(1 + x[2] ** 2)
    f2 = 20 * math.sqrt(16 + x[2] ** 2) / (x[2] * x[0])
    g1 = f1 - 0.1
    g2 = f2 - 10000
    g3 = 80 * math.sqrt(1 + x[2] ** 2) / (x[2] * x[1]) - 10000
    return f1, f2, g1, g2, g3


# --------------------------- https://pymoo.org/problems/ ------------------------------------------------------------

# VEGA - don't work: num_generations=100 - 4.40s, num_generations=1000 - 40.27s
# FFGA - acceptable work: num_generations=100 - 4.69s, num_generations=1000 - 48.79s
def zdt1(x: list) -> tuple:
    """
    conditions=['optimize']*2
    gens=((0, 1, 0.01),)*30
    """
    f1 = x[0]
    g = 1 + 9 / (len(x) - 1) * sum(x[1:])
    h = 1 - math.sqrt(f1 / g)
    f2 = g * h

    return f1, f2


# VEGA - don't work: num_generations=100 - 3.90s, num_generations=1000 - 39.75s
# FFGA - acceptable work: num_generations=100 - 4.10s, num_generations=1000 - 42.52s
def zdt2(x: list) -> tuple:
    """
    conditions=['optimize']*2
    gens=((0, 1, 0.01),)*30
    """
    f1 = x[0]
    g = 1 + 9 / (len(x) - 1) * sum(x[1:])
    h = 1 - (f1 / g) ** 2
    f2 = g * h

    return f1, f2


# VEGA - don't work: num_generations=100 - 4.32s, num_generations=1000 - 42.27s
# FFGA - acceptable work: num_generations=100 - 4.71s, num_generations=1000 - 50.00s
def zdt3(x: list) -> tuple:
    """
    conditions=['optimize']*2
    gens=((0, 1, 0.01),)*30
    """
    f1 = x[0]
    g = 1 + 9 / (len(x) - 1) * sum(x[1:])
    h = 1 - math.sqrt(f1 / g) - (f1 / g) * math.sin(10 * math.pi * f1)
    f2 = g * h

    return f1, f2


# VEGA - don't work: num_generations=100 - 5.63s, num_generations=1000 - 55.81s
# FFGA - don't work: num_generations=100 - 5.55s, num_generations=1000 - 56.13s
def zdt4(x: list) -> tuple:
    """
    conditions=['optimize']*2
    gens=tuple([(0, 1, 0.01)]+[(-10, 10, 0.01)]*10)
    """
    f1 = x[0]

    s = 0
    for xi in x[1:]:
        s += xi ** 2 - 10 * math.cos(4 * math.pi * xi)

    g = 1 + 10 * (len(x) - 1) + s
    h = 1 - math.sqrt(f1 / g)
    f2 = g * h

    return f1, f2


# VEGA - acceptable work: num_generations=100 - 5.19s, num_generations=1000 - 52.12s
# FFGA - acceptable work: num_generations=100 - 8.00s, num_generations=200 - 126.20s
def zdt5(x: list) -> tuple:
    """
    conditions=['optimize']*2
    gens=tuple([(0, 2**30-1, 1)]+[(0, 2**5-1, 1)]*10)
    """
    u = lambda x: bin(int(x)).count('1')

    f1 = 1 + u(x[0])

    g = 0
    for xi in x[1:]:
        ui = u(xi)
        vi = ui + 2 if ui < 5 else 1
        g += vi
    h = 1 / f1
    f2 = g * h

    return f1, f2


# VEGA - don't work: num_generations=100 - 2.96s, num_generations=1000 - 26.15s
# FFGA - acceptable work: num_generations=100 - 2.99s, num_generations=300 - 68.68s
def zdt6(x: list) -> tuple:
    """
    conditions=['optimize']*2
    gens=((0, 1, 0.01),)*10
    """
    f1 = 1 - math.exp(-4 * x[0]) * (math.sin(6 * math.pi * x[0])) ** 6

    g = 1 + 9 * (sum(x[1:]) / 9) ** 0.25
    h = 1 - (f1 / g) ** 2
    f2 = g * h

    return f1, f2


# VEGA - acceptable work: num_generations=100 - 2.49s, num_generations=1000 - 25.59s
# FFGA - acceptable work: num_generations=10 - 7.12s, num_generations=20 - 47.24s
# (контур фронта Парето похож, но вместо линии расходящийся шлейф)
def bnh(x: list) -> tuple:
    """
    conditions=['optimize']*2 + ['<=']*2
    gens=tuple([(0, 5, 0.01)] + [(0, 3, 0.01)])
    """
    f1 = 4 * x[0] ** 2 + 4 * x[1] ** 2
    f2 = (x[0] - 5) ** 2 + (x[1] - 5) ** 2

    c1 = (x[0] - 5) ** 2 + x[1] ** 2 - 25
    c2 = 7.7 - (x[0] - 8) ** 2 - (x[1] + 3) ** 2

    return f1, f2, c1, c2


# VEGA - don't work: num_generations=100 - 4.89s, num_generations=1000 - 48.61s
# FFGA - don't work: num_generations=10 - 12.27s, num_generations=20 - 90.05s
def osy(x: list) -> tuple:
    """
    conditions=['optimize']*2 + ['<=']*6
    gens=tuple([(0, 10, 0.01)]*2 + [(1, 5, 0.01)] + [(0, 6, 0.01)] + [(1, 5, 0.01)] + [(0, 10, 0.01)])
    """
    f1 = -(25 * (x[0] - 2) ** 2 + (x[1] - 2) ** 2 + (x[2] - 1) ** 2 + (x[3] - 4) ** 2 + (x[4] - 1) ** 2)
    f2 = x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2 + x[4] ** 2 + x[5] ** 2

    c = [0] * 6
    c[0] = (x[0] + x[1] - 2) / 2
    c[1] = (6 - x[0] - x[1]) / 6
    c[2] = (2 - x[1] + x[0]) / 2
    c[3] = (2 - x[0] + 3 * x[1]) / 2
    c[4] = (4 - (x[2] - 3) ** 2 - x[3]) / 4
    c[5] = ((x[4] - 3) ** 2 + x[5] - 4) / 4
    c = [-i for i in c]

    return tuple([f1] + [f2] + c)


# VEGA - don't work: num_generations=100 - 2.03s, num_generations=1000 - 19.77s
# FFGA - don't work: num_generations=10 - 5.87s, num_generations=20 - 37.76s
def tnk(x: list) -> tuple:
    """
    conditions=['optimize']*2 + ['<=']*2
    gens=((0, math.pi, 0.01),)*2
    """
    f1 = x[0]
    f2 = x[1]

    c1 = x[0] ** 2 + x[1] ** 2 - 1 - 0.1 * math.cos(16 * math.atan2(x[0], x[1]))
    c2 = (x[0] - 0.5) ** 2 + (x[1] - 0.5) ** 2 - 0.5

    return f1, f2, -c1, c2


# --------------------- https://desdeo-problem.readthedocs.io/en/latest/problems.html ----------------------------------

def carside(x: list) -> tuple:
    """
    conditions=['optimize']*3 + ['<=']*10
    gens=tuple([(0.5, 1.5, 0.01)] + [(0.45, 1.35, 0.01)] + [(0.5, 1.5, 0.01)]*2 + [(0.875, 2.625, 0.01)] + [(0.4, 1.2, 0.01)]*2)
    """
    f1 = 1.98 + 4.9 * x[0] + 6.67 * x[1] + 6.98 * x[2] + 4.01 * x[3] + 1.78 * x[4] + 0.00001 * x[5] + 2.73 * x[6]
    f2 = 4.72 - 0.5 * x[3] - 0.19 * x[1] * x[2]

    v_mbp = 10.58 - 0.674 * x[0] * x[1] - 0.67275 * x[1]
    v_fd = 16.45 - 0.489 * x[2] * x[6] - 0.843 * x[4] * x[5]
    f3 = 0.5 * v_mbp * v_fd

    g = [0] * 10
    g[0] = 1 - 1.16 + 0.3717 * x[1] * x[3] + 0.0092928 * x[2]
    g[1] = 0.32 - 0.261 + 0.0159 * x[0] * x[1] + 0.06486 * x[0] + 0.019 * x[1] * x[6] - 0.0144 * x[2] * x[
        4] - 0.0154464 * x[5]
    g[2] = 0.32 - 0.214 - 0.00817 * x[4] + 0.045195 * x[0] + 0.0135168 * x[0] - 0.03099 * x[1] * x[5] + 0.018 * x[1] * \
           x[6] - \
           0.007176 * x[2] - 0.023232 * x[2] + 0.00364 * x[4] * x[5] + 0.018 * x[1] ** 2
    g[3] = 0.32 - 0.74 + 0.61 * x[1] + 0.031296 * x[2] + 0.031872 * x[6] - 0.227 * x[1] ** 2
    g[4] = 32 - 28.98 - 3.818 * x[2] + 4.2 * x[0] * x[1] - 1.27296 * x[5] + 2.68065 * x[6]
    g[5] = 32 - 33.86 - 2.95 * x[2] + 5.057 * x[0] * x[1] + 3.795 * x[1] + 3.4431 * x[6] - 1.45728
    g[6] = 32 - 46.36 + 9.9 * x[1] + 4.4505 * x[0]
    g[7] = 4 - f2
    g[8] = 9.9 - v_mbp
    g[9] = 15.7 - v_fd
    g = [-i for i in g]

    return tuple([f1] + [f2] + [f3] + g)


# VEGA - don't work: num_generations=100 - 5.59s, num_generations=1000 - 55.35s
# FFGA - don't work: num_generations=10 - 10.02s, num_generations=20 - 67.79s
def welded_beam(x: list) -> tuple:
    """
    conditions=['optimize']*2 + ['<=']*4
    gens=tuple([(0.125, 5, 0.01)] + [(0.1, 10, 0.01)]*2 + [(0.125, 5, 0.01)])
    """
    P = 6000
    E = 30 * 10 ** 6
    tau_max = 13600
    L = 14
    G = 12 * 10 ** 6
    sigma_max = 30000

    M = P * (L + x[1] / 2)
    R = math.sqrt(x[1] ** 2 / 4 + ((x[0] + x[2]) / 2) ** 2)
    J = 2 * (math.sqrt(2) * x[0] * x[1] * (x[1] ** 2 / 12 + ((x[0] + x[2]) / 2) ** 2))
    tau_1 = P / (math.sqrt(2) * x[0] * x[1])
    tau_2 = M * R / J
    tau = math.sqrt(tau_1 ** 2 + 2 * tau_1 * tau_2 * x[1] / (2 * R) + tau_2 ** 2)
    sigma = 6 * P * L / (x[3] * x[2] ** 2)
    Pc = 4.013 * E * math.sqrt(x[2] ** 2 * x[3] ** 6 / 36) / (L ** 2) * (1 - x[2] / (2 * L) * math.sqrt(E / (4 * G)))

    f1 = 1.10471 * x[0] ** 2 * x[1] + 0.04811 * x[2] * x[3] * (14 + x[1])
    f2 = 4 * P * L ** 3 / (E * x[3] * x[2] ** 3)

    g = [0] * 4
    g[0] = tau_max - tau
    g[1] = sigma_max - sigma
    g[2] = x[3] - x[1]
    g[3] = Pc - P
    g = [-i for i in g]

    return tuple([f1] + [f2] + g)


# ------------------------------- https://arxiv.org/pdf/1612.07603v2 ---------------------------------------------------
# DAS-CMOP1-9 (page 8)

# VEGA - don't work: num_generations=100 - 24.83s, num_generations=1000 - 247.32s
# FFGA - don't work: num_generations=10 - 74.30s, num_generations=20 - 525.21s
def das_cmop1(x: list) -> tuple:
    """
    conditions=['optimize']*2 + ['<=']*12
    gens=((0, 1, 0.01),)*30
    """
    eta = 0.25
    zeta = 0.25
    gamma = 0.25

    g1 = 0
    g2 = 0
    for i in range(1, len(x)):
        if i % 2 != 0:
            g1 += (x[i] - math.sin(0.5 * math.pi * x[0])) ** 2
        else:
            g2 += (x[i] - math.cos(0.5 * math.pi * x[0])) ** 2

    f1 = x[0] + g1
    f2 = 1 - x[0] ** 2 + g2

    a = 20
    d = 0.5
    b = 2 * eta - 1
    e = d - math.log(zeta)
    r = gamma / 2
    ak = 0.3
    bk = 1.2
    theta = -0.25 * math.pi
    pk = [0, 1, 0, 1, 2, 0, 1, 2, 3]
    qk = [1.5, 0.5, 2.5, 1.5, 0.5, 3.5, 2.5, 1.5, 0.5]

    c = [0] * 12
    c[0] = math.sin(a * math.pi * x[0]) - b
    c[1] = (e - g1) * (g1 - d)
    c[2] = (e - g2) * (g2 - d)
    for k in range(len(pk)):
        c[k + 3] = ((f1 - pk[k]) * math.cos(theta) - (f2 - qk[k]) * math.sin(theta)) ** 2 / ak + \
                   ((f1 - pk[k]) * math.sin(theta) + (f2 - qk[k]) * math.cos(theta)) ** 2 / bk - r
    c = [-i for i in c]

    return tuple([f1] + [f2] + c)


# VEGA - don't work: num_generations=100 - 19.17s, num_generations=1000 - 187.61s
# FFGA - don't work: num_generations=10 - 51.92s, num_generations=20 - 381.10s
def das_cmop7(x: list) -> tuple:
    """
    conditions=['optimize']*3 + ['<=']*7
    gens=((0, 1, 0.01),)*30
    """
    eta = 0.75
    zeta = 0.75
    gamma = 0.75

    s = 0
    for xi in x[2:]:
        s += (xi - 0.5) ** 2 - math.cos(20 * math.pi * (xi - 0.5))
    g = (len(x) - 2) + s

    f = [0] * 3
    f[0] = x[0] * x[1] + g
    f[1] = x[1] * (1 - x[0]) + g
    f[2] = 1 - x[1] + g

    a = 20
    d = 0.5
    b = 2 * eta - 1
    e = d - math.log(zeta)
    r = gamma / 2

    c = [0] * 7
    c[0] = math.sin(a * math.pi * x[0]) - b
    c[1] = math.cos(a * math.pi * x[1]) - b
    c[2] = (e - g) * (g - d)
    for k in range(3):
        s = 0
        for j in range(3):
            if j != k:
                s += f[j] ** 2
        c[k + 3] = s + (f[k] - 1) ** 2 - r ** 2
    s = 0
    for j in range(3):
        s += (f[j] - 1 / math.sqrt(3)) ** 2
    c[6] = s - r ** 2
    c = [-i for i in c]

    return tuple(f + c)


# https://staff.fmi.uvt.ro/~daniela.zaharie/ma2019/Projects/ResearchPapers/MultiobjectiveOptimization/MOEA_TestSuites_2019.pdf
# MW1-14 (page 6)

# VEGA - don't work: num_generations=100 - 8.24s, num_generations=1000 - 85.50s
# FFGA - don't work: num_generations=10 - 9.52s, num_generations=20 - 46.52s
def mw3(x: list) -> tuple:
    """
    conditions=['optimize']*2 + ['<=']*2
    gens=((0, 1, 0.01),)*15
    """
    s = 0
    for i in range(1, len(x)):
        s += 2 * (x[i] + (x[i - 1] - 0.5) ** 2 - 1) ** 2
    g = 1 + s

    f1 = x[0]
    f2 = g * (1 - f1 / g)

    l = math.sqrt(2) * (f2 - f1)
    c1 = 1.05 - f1 - f2 + 0.45 * math.sin((0.75 * math.pi * l) ** 6)
    c2 = 0.85 - f1 - f2 + 0.3 * math.sin((0.75 * math.pi * l) ** 2)

    return f1, f2, -c1, c2


# ---------- https://ls11-www.cs.tu-dortmund.de/people/rudolph/publications/papers/emo2007.pdf -------------------------
# SYM-PART (page 3)

# VEGA - acceptable work: num_generations=100 - 1.45s, num_generations=1000 - 15.22s
# FFGA - acceptable work: num_generations=10 - 0.46s, num_generations=40 - 48.19s
def sym_part(x: list) -> tuple:
    """
    conditions=['optimize']*2
    gens=((0, 1, 0.01),)*2
    """
    f1 = (x[0] + 1) ** 2 + x[1] ** 2
    f2 = (x[0] - 1) ** 2 + x[1] ** 2

    return f1, f2


# ------------------------------------- some simple quadrics -----------------------------------------------------------

def quadratic(x: list) -> tuple:
    """
    conditions=['optimize']*3
    gens=tuple([(0, 100, 0.01)]*2)
    """
    return (x[0] - 2)**2 + (x[1] - 2)**2, x[0]*x[0] + x[1]*x[1], (x[0] - 4)**2 + (x[1] - 1)**2


# moga = FFGA(num_generations=10,
#             num_individ=100,
#             gens=((0.00001, 100, 0.01), (0.00001, 100, 0.01), (1, 3, 0.01)),
#             obj_function=rcm03,
#             # obj_value=0,
#             conditions=['optimize']*2 + ['<=']*3,
#             penalty=DynamicPenalty(),
#             is_gray=False,
#             mutation_lvl=0.05,
#             early_stop=None)

moga = FFGA(num_generations=1000,
            num_individ=20,
            conditions=['optimize'] * 2 + ['<='] * 6,
            gens=tuple([(0, 10, 0.01)] * 2 + [(1, 5, 0.01)] + [(0, 6, 0.01)] + [(1, 5, 0.01)] + [(0, 10, 0.01)]),
            obj_function=osy,
            is_gray=False,
            mutation_lvl=0.05,
            early_stop=None)

with timex():
    ga_data = moga.optimize()
# ga_data_conditions = bin_ga_conditions.optimize()
