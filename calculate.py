import math
import numpy
import scipy
import input

def stock_cmp(stock1, stock2):
    return (stock1["Date"] - stock2["Date"]).days < 0

def stock_sort(stocks):
    temp = {}
    for i in range(len(stocks) - 1):
        mini = i
        for j in range(i + 1, len(stocks)):
            if False == stock_cmp(stocks[mini], stocks[j]):
                mini = j
            if mini != i:
                temp = stocks[mini]
                stocks[mini] = stocks[i]
                stocks[i] = temp

def sub_bsda(s, x, t, r, d, ns, nw, sig, gama, w):
    sart = s * math.exp(-d * t) + w * nw / ns
    d_1 = (numpy.log(sart / x) + (r + sig**2 / 2) * t) / (sig * t**0.5)
    d_2 = d_1 - sig * t**0.5
    result = (ns / (ns / gama + nw)) * \
               (sart * scipy.stats.norm.cdf(d_1) - \
                x * math.exp(-r * t) * scipy.stats.norm.cdf(d_2))
    return result

def bsda(s, x, t, r, d, ns, nw, sig, gama):
    w_old = 0
    y_old = w_old - sub_bsda(s, x, t, r, d, ns, nw, sig, gama, w_old)

    w = 1
    y = w - sub_bsda(s, x, t, r, d, ns, nw, sig, gama, w)

    while abs(y) / w >= 1e-8:
        w_temp = w
        w = w - y * (w - w_old) / (y - y_old)
        y_old = y
        y = w - sub_bsda(s, x, t, r, d, ns, nw, sig, gama, w)
        w_old = w_temp

    return w

def sub_ukhov(s, x, t, r, d, ns, nw, sigstar, gama, v, sig):
    sart = s * math.exp(-d * t)
    d_1 = (numpy.log(gama * v / (ns * x)) + (r + sig**2 / 2) * t) / \
          (sig * t**0.5)
    d_2 = d_1 - sig * t**0.5
    delta_s = (ns + gama * nw - nw * gama * scipy.stats.norm.cdf(d_1)) / \
              (ns * (ns + gama * nw))
    v = s * ns + nw * (1/ (ns + gama * nw)) * \
        (gama * v * scipy.stats.norm.cdf(d_1) - \
         math.exp(-r * t) * ns * x * scipy.stats.norm.cdf(d_2))
    sig = sigstar / (v * delta_s / s)
    temp = {}
    temp["V"] = v
    temp["sig"] = sig
    return temp

def Ukhov(s, x, t, r, d, ns, nw, sigstar, gama):
    v_old = 0
    sig_old = 0
    temp = sub_ukhov(s, x, t, r, d, ns, nw, sigstar, gama, v_old, sig_old)
    y_old = v_old - temp["V"]

    v = 1
    sig = 1
    temp = sub_ukhov(s, x, t, r, d, ns, nw, sigstar, gama, v, sig)
    y = v - temp["V"]

    while abs(y) / v >= 1e-8:
        v_temp = v
        sig_temp = sig
        v = v - y * (v - v_old) / (y - y_old)
        sig = sig - y * (sig - sig_old) / (y - y_old)
        y_old = y
        temp = sub_ukhov(s, x, t, r, d, ns, nw, sigstar, gama, v, sig)
        y = v - temp["V"]
        v_old = v_temp
        sig_old = sig_temp

    w = (v - s * ns) / nw
    return w
