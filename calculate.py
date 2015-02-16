import math
import numpy
import scipy
import statistics
from input import WARRANT_LIST, R

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

def ukhov(s, x, t, r, d, ns, nw, sigstar, gama):
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

def calculate(warrant):
    bs_result = []
    bsda_result = []
    nw_result = []
    ukhov_result = []

    avg_bs = []
    avg_bsda = []
    avg_nw = []
    avg_ukhov = []

    # get the data that match the date
    stock_match = []
    for day in warrant.target_stock.everyday_price:
        if (day["Date"] - warrant.end_date).days < 0 and \
           (warrant.start_date - day["Date"]).days <= 365:
            stock_match.append(day)

    # sort the date
    stock_sort(stock_match)

    # calculate the start time
    for i in range(len(stock_match)):
        if (stock_match[i]["Date"] - warrant.start_date).days >= 0:
            start_node = i
            break

    for i in range(start_node + 1, len(stock_match)):
        temp = []
        for j in range(i - start_node, i):
            temp.append(numpy.log(stock_match[j]["price"] / \
                                  stock_match[j - 1]["price"]))
            sigma = numpy.var(temp)**0.5

        t = (warrant.end_date - stock_match[i]["Date"]).days
        d_1 = (numpy.log(stock_match[i]["price"] / warrant.price) + \
               (R + 0.5 * sigma * sigma) * t) / (sigma * math.pow(t, 0.5))
        d_2 = d_1 - sigma * math.pow(t, 0.5)

        #calculate the BS model
        if warrant.species == "C":
            result = stock_match[i]["price"] * scipy.stats.norm.cdf(d_1) \
                     - warrant.price * math.exp(-R * t) * \
                     scipy.stats.norm.cdf(d_2)
            tempp = {}
            tempp["price"] = result
            tempp["Date"] = stock_match[i]["Date"]
            bs_result.append(tempp)

        for warrant_r in warrant.everyday_price:
            if (tempp["Date"] - warrant_r["Date"]).days == 0:
                avg_bs.append(abs(tempp["price"] - warrant_r["price"]) / \
                              warrant_r["price"])

        # calculate the NW model
        m = warrant.amount
        n = warrant.target_stock.amount
        y = warrant.portion

        if warrant.species == "C":
            result = (n / (n + m * y)) * result
            tempp = {}
            tempp["price"] = result
            tempp["Date"] = stock_match[i]["Date"]
            nw_result.append(tempp)

        for warrant_r in warrant.everyday_price:
            if (tempp["Date"] - warrant_r["Date"]).days == 0:
                avg_nw.append(abs(tempp["price"] - warrant_r["price"]) / \
                              warrant_r["price"])

        # calculate the BSDA model result
        if warrant.species == "C":
            result = bsda(stock_match[i]["price"], warrant.price,
                          t, R, 0, n, m, sigma, y)
            tempp = {}
            tempp["price"] = result
            tempp["Date"] = stock_match[i]["Date"]
            bsda_result.append(tempp)

        for warrant_r in warrant.everyday_price:
            if (tempp["Date"] - warrant_r["Date"]).days == 0:
                avg_bsda.append(abs(tempp["price"] - warrant_r["price"]) \
                                / warrant_r["price"])

        #calculate the Ukhov model result
        if warrant.species == "C":
            result = ukhov(stock_match[i]["price"], warrant.price,
                           t, R, 0, n, m, sigma, y)
            tempp = {}
            tempp["price"] = result
            tempp["Date"] = stock_match[i]["Date"]
            ukhov_result.append(tempp)

        for warrant_r in warrant.everyday_price:
            if (tempp["Date"] - warrant_r["Date"]).days == 0:
                avg_ukhov.append(abs(tempp["price"] - warrant_r["price"])\
                                 / warrant_r["price"])

    avg = statistics.mean(avg_bs)
    print("{} \t {} \t {}".format("BS", warrant.code, avg))

    avg = statistics.mean(avg_nw)
    print("{} \t {} \t {}".format("NW", warrant.code, avg))

    avg = statistics.mean(avg_bsda)
    print("{} \t {} \t {}".format("BSDA", warrant.code, avg))

    avg = statistics.mean(avg_ukhov)
    print("{} \t {} \t {}".format("Ukhov", warrant.code, avg))
