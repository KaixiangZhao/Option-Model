import math
import numpy
import scipy.stats
import scipy.optimize
from input import WARRANT_LIST, R, input_all, print_warrant_price

TOTAL_BS = []
TOTAL_NW = []
TOTAL_BSDA = []
TOTAL_UKHOV = []

f_bs_all = open('reg_data/bs_all.txt', 'w')
f_nw_all = open('reg_data/nw_all.txt', 'w')
f_bsda_all = open('reg_data/bsda_all.txt', 'w')
f_ukhov_all = open('reg_data/ukhov_all.txt', 'w')

f_bs_in = open('reg_data/bs_in.txt', 'w')
f_nw_in = open('reg_data/nw_in.txt', 'w')
f_bsda_in = open('reg_data/bsda_in.txt', 'w')
f_ukhov_in = open('reg_data/ukhov_in.txt', 'w')
f_bs_out = open('reg_data/bs_out.txt', 'w')
f_nw_out = open('reg_data/nw_out.txt', 'w')
f_bsda_out = open('reg_data/bsda_out.txt', 'w')
f_ukhov_out = open('reg_data/ukhov_out.txt', 'w')
f_bs_at = open('reg_data/bs_at.txt', 'w')
f_nw_at = open('reg_data/nw_at.txt', 'w')
f_bsda_at = open('reg_data/bsda_at.txt', 'w')
f_ukhov_at = open('reg_data/ukhov_at.txt', 'w')


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

def get_sigma(date, start_node, stock_match):
    temp = []
    for j in range(date - start_node, date):
        temp.append(numpy.log(stock_match[j]["price"] / stock_match[j - 1]["price"]))
    return numpy.var(temp)**0.5

def get_sigma_impldvol(date, stock_match, warrant):
    for everyday_warrant in warrant.everyday_price:
        if (everyday_warrant['Date'] - stock_match[date]['Date']).days == -1:
            return everyday_warrant['impldvol']
    return 0

# def classify(warrant_list):
#     classify_list = []
#     for warrant in warrant_list:
#         temp = [[[] for j in range(3)] for i in range(6)]
#         price_temp = [[[] for j in range(3)] for i in range(6)]
#         avg = [[0 for j in range(3)] for i in range(6)]
#         for everyday_warrant in warrant.everyday_price:
#             for stock in warrant.target_stock.everyday_price:
#                 if (everyday_warrant['Date'] - stock['Date']).days == 0:
#                     ratio = (stock['price'] - everyday_warrant['price']) / warrant.price
#
#             if ratio < 0.94 and (warrant.end_date - everyday_warrant['Date']).days < 60:
#                 temp[0][0].append(everyday_warrant['Date'])
#                 price_temp[0][0].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio < 0.94 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
#                     (warrant.end_date - everyday_warrant['Date']).days <= 180:
#                 temp[0][1].append(everyday_warrant['Date'])
#                 price_temp[0][1].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio < 0.94 and (warrant.end_date - everyday_warrant['Date']).days > 180:
#                 temp[0][2].append(everyday_warrant['Date'])
#                 price_temp[0][2].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 0.94 and ratio < 0.97 and (warrant.end_date - everyday_warrant['Date']).days < 60:
#                 temp[1][0].append(everyday_warrant['Date'])
#                 price_temp[1][0].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 0.94 and ratio < 0.97 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
#                     (warrant.end_date - everyday_warrant['Date']).days <= 180:
#                 temp[1][1].append(everyday_warrant['Date'])
#                 price_temp[1][1].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 0.94 and ratio < 0.97 and (warrant.end_date - everyday_warrant['Date']).days > 180:
#                 temp[1][2].append(everyday_warrant['Date'])
#                 price_temp[1][2].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 0.97 and ratio < 1.00 and (warrant.end_date - everyday_warrant['Date']).days < 60:
#                 temp[2][0].append(everyday_warrant['Date'])
#                 price_temp[2][0].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 0.97 and ratio < 1.00 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
#                     (warrant.end_date - everyday_warrant['Date']).days <= 180:
#                 temp[2][1].append(everyday_warrant['Date'])
#                 price_temp[2][1].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 0.97 and ratio < 1.00 and (warrant.end_date - everyday_warrant['Date']).days > 180:
#                 temp[2][2].append(everyday_warrant['Date'])
#                 price_temp[2][2].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 1.00 and ratio < 1.03 and (warrant.end_date - everyday_warrant['Date']).days < 60:
#                 temp[3][0].append(everyday_warrant['Date'])
#                 price_temp[3][0].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 1.00 and ratio < 1.03 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
#                     (warrant.end_date - everyday_warrant['Date']).days <= 180:
#                 temp[3][1].append(everyday_warrant['Date'])
#                 price_temp[3][1].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 1.00 and ratio < 1.03 and (warrant.end_date - everyday_warrant['Date']).days > 180:
#                 temp[3][2].append(everyday_warrant['Date'])
#                 price_temp[3][2].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 1.03 and ratio < 1.06 and (warrant.end_date - everyday_warrant['Date']).days < 60:
#                 temp[4][0].append(everyday_warrant['Date'])
#                 price_temp[4][0].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 1.03 and ratio < 1.06 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
#                     (warrant.end_date - everyday_warrant['Date']).days <= 180:
#                 temp[4][1].append(everyday_warrant['Date'])
#                 price_temp[4][1].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 1.03 and ratio < 1.06 and (warrant.end_date - everyday_warrant['Date']).days > 180:
#                 temp[4][2].append(everyday_warrant['Date'])
#                 price_temp[4][2].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 1.06 and (warrant.end_date - everyday_warrant['Date']).days < 60:
#                 temp[5][0].append(everyday_warrant['Date'])
#                 price_temp[5][0].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 1.06 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
#                     (warrant.end_date - everyday_warrant['Date']).days <= 180:
#                 temp[5][1].append(everyday_warrant['Date'])
#                 price_temp[5][1].append(everyday_warrant.get('sigma_bs', 0.5))
#             elif ratio >= 1.06 and (warrant.end_date - everyday_warrant['Date']).days > 180:
#                 temp[5][2].append(everyday_warrant['Date'])
#                 price_temp[5][2].append(everyday_warrant.get('sigma_bs', 0.5))
#
#         for i in range(0, 6):
#             for j in range(0, 3):
#                 avg[i][j] = numpy.mean(price_temp[i][j])
#         tmp = {}
#         tmp[warrant.code] = avg
#         classify_list.append(tmp)
#
#     return classify_list

def classify(warrant_list):
    for warrant in warrant_list:
        for everyday_warrant in warrant.everyday_price:
            for stock in warrant.target_stock.everyday_price:
                if (everyday_warrant['Date'] - stock['Date']).days == 0:
                    ratio = (stock['price'] - everyday_warrant['price']) / warrant.price
                    break
            if ratio < 0.94 and (warrant.end_date - everyday_warrant['Date']).days < 60:
                everyday_warrant['class'] = 0
            elif ratio < 0.94 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
                    (warrant.end_date - everyday_warrant['Date']).days <= 180:
                everyday_warrant['class'] = 1
            elif ratio < 0.94 and (warrant.end_date - everyday_warrant['Date']).days > 180:
                everyday_warrant['class'] = 2
            elif ratio >= 0.94 and ratio < 0.97 and (warrant.end_date - everyday_warrant['Date']).days < 60:
                everyday_warrant['class'] = 3
            elif ratio >= 0.94 and ratio < 0.97 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
                    (warrant.end_date - everyday_warrant['Date']).days <= 180:
                everyday_warrant['class'] = 4
            elif ratio >= 0.94 and ratio < 0.97 and (warrant.end_date - everyday_warrant['Date']).days > 180:
                everyday_warrant['class'] = 5
            elif ratio >= 0.97 and ratio < 1.00 and (warrant.end_date - everyday_warrant['Date']).days < 60:
                everyday_warrant['class'] = 6
            elif ratio >= 0.97 and ratio < 1.00 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
                    (warrant.end_date - everyday_warrant['Date']).days <= 180:
                everyday_warrant['class'] = 7
            elif ratio >= 0.97 and ratio < 1.00 and (warrant.end_date - everyday_warrant['Date']).days > 180:
                everyday_warrant['class'] = 8
            elif ratio >= 1.00 and ratio < 1.03 and (warrant.end_date - everyday_warrant['Date']).days < 60:
                everyday_warrant['class'] = 9
            elif ratio >= 1.00 and ratio < 1.03 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
                    (warrant.end_date - everyday_warrant['Date']).days <= 180:
                everyday_warrant['class'] = 10
            elif ratio >= 1.00 and ratio < 1.03 and (warrant.end_date - everyday_warrant['Date']).days > 180:
                everyday_warrant['class'] = 11
            elif ratio >= 1.03 and ratio < 1.06 and (warrant.end_date - everyday_warrant['Date']).days < 60:
                everyday_warrant['class'] = 12
            elif ratio >= 1.03 and ratio < 1.06 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
                    (warrant.end_date - everyday_warrant['Date']).days <= 180:
                everyday_warrant['class'] = 13
            elif ratio >= 1.03 and ratio < 1.06 and (warrant.end_date - everyday_warrant['Date']).days > 180:
                everyday_warrant['class'] = 14
            elif ratio >= 1.06 and (warrant.end_date - everyday_warrant['Date']).days < 60:
                everyday_warrant['class'] = 15
            elif ratio >= 1.06 and (warrant.end_date - everyday_warrant['Date']).days >=60 and \
                    (warrant.end_date - everyday_warrant['Date']).days <= 180:
                everyday_warrant['class'] = 16
            elif ratio >= 1.06 and (warrant.end_date - everyday_warrant['Date']).days > 180:
                everyday_warrant['class'] = 17

def classify_result():
    result = [[0 for j in range(3)] for i in range(6)]
    for warrant in WARRANT_LIST:
        for everyday_warrant in warrant.everyday_price:
            i = int(everyday_warrant['class'] / 3)
            j = everyday_warrant['class'] % 3
            result[i][j] += 1
    return result

def convergence():
    convergence_list = []
    for warrant in WARRANT_LIST:
        stock_match = []
        for day in warrant.target_stock.everyday_price:
            if (day["Date"] - warrant.end_date).days <= 0 and \
               (warrant.start_date - day["Date"]).days <= 0:
                stock_match.append(day)

        if (stock_match[0]['Date'] - warrant.end_date).days == 0:
            temp = {}
            temp['code'] = warrant.code
            temp['status'] = 0
            temp['error'] = abs(warrant.price + warrant.everyday_price[-1]['price'] -
                                stock_match[0]['price']) / stock_match[0]['price']
            convergence_list.append(temp)
        elif (stock_match[0]['Date'] - warrant.end_date).days > 0:
            temp = {}
            temp['code'] = warrant.code
            temp['status'] = 1
            temp['error'] = abs(warrant.price + warrant.everyday_price[-1]['price'] -
                                stock_match[0]['price']) / stock_match[0]['price']
            convergence_list.append(temp)
        else:
            temp = {}
            temp['code'] = warrant.code
            temp['status'] = 2
            temp['error'] = abs(warrant.price + warrant.everyday_price[-1]['price'] -
                                stock_match[0]['price']) / stock_match[0]['price']
            convergence_list.append(temp)
    return convergence_list

def calculate_sigma(warrant_list):
    for warrant in warrant_list:
        stock_match = []
        for day in warrant.target_stock.everyday_price:
            if (day["Date"] - warrant.end_date).days <= 0 and \
               (warrant.start_date - day["Date"]).days <= 0:
                stock_match.append(day)

        stock_sort(stock_match)

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

    f_bs = open('reg_data/bs_{}.txt'.format(warrant.code), 'w')
    f_nw = open('reg_data/nw_{}.txt'.format(warrant.code), 'w')
    f_bsda = open('reg_data/bsda_{}.txt'.format(warrant.code), 'w')
    f_ukhov = open('reg_data/ukhov_{}.txt'.format(warrant.code), 'w')

    # get the data that match the date
    stock_match = []
    for day in warrant.target_stock.everyday_price:
        if (day["Date"] - warrant.end_date).days <= 0 and \
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
        sigma = get_sigma(i, start_node, stock_match)
        #sigma = get_sigma_impldvol(i, stock_match, warrant)

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
                s_k = stock_match[i]['price'] / (warrant.price * math.exp(-R * t))
                avg_bs.append(abs(tempp["price"] - warrant_r["price"]) / warrant_r["price"])
                f_bs.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                f_bs_all.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                if s_k >= 1.1:
                    f_bs_in.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                elif s_k < 0.9:
                    f_bs_out.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                else:
                    f_bs_at.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))

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
                s_k = stock_match[i]['price'] / (warrant.price * math.exp(-R * t))
                avg_nw.append(abs(tempp["price"] - warrant_r["price"]) / warrant_r["price"])
                f_nw.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                f_nw_all.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                if s_k >= 1.1:
                    f_nw_in.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                elif s_k < 0.9:
                    f_nw_out.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                else:
                    f_nw_at.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))

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
                s_k = stock_match[i]['price'] / (warrant.price * math.exp(-R * t))
                avg_bsda.append(abs(tempp["price"] - warrant_r["price"]) / warrant_r["price"])
                f_bsda.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                f_bsda_all.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                if s_k >= 1.1:
                    f_bsda_in.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                elif s_k < 0.9:
                    f_bsda_out.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                else:
                    f_bsda_at.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))

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
                s_k = stock_match[i]['price'] / (warrant.price * math.exp(-R * t))
                avg_ukhov.append(abs(tempp["price"] - warrant_r["price"]) / warrant_r["price"])
                f_ukhov.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                f_ukhov_all.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                if s_k >= 1.1:
                    f_ukhov_in.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                elif s_k < 0.9:
                    f_ukhov_out.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))
                else:
                    f_ukhov_at.write('{} {} {} {}\n'.format((tempp["price"] - warrant_r["price"]) / warrant_r["price"],
                           sigma,
                           s_k,
                           t * 1.0 / 365.0))

    TOTAL_BS.extend(avg_bs)
    TOTAL_NW.extend(avg_nw)
    TOTAL_BSDA.extend(avg_bsda)
    TOTAL_UKHOV.extend(avg_ukhov)


    mean_bs = numpy.mean(avg_bs)
    mean_nw = numpy.mean(avg_nw)
    mean_bsda = numpy.mean(avg_bsda)
    mean_ukhov = numpy.mean(avg_ukhov)
    std_bs = numpy.std(avg_bs)
    std_nw = numpy.std(avg_nw)
    std_bsda = numpy.std(avg_bsda)
    std_ukhov = numpy.std(avg_ukhov)

    print('{}& {}& {}& {}& {}& {}& {}& {}& {}'.format(warrant.code, len(warrant.everyday_price),
          mean_bs, std_bs, mean_nw, std_nw, mean_bsda, std_bsda, mean_ukhov, std_ukhov))

    f_bs.close()
    f_nw.close()
    f_bsda.close()
    f_ukhov.close()

    return (bs_result, nw_result, bsda_result, ukhov_result)

def get_result(code):
    for warrant in WARRANT_LIST:
        if warrant.code == code:
            return calculate(warrant)
    print("No such warrant!\n")

def get_all_result():
    for warrant in WARRANT_LIST:
        get_result(warrant.code)

    mean_bs = numpy.mean(TOTAL_BS)
    mean_nw = numpy.mean(TOTAL_NW)
    mean_bsda = numpy.mean(TOTAL_BSDA)
    mean_ukhov = numpy.mean(TOTAL_UKHOV)
    std_bs = numpy.std(TOTAL_BS)
    std_nw = numpy.std(TOTAL_NW)
    std_bsda = numpy.std(TOTAL_BSDA)
    std_ukhov = numpy.std(TOTAL_UKHOV)

    print("All", len(TOTAL_BS),
          mean_bs, std_bs, mean_nw, std_nw, mean_bsda, std_bsda, mean_ukhov, std_ukhov)


if __name__ == '__main__':
    input_all()
    #get_result('031001')
    get_all_result()
    #calculate_sigma(WARRANT_LIST)
    #classify(WARRANT_LIST)
    #print(classify_result())
    #print(convergence())
    f_bs_all.close()
    f_nw_all.close()
    f_bsda_all.close()
    f_ukhov_all.close()
    f_bs_in.close()
    f_nw_in.close()
    f_bsda_in.close()
    f_ukhov_in.close()
    f_bs_out.close()
    f_nw_out.close()
    f_bsda_out.close()
    f_ukhov_out.close()
    f_bs_at.close()
    f_nw_at.close()
    f_bsda_at.close()
    f_ukhov_at.close()


