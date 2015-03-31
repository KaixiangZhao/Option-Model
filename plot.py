import matplotlib.pyplot as plt
import numpy
import math
import scipy.stats
from calculate import WARRANT_LIST, calculate, input_all, stock_sort, R

def plot(warrant, bs_result, nw_result, bsda_result, ukhov_result):
    bs_plot = []
    nw_plot = []
    bsda_plot = []
    ukhov_plot = []
    real_plot = []

    bs_lag = []
    nw_lag = []
    bsda_lag = []
    ukhov_lag = []

    bs_avg = []
    nw_avg = []
    bsda_avg = []
    ukhov_avg = []
    real_avg = []

    for i in bs_result:
        bs_plot.append(i["price"])

    for i in nw_result:
        nw_plot.append(i["price"])

    for i in bsda_result:
        bsda_plot.append(i["price"])

    for i in ukhov_result:
        ukhov_plot.append(i["price"])

    for i in warrant.everyday_price:
        real_plot.append(i["price"])

    real_plot.remove(real_plot[0])

    for i in range(len(real_plot)):
        bs_lag.append(abs(bs_plot[i] - real_plot[i]))
        nw_lag.append(abs(nw_plot[i] - real_plot[i]))
        bsda_lag.append(abs(bsda_plot[i] - real_plot[i]))
        ukhov_lag.append(abs(ukhov_plot[i] - real_plot[i]))


#    plt.plot(range(len(real_plot)), real_plot, range(len(bs_plot)), bs_plot)
    plt.plot(range(len(real_plot)), real_plot,
             range(len(bs_plot)), bs_plot,
             range(len(nw_plot)), nw_plot,
             range(len(bsda_plot)), bsda_plot,
             range(len(ukhov_plot)), ukhov_plot,
             range(len(bs_lag)), bs_lag,
             range(len(nw_lag)), nw_lag,
             range(len(bsda_lag)), bsda_lag,
             range(len(ukhov_lag)), ukhov_lag)
    plt.ylabel("price")
    plt.show()

def get_plot(code):
    for warrant in WARRANT_LIST:
        if warrant.code == code:
            (bs_result, nw_result,
             bsda_result, ukhov_result) = calculate(warrant)
            plot(warrant, bs_result, nw_result, bsda_result, ukhov_result)

def plot_sigma(warrant):
    sigma_list = numpy.linspace(0, 0.5, num=101)
    result_list = []
    stock_match = []
    for day in warrant.target_stock.everyday_price:
        if (day["Date"] - warrant.end_date).days < 0 and \
           (warrant.start_date - day["Date"]).days <= 365:
            stock_match.append(day)

    stock_sort(stock_match)

    for sigma in sigma_list:
        t = (warrant.end_date - stock_match[300]["Date"]).days
        d_1 = (numpy.log(stock_match[300]["price"] / warrant.price) + \
               (R + 0.5 * sigma * sigma) * t) / (sigma * math.pow(t, 0.5))
        d_2 = d_1 - sigma * math.pow(t, 0.5)
        result = stock_match[300]["price"] * scipy.stats.norm.cdf(d_1) \
                     - warrant.price * math.exp(-R * t) * \
                     scipy.stats.norm.cdf(d_2)
        result_list.append(result)

    print(stock_match[1]["price"])
    plt.plot(sigma_list, result_list)
    plt.show()


if __name__ == '__main__':
    input_all()
    get_plot('031005')
    #plot_sigma(WARRANT_LIST[0])
