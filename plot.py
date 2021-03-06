import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pylab
import numpy
import math
import scipy.stats
from calculate import WARRANT_LIST, calculate, input_all, stock_sort, R, sub_ukhov, ukhov

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


    pylab.plot(range(len(real_plot)), real_plot, color='blue', linewidth=1, linestyle='-', label='real')
    pylab.plot(range(len(ukhov_plot)), ukhov_plot, color='green', linewidth=1, linestyle='-', label='bs')
    # pylab.plot(range(len(real_plot)), real_plot,
    #          range(len(bs_plot)), bs_plot,
    #          range(len(nw_plot)), nw_plot,
    #          range(len(bsda_plot)), bsda_plot,
    #          range(len(ukhov_plot)), ukhov_plot,
    #          range(len(bs_lag)), bs_lag,
    #          range(len(nw_lag)), nw_lag,
    #          range(len(bsda_lag)), bsda_lag,
    #          range(len(ukhov_lag)), ukhov_lag)
    pylab.ylabel("price")
    pylab.legend(loc='upper left')
    pylab.show()

def get_plot(code):
    for warrant in WARRANT_LIST:
        if warrant.code == code:
            (bs_result, nw_result,
             bsda_result, ukhov_result) = calculate(warrant)
            plot(warrant, bs_result, nw_result, bsda_result, ukhov_result)

def plot_sigma(warrant):
    sigma_list = numpy.linspace(0, 1, num=101)
    result_list = []
    stock_match = []
    for day in warrant.target_stock.everyday_price:
        if (day["Date"] - warrant.end_date).days < 0 and \
           (warrant.start_date - day["Date"]).days <= 0:
            stock_match.append(day)

    stock_sort(stock_match)

    for sigma in sigma_list:
        t = (warrant.end_date - stock_match[20]["Date"]).days
        d_1 = (numpy.log(stock_match[20]["price"] / warrant.price) + \
               (R + 0.5 * sigma * sigma) * t) / (sigma * math.pow(t, 0.5))
        d_2 = d_1 - sigma * math.pow(t, 0.5)
        result = stock_match[20]["price"] * scipy.stats.norm.cdf(d_1) \
                     - warrant.price * math.exp(-R * t) * \
                     scipy.stats.norm.cdf(d_2)
        result_list.append(result)

    print(stock_match[20]["price"])
    print(result_list)
    plt.plot(sigma_list, result_list)
    plt.show()

def ukhov_plot():
    w, d1, v_list, sig_list = ukhov(17.5, 7.0, 361, 8.09e-05, 0, 50000000.0,
                                    149522567.0, 0.0350827369223, 1.0)
    v_list.remove(v_list[0])
    v_list.remove(v_list[0])
    v_list.remove(v_list[0])
    v_list.remove(v_list[0])
    sig_list.remove(sig_list[0])
    sig_list.remove(sig_list[0])
    sig_list.remove(sig_list[0])
    sig_list.remove(sig_list[0])
    n = 256
    v = numpy.linspace(0.5*17.5*149522567 , 1.5*17.5*149522567, n)
    sig = numpy.linspace(0, 0.5, n)
    X, Y = numpy.meshgrid(v, sig)
    pylab.contourf(X, Y,
                   ((abs(sub_ukhov(17.5, 7.0, 361, 8.09e-05, 0, 50000000.0,
                             149522567.0, 0.0350827369223, 1.0, X, Y)[0]["V"] - X))/X)**2 +
                   abs(sub_ukhov(17.5, 7.0, 361, 8.09e-05, 0, 50000000.0,
                             149522567.0, 0.0350827369223, 1.0, X, Y)[0]["sig"] - Y)**2,
                   50, alpha=.75, cmap='jet')
    C = pylab.contour(X, Y,
                      ((abs(sub_ukhov(17.5, 7.0, 361, 8.09e-05, 0, 50000000.0,
                             149522567.0, 0.0350827369223, 1.0, X, Y)[0]["V"] - X))/X)**2 +
                      abs(sub_ukhov(17.5, 7.0, 361, 8.09e-05, 0, 50000000.0,
                             149522567.0, 0.0350827369223, 1.0, X, Y)[0]["sig"] - Y)**2,
                      50, colors='black', linewidth=.5)
    #pylab.scatter(v_list, sig_list, s=100, alpha=.75, marker=(5,2), c='black')
    pylab.clabel(C, inline=1, fontsize=10)
    pylab.show()

    # pylab.contourf(X, Y,
    #                abs(sub_ukhov(17.5, 7.0, 361, 8.09e-05, 0, 50000000.0,
    #                          149522567.0, 0.0350827369223, 1.0, X, Y)[0]["sig"] - Y),
    #                20, alpha=.75, cmap='jet')
    # CC = pylab.contour(X, Y,
    #                   abs(sub_ukhov(17.5, 7.0, 361, 8.09e-05, 0, 50000000.0,
    #                          149522567.0, 0.0350827369223, 1.0, X, Y)[0]["sig"] - Y),
    #                   20, colors='black', linewidth=.5)
    # pylab.clabel(CC, inline=1, fontsize=10)
    # pylab.show()

if __name__ == '__main__':
    input_all()
    get_plot('031001')
    #plot_sigma(WARRANT_LIST[0])
    #ukhov_plot()
