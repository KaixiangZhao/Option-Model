import matplotlib.pyplot as plt
from calculate import WARRANT_LIST, calculate, input_all

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

if __name__ == '__main__':
    input_all()
    get_plot('031001')
