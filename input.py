import xlrd
import time
import datetime

WARRANT_LIST = []
STOCK_LIST = []

R = 8.09e-5

class Warrant:
    def __init__(self, code, target_stock, species, start_date, end_date,
                 price, portion, amount, everyday_price):
        self.code = code
        self.target_stock = target_stock
        self.species = species
        self.start_date = start_date
        self.end_date = end_date
        self.price = price
        self.portion = portion
        self.amount = amount
        self.everyday_price = everyday_price

    def print_information(self):
        print('{} \t {} \t {} \t {} \t {} \t {} \t {} \t {}'.format(
            self.code, self.target_stock.code, self.species,
            self.start_date, self.end_date, self.price, self.portion,
            self.amount
        ))

    def print_everyday_price(self):
        print('{} \t {}'.format(self.code, self.everyday_price))

class Stock:
    def __init__(self, code, amount, everyday_price):
        self.code = code
        self.amount = amount
        self.everyday_price = everyday_price

    def print_information(self):
        print('{} \t {}'.format(self.code, self.amount))

    def print_everyday_price(self):
        print('{} \t {}'.format(self.code, self.everyday_price))

def input_stock():
    excel_stock = xlrd.open_workbook("Stock.xls")
    sheet_stock = excel_stock.sheets()[0]

    excel_n = xlrd.open_workbook("N.xls")
    sheet_n = excel_n.sheets()[0]

    for i in range(len(sheet_stock.col_values(0)) - 3):
        temp = {}
        is_already_in = False
        day_time = time.strptime(sheet_stock.row_values(i + 3)[1],
                                 "%Y-%m-%d")
        temp["Date"] = datetime.datetime(
            day_time[0], day_time[1], day_time[2])
        temp["price"] = sheet_stock.row_values(i + 3)[2]
        for j in range(len(sheet_n.col_values(0)) - 3):
            if sheet_stock.row_values(i + 3)[0] == \
               sheet_n.row_values(j + 3)[0]:
                amount = sheet_n.row_values(j + 3)[2] * 1000

        for j in STOCK_LIST:
            if sheet_stock.row_values(i + 3)[0] == j.code:
                j.everyday_price.append(temp)
                is_already_in = True
                break

        if is_already_in == False:
            STOCK_LIST.append(Stock(sheet_stock.row_values(i + 3)[0],
                                    amount,
                                    [temp]))

def input_warrant():
    excel_warrant = xlrd.open_workbook("Warrant.xls")
    sheet_warrant = excel_warrant.sheets()[0]

    excel_m = xlrd.open_workbook("M.xls")
    sheet_m = excel_m.sheets()[0]

    excel_warrant_real = xlrd.open_workbook("Warrant_real.xls")
    sheet_warrant_real = excel_warrant_real.sheets()[0]

    for i in range(len(sheet_warrant.col_values(0)) - 1):
        if sheet_warrant.row_values(i + 1)[4] == "E":
            start_time = time.strptime(
                sheet_warrant.row_values(i + 1)[6],
                "%Y-%m-%d")
            end_time = time.strptime(sheet_warrant.row_values(i + 1)[7],
                                     "%Y-%m-%d")
            for j in range(len(sheet_m.col_values(0)) - 3):
                if sheet_warrant.row_values(i + 1)[0] == \
                   sheet_m.row_values(j + 3)[0]:
                    amount = sheet_m.row_values(j + 3)[2]

            everyday_price = []
            for j in range(len(sheet_warrant_real.col_values(0)) - 3):
                if sheet_warrant_real.row_values(j + 3)[0] == \
                   sheet_warrant.row_values(i + 1)[0]:
                    temp = {}
                    day_time = time.strptime(
                        sheet_warrant_real.row_values(i + 3)[2],
                        "%Y-%m-%d")
                    temp["Date"] = datetime.datetime(
                        day_time[0], day_time[1], day_time[2])
                    temp["price"] = sheet_warrant_real.row_values(i+3)[3]
                    everyday_price.append(temp)

            for stock in STOCK_LIST:
                if stock.code == sheet_warrant.row_values(i + 1)[2]:
                    target_stock = stock
            WARRANT_LIST.append(
                Warrant(sheet_warrant.row_values(i + 1)[0],
                        target_stock,
                        sheet_warrant.row_values(i + 1)[5],
                        datetime.datetime(start_time[0], start_time[1],
                                          start_time[2]),
                        datetime.datetime(end_time[0], end_time[1],
                                          end_time[2]),
                        sheet_warrant.row_values(i + 1)[9],
                        sheet_warrant.row_values(i + 1)[10],
                        amount,
                        everyday_price))

def print_warrant_information():
    for i in WARRANT_LIST:
        i.print_information()

def print_warrant_price():
    for i in WARRANT_LIST:
        i.print_everyday_price()

def print_stock_information():
    for i in STOCK_LIST:
        i.print_information()

def print_stock_price():
    for i in STOCK_LIST:
        i.print_everyday_price()

def input_all():
    input_stock()
    input_warrant()

if __name__ == '__main__':
    input_stock()
    input_warrant()
