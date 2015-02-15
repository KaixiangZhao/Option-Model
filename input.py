import xlrd
import time
import datetime

WARRANT_LIST = []

class Warrant:
    def __init__(self, code, target_stock, species, start_date, end_date,
                 price, portion, amount):
        self.code = code
        self.target_stock = target_stock
        self.species = species
        self.start_date = start_date
        self.end_date = end_date
        self.price = price
        self.portion = portion
        self.amount = amount

    def print_information(self):
        print('{} \t {} \t {} \t {} \t {} \t {} \t {} \t {}'.format(
            self.code, self.target_stock, self.species,
            self.start_date, self.end_date, self.price, self.portion,
            self.amount
        ))

def input_warrant():
    excel_warrant = xlrd.open_workbook("Warrant.xls")
    sheet_warrant = excel_warrant.sheets()[0]

    excel_m = xlrd.open_workbook("M.xls")
    sheet_m = excel_m.sheets()[0]

    for i in range(len(sheet_warrant.col_values(0)) - 1):
        if sheet_warrant.row_values(i + 1)[4] == "E":
            start_time = time.strptime(sheet_warrant.row_values(i + 1)[6],
                                       "%Y-%m-%d")
            end_time = time.strptime(sheet_warrant.row_values(i + 1)[7],
                                     "%Y-%m-%d")
            for j in range(len(sheet_m.col_values(0)) - 3):
                if sheet_warrant.row_values(i + 1)[0] == \
                   sheet_m.row_values(j + 3)[0]:
                    amount = sheet_m.row_values(j + 3)[2]
            WARRANT_LIST.append(
                Warrant(sheet_warrant.row_values(i + 1)[0],
                        sheet_warrant.row_values(i + 1)[2],
                        sheet_warrant.row_values(i + 1)[5],
                        datetime.datetime(start_time[0], start_time[1],
                                          start_time[2]),
                        datetime.datetime(end_time[0], end_time[1],
                                          end_time[2]),
                        sheet_warrant.row_values(i + 1)[9],
                        sheet_warrant.row_values(i + 1)[10],
                        amount))

def print_all():
    for i in WARRANT_LIST:
        i.print_information()
