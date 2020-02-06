# Дан файл с данными по генной терапии - genetherapy.csv.
# Задача: Сравнить эффективность четырех различных типов терапии, представленных в таблице.

import csv
import pandas as pd
import statistics
import scipy.stats as sp


class Anova:

    def __init__(self, file_for_analyze):
        self.groups = {}
        self.file = file_for_analyze
        self.sum_mean = 0
        self.value_ssb, self.value_ssw, self.sst = 0, 0, 0

    def open_file(self):
        with open(self.file, 'r', newline='') as csv_file:
            reader = csv.DictReader(csv_file, delimiter=',')
            for val in reader:
                if val["Therapy"] in self.groups:
                    self.groups[val["Therapy"]]["mean"] += [(int(val["expr"]))]
                else:
                    self.groups[val["Therapy"]] = {'mean': [(int(val["expr"]))]}
        print(self.groups)
        self.calculate_mean()

    def calculate_mean(self):
        total_mean, n = 0, 0
        for group, values in self.groups.items():
            summ = sum(values["mean"])
            lenght = len(values["mean"])
            total_mean += summ
            n += lenght
            self.groups[group] = {'df': lenght}
            self.groups[group]["mean"] = summ / lenght
        print(self.groups)
        self.ssb(total_mean / n)

    def ssb(self, mean_gr):
        for i in self.groups.values():
            self.value_ssb += i["df"] * ((i["mean"] - mean_gr) ** 2)

    def ssw(self):
        pass

    def sst(self):
        pass

    def f_value(self):
        pass


if __name__ == '__main__':
    my_statistic = Anova(file_for_analyze='genetherapy.csv')
    my_statistic.open_file()

# ===== on PANDAS =====
# data = pd.read_csv('genetherapy.csv', sep=',')
#
# groups = pd.unique(data.Therapy.values)
# dict_data = {group: data['expr'][data.Therapy == group] for group in groups}
#
# F_value, p_value = sp.f_oneway(dict_data["A"], dict_data["B"], dict_data["C"], dict_data["D"])
# print(F_value, p_value)

# TODO задействовать collections для работы со словарями
# TODO после освоения расчетов статистических перенести работу на статистические пакеты
