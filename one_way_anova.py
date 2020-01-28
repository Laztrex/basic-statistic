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
        self.means = {}
        self.sum_mean = 0

    def open_file(self):
        with open(self.file, 'r', newline='') as csv_file:
            reader = csv.DictReader(csv_file, delimiter=',')
            for val in reader:
                if val["Therapy"] in self.groups:
                    self.groups[val["Therapy"]].append(int(val["expr"]))
                    self.means[val["Therapy"]] = 0
                else:
                    self.groups[val["Therapy"]] = [int(val["expr"])]
        self.calculate_mean()

    def calculate_mean(self):
        total_mean_groups = []
        for group, values in self.groups.items():
            total_mean_groups += values
            self.means[group] = statistics.mean(values)
        self.sum_mean = sum(total_mean_groups) / len(total_mean_groups)
        print(self.sum_mean)
        print(self.means)
        print(f'Вся информация о таблице: {self.groups}')
        self.sst()

    def sst(self):
        for group, values in self.groups.items():
            print(sum([(x - self.sum_mean) ** 2 for x in values]))

    #
    # def calc_sum_sq(self):


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

