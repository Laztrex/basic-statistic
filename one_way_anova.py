# Дан файл с данными по генной терапии - genetherapy.csv.
# Задача: Сравнить эффективность четырех различных типов терапии, представленных в таблице.

import csv
import pandas as pd
import statistics


class Anova:

    def __init__(self, file_for_analyze):
        self.groups = {}
        self.file = file_for_analyze

    def open_file(self):
        with open(self.file, 'r', newline='') as csv_file:
            reader = csv.DictReader(csv_file, delimiter=',')
            for val in reader:
                if val["Therapy"] in self.groups:
                    self.groups[val["Therapy"]].append(int(val["expr"]))
                else:
                    self.groups[val["Therapy"]] = []
        self.calculate()

    def calculate(self):
        print(statistics.mean(self.groups["A"]))
        print(statistics.mean(self.groups["B"]))
        print(statistics.mean(self.groups["C"]))
        print(statistics.mean(self.groups["D"]))
        print(f'Вся информация о таблице: {self.groups}')


if __name__ == '__main__':
    my_statistic = Anova(file_for_analyze='genetherapy.csv')
    my_statistic.open_file()

# ===== on PANDAS =====
# data = pd.read_csv('genetherapy.csv', sep=',')
# print(data["Therapy"])
#
# print(pd.unique(data.Therapy.values))
