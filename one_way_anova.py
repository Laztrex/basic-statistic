# Дан файл с данными по генной терапии - genetherapy.csv.
# Задача: Сравнить эффективность четырех различных типов терапии, представленных в таблице.

import csv
from math import sqrt
import pandas as pd
import statistics
import scipy.stats as sp


class Anova:

    def __init__(self, file_for_analyze):
        self.groups = {}
        self.file = file_for_analyze
        self.sum_mean = 0
        self.value_ssb, self.value_ssw, self.sst, self.f, self.p = 0, 0, 0, 0, 0

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
            mean = summ / lenght
            self.groups[group]["mean"] = mean

            self.groups[group].update({'sd': sqrt(self.ssw(values["mean"], mean) / lenght - 1)})
        print(self.groups)
        self.ssb(total_mean / n)
        self.f = self.f_value(n)
        print(self.value_ssw)
        print(self.value_ssw / (n - len(self.groups)))
        self.beautiful_made_table()

    def ssb(self, mean_gr):
        for i in self.groups.values():
            self.value_ssb += i["df"] * ((i["mean"] - mean_gr) ** 2)
        print(f'ssb - {self.value_ssb}')

    def ssw(self, values, mean_group):
        p = 0
        for i in values:
            val = (i - mean_group) ** 2
            p += val
            self.value_ssw += val
        return p

    def p_value(self, f, dfb, dfw):
        self.p = sp.f.sf(f, dfb, dfw)
        return self.p

    def f_value(self, n):
        print(f'ssw {self.value_ssw}')
        f = (self.value_ssb / (len(self.groups) - 1)) / (self.value_ssw / (n - len(self.groups)))
        print(f'f-value - {f}')
        print(self.value_ssb / (len(self.groups) - 1))
        print(f'p-value - {self.p_value(f, len(self.groups) - 1, n - len(self.groups))}')
        return f

    def beautiful_made_table(self):
        print(f'{"+":-<9}{"+":-<5}{"+":-<20}{"+":-<6}{"+":->{24}}\n'
              f'| {"группа":^} | {"df":^} | {"mean":^{17}} | {"sd":^{26}} |\n'
              f'{"+":-<9}{"+":-<5}{"+":-<20}{"+":-<6}{"+":->{24}}')
        for i, j in self.groups.items():
            print(f'| {i:^6} | {j["df"]:^{2}} | {j["mean"]:^{17}} | {j["sd"]:^{26}} |')
        print(f'{"+":-<9}{"+":-^}{"+":->{20}}{"+":->{34}}\n'
              f'| {"f-value":^}| {self.f:^} | {"p-value":^} | {self.p:^}|\n'
              f'{"+":-<9}{"+":-^}{"+":->{20}}{"+":->{34}}')

    def gracefully_with_stat_packages(self):
        data = pd.read_csv('genetherapy.csv', sep=',')

        groups = pd.unique(data.Therapy.values)
        dict_data = {group: data['expr'][data.Therapy == group] for group in groups}

        F_value, p_value = sp.f_oneway(dict_data["A"], dict_data["B"], dict_data["C"], dict_data["D"])
        print(F_value, p_value)


class MultiAnova(Anova):

    def __init__(self, file_for_analyze):
        super().__init__(file_for_analyze=file_for_analyze)
        pass


if __name__ == '__main__':
    my_statistic = Anova(file_for_analyze='genetherapy.csv')
    my_statistic.open_file()

# ===== on PANDAS =====


# TODO задействовать collections для работы со словарями
# TODO после освоения расчетов статистических перенести работу на статистические пакеты
# TODO Надо каждый метод сделать независимым, масштабируемым. С возможностью применять многофакторный анализ
