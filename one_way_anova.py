# Дан файл с данными по генной терапии - genetherapy.csv.
# Задача: Сравнить эффективность четырех различных типов терапии, представленных в таблице.

import csv
from math import sqrt
import pandas as pd
import statistics
import scipy.stats as sp
from decimal import *
import collections


class Anova:

    def __init__(self, file_for_analyze, indep_var="Therapy", dep_var="expr"):
        getcontext().prec = 35
        self.groups = {}
        self.file = file_for_analyze
        self.indep_var = indep_var
        self.dep_var = dep_var
        self.sum_mean = 0
        self.ssb, self.ssw, self.sst, self.f, self.p = Decimal(0), Decimal(0), Decimal(0), Decimal(0), 0

    def run(self):
        self.open_file()
        self.calculate()

    def open_file(self):
        with open(self.file, 'r', newline='') as csv_file:
            reader = csv.DictReader(csv_file, delimiter=',')
            self.writer(reader)

    def writer(self, data):
        for val in data:
            if val[self.indep_var] in self.groups:
                self.groups[val[self.indep_var]]["mean"] += [Decimal((val[self.dep_var]))]
            else:
                self.groups[val[self.indep_var]] = {'mean': [Decimal(val[self.dep_var])]}

    def calculate(self):
        total_mean, n = 0, 0
        for group, values in self.groups.items():
            summ = Decimal(sum(values["mean"]))
            lenght = len(values["mean"])
            total_mean += summ
            n += lenght
            self.groups[group] = {'df': lenght}
            mean = Decimal(summ / lenght)
            self.groups[group]["mean"] = mean
            self.groups[group].update({'sd': Decimal(sqrt(self.calc_ssw(values["mean"], mean) / lenght - 1))})

        print(self.groups)
        print(self.calc_ssb(mean_gr=total_mean / n))
        print(self.ssw)
        print(self.ssw / (n - len(self.groups)))
        print(self.f_value(n))
        print(self.p_value(int(self.f), len(self.groups) - 1, n - len(self.groups)))
        self.beautiful_made_table()

    def calc_ssb(self, mean_gr):
        for i in self.groups.values():
            self.ssb += i["df"] * ((i["mean"] - mean_gr) ** 2)
        return f'ssb - {self.ssb}'

    def calc_ssw(self, values, mean_group):
        p = 0
        for i in values:
            val = (i - mean_group) ** 2
            p += val
            self.ssw += val
        return p

    def p_value(self, f, dfb, dfw):
        self.p = sp.f.sf(f, dfb, dfw)
        return f'p-value - {self.p}'

    def f_value(self, n):
        self.f = (self.ssb / (len(self.groups) - 1)) / (self.ssw / (n - len(self.groups)))
        return f'f-value = {self.f}'

    def beautiful_made_table(self):
        print(f'{"+":-<9}{"+":-<5}{"+":-<20}{"+":-<6}{"+":->{24}}\n'
              f'| {"группа":^} | {"df":^} | {"mean":^{17}} | {"sd":^{26}} |\n'
              f'{"+":-<9}{"+":-<5}{"+":-<20}{"+":-<6}{"+":->{24}}')
        for i, j in self.groups.items():
            print(f'| {i:^6} | {j["df"]:^{2}} | '
                  f'{j["mean"].quantize(Decimal("1.000")):^{17}} | {j["sd"].quantize(Decimal("1.000")):^{26}} |')
        print(f'{"+":-<9}{"+":-^}{"+":->{20}}{"+":->{34}}\n'
              f'| {"f-value":^}| {self.f.quantize(Decimal("1.000")):^} | {"p-value":^} | {round(self.p, 4):^}|\n'
              f'{"+":-<9}{"+":-^}{"+":->{20}}{"+":->{34}}')

    def gracefully_with_stat_packages(self):
        """on pandas"""
        data = pd.read_csv('genetherapy.csv', sep=',')

        groups = pd.unique(data.Therapy.values)
        dict_data = {group: data['expr'][data.Therapy == group] for group in groups}

        F_value, p_value = sp.f_oneway(dict_data["A"], dict_data["B"], dict_data["C"], dict_data["D"])
        print(F_value, p_value)


class MultiAnova(Anova):

    def __init__(self, file_for_analyze, *args):
        super().__init__(file_for_analyze=file_for_analyze)
        self.indep_var = args

    def run(self):
        self.open_file()
        print(self.groups)
        for i, j in self.groups.items():
            print(i, j)

    def writer(self, data):
        for val in data:
            if val[self.indep_var[0]] in self.groups and val[self.indep_var[1]] in self.groups:
                self.groups[val[self.indep_var[0]]][val[self.indep_var[1]]]["mean"] += ([Decimal(val[self.dep_var])])
            else:
                try:
                    self.groups[val[self.indep_var[0]]][val[self.indep_var[1]]]['mean'] += [Decimal(val[self.dep_var])]
                except KeyError:
                    self.groups[val[self.indep_var[0]]] = {
                        val[self.indep_var[1]]: {'mean': [Decimal(val[self.dep_var])]}}
                #     if val[self.indep_var[1]] in self.groups[self.indep_var[0]]:
                #         self.groups[val[self.indep_var[0]]].update \
                #             ({val[self.indep_var[1]]: {'mean': [Decimal(val[self.dep_var])]}})
                # else:
                #     self.groups[val[self.indep_var[0]]] = {
                #         val[self.indep_var[1]]: {'mean': [Decimal(val[self.dep_var])]}}


if __name__ == '__main__':
    # my_statistic = Anova(file_for_analyze='genetherapy.csv')
    # my_statistic.run()
    my = MultiAnova('atherosclerosis.csv', 'age', 'dose')
    my.run()

# TODO после освоения расчетов статистических перенести работу на статистические пакеты
# TODO Надо каждый метод сделать независимым, масштабируемым. С возможностью применять многофакторный анализ
