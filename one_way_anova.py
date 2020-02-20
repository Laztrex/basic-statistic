# Дан файл с данными по генной терапии - genetherapy.csv.
# Задача: Сравнить эффективность четырех различных типов терапии, представленных в таблице.

import csv
from math import sqrt
import pandas as pd
import statistics
import scipy.stats as sp
from decimal import *
from termcolor import cprint
import collections
from collections import Counter
from pandas import read_csv


class Anova:

    def __init__(self, file_for_analyze, indep_var="Therapy", dep_var="expr"):
        getcontext().prec = 35
        self.groups = {}
        self.file = file_for_analyze
        self.indep_var = indep_var
        self.dep_var = dep_var
        self.sum_mean = 0
        self.multi = False
        self.dict_set = []
        self.ssb, self.ssw, self.sst, self.f, self.p = Decimal(0), Decimal(0), Decimal(0), Decimal(0), 0

    def run(self):
        self.open_file()
        if not self.multi:
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

    def calculate(self, subtree=None):
        if subtree is None:
            subtree = self.groups
        total_mean, n = 0, 0

        for group, values in subtree.items():
            if self.multi:
                summ = Decimal(sum(values))
                lenght = len(values)
                total_mean += summ
                n += lenght
                mean = summ / lenght
                subtree[group] = mean
            else:
                summ = Decimal(sum(values["mean"]))
                lenght = len(values["mean"])
                total_mean += summ
                n += lenght
                subtree[group] = {'df': lenght}
                mean = summ / lenght
                subtree[group]["mean"] = mean
                subtree[group].update({'sd': Decimal(sqrt(self.calc_ssw(values["mean"], mean) / lenght - 1))})

        self.total_mean = total_mean / n
        print(subtree)
        print(self.calc_ssb(subtree=subtree, mean_gr=total_mean / n))
        print(f'mean Sq ssb - {self.ssb / (len(subtree) - 1)}')
        print(f'ssw - {self.ssw}')
        print(f'Mean Sq ssw - {self.ssw / (n - len(subtree))}')
        # print(self.f_value(subtree, n))
        # print(self.p_value(int(self.f), len(subtree) - 1, n - len(subtree)))
        # self.beautiful_made_table()

    def calc_ssb(self, subtree, mean_gr):
        for i in subtree.values():
            self.ssb += i["df"] * ((i["mean"] - mean_gr) ** 2)
        return f'ssb - {self.ssb}'

    def calc_ssw(self, values, mean_group, mean_total=None):
        p = 0
        if mean_total:
            pass
        else:
            for i in values:
                val = (i - mean_group) ** 2
                p += val
                self.ssw += val
            return p

    def p_value(self, f, dfb, dfw):
        self.p = sp.f.sf(f, dfb, dfw)
        return f'p-value - {self.p}'

    def f_value(self, subtree, n):
        self.f = (self.ssb / (len(subtree) - 1)) / (self.ssw / (n - len(subtree)))
        return f'f-value = {self.f}'

    # def beautiful_made_table(self):
    #     print(f'{"+":-<9}{"+":-<5}{"+":-<20}{"+":-<6}{"+":->{24}}\n'
    #           f'| {"группа":^} | {"df":^} | {"mean":^{17}} | {"sd":^{26}} |\n'
    #           f'{"+":-<9}{"+":-<5}{"+":-<20}{"+":-<6}{"+":->{24}}')
    #     for i, j in self.groups.items():
    #         print(f'| {i:^6} | {j["df"]:^{2}} | '
    #               f'{j["mean"].quantize(Decimal("1.000")):^{17}} | {j["sd"].quantize(Decimal("1.000")):^{26}} |')
    #     print(f'{"+":-<9}{"+":-^}{"+":->{20}}{"+":->{34}}\n'
    #           f'| {"f-value":^}| {self.f.quantize(Decimal("1.000")):^} | {"p-value":^} | {round(self.p, 4):^}|\n'
    #           f'{"+":-<9}{"+":-^}{"+":->{20}}{"+":->{34}}')

    def gracefully_with_stat_packages(self):
        """on pandas"""
        data = pd.read_csv('genetherapy.csv', sep=',')

        groups = pd.unique(data.Therapy.values)
        dict_data = {group: data['expr'][data.Therapy == group] for group in groups}

        F_value, p_value = sp.f_oneway(dict_data["A"], dict_data["B"], dict_data["C"], dict_data["D"])
        print(F_value, p_value)


class MultiAnova:

    def __init__(self, file_for_analyze, dep, repeat=3):
        self.file_with_data = file_for_analyze
        self.dep_var = dep
        self.factors = None
        self.data = None
        self.indept_list = []
        self.df_a, self.df_b = 0, 0
        self.group = collections.defaultdict(list)
        self.repeat = repeat
        self.dispers_val = []
        self.ssa_ssb = 0
        self.ssw_val = 0

    def run(self):
        self.open_file()
        self.calculate()

    def open_file(self):
        self.data = read_csv(self.file_with_data)

    def consecutive(self, my_iter, group):
        """Добавляет в словарь среднее для каждой из групп"""
        copy_list = self.indept_list.copy()
        copy_list.remove(group)
        for i in my_iter:
            if i == 1 or i == 2:
                self.group[(group, i)] += [
                    self.data[self.data[group] == i]
                    [self.data[self.data[group] == i][copy_list[0]] == d][self.dep_var].mean()
                    for d in self.data[self.data[group] == i][copy_list[0]]
                ]
        # TODO упростить выражение

    def calculate(self):
        self.indept_list = list(self.data.keys())
        self.indept_list.remove(self.dep_var)
        self.factors = dict(map(self.ssx, self.indept_list))
        print(self.df_a, self.df_b)
        for i in self.indept_list:
            self.consecutive(my_iter=self.data[i].unique(), group=i)

        if self.repeat == 1:
            self.ssw_val = self.ss_versus(*self.factors.values(), self.sst(), 0)
            self.repeat += 1
            self.dispers_val = (self.dispersia(*self.factors.values(), self.ssw_val, 0))
        else:
            self.ssw_val = self.ssw()
            self.ssa_ssb = self.ss_versus(*self.factors.values(), self.sst(), self.ssw_val)
            self.dispers_val = (self.dispersia(*self.factors.values(), self.ssa_ssb, self.ssw_val))

        print(f'ssw {self.ssw_val}')
        print(f'ssa_ssb {self.ssa_ssb}')
        print(self.dispers_val)
        print(self.f_values(*self.dispers_val))

    def ssx(self, group_mean):
        """вычисление SSa (объяснённая влиянием фактора A сумма квадратов отклонений)
         и SSb (объяснённая влиянием фактора B сумма квадратов отклонений)"""
        if not self.df_a:
            self.df_a = len(self.data[group_mean].unique()) - 1
        else:
            self.df_b = len(self.data[group_mean].unique()) - 1

        ssb = sum([(self.data[self.data[group_mean] == j][self.dep_var].mean() -
                    self.data[self.dep_var].mean()) ** 2 for j in self.data[group_mean]])
        print(f'SSB ({group_mean}) - {ssb}')
        return group_mean, ssb

    def sst(self):
        """Общая сумма SS (квадратов отклонений)"""
        sst = sum((self.data[self.dep_var] - self.data[self.dep_var].mean()) ** 2)
        print(f'SS - {sst}')
        return sst

    def ss_versus(self, ssa, ssb, sst, ssq_w):
        """SSab - объяснённая влиянием взаимодействия факторов A и B сумма квадратов отклонений"""
        return sst - ssa - ssb - ssq_w

    def ssw(self):
        """Необъяснённая сумма квадратов отклонений или сумма квадратов отклонений ошибки"""
        s = 0
        for i, j in self.group.items():
            s += sum((self.data[self.data[i[0]] == i[1]][self.dep_var] - j) ** 2)
        return s

    def dispersia(self, ssa, ssb, ssab, ssq):
        return ssa / self.df_a, ssb / self.df_b, ssab / (self.df_a * self.df_b), ssq / (
                    (self.df_a + 1) * (self.df_b + 1) * (self.repeat - 1))

    def f_values(self, msa, msb, msab, msq):
        m = len(self.data[self.dep_var].values) - len(self.data.keys())
        dell = msq if msq else msab / m
        return msa / dell, msb / dell, msab / dell


if __name__ == '__main__':
    # my_statistic = Anova(file_for_analyze='genetherapy.csv')
    # my_statistic.run()
    my = MultiAnova('test_sample.csv', 'expr', 3)
    my.run()

# TODO после освоения расчетов статистических перенести работу на статистические пакеты
