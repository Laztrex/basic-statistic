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
from itertools import chain


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
                # subtree[group].update({'sd': Decimal(sqrt(self.calc_ssw(values["mean"], mean) / lenght - 1))})

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
        pass
        # p = 0
        # if mean_total:
        #     pass
        # else:
        #     for i in values:
        #         val = (i - mean_group) ** 2
        #         p += val
        #         self.ssw += val
        #     return p

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


class MultiAnova(Anova):

    def __init__(self, file_for_analyze, dep, repeat=1):
        self.repeat = repeat
        super().__init__(file_for_analyze=file_for_analyze)
        self.dep_var = dep
        self.multi = True
        self.one_dict = collections.OrderedDict()
        self.two_dict = collections.OrderedDict()
        self.a, self.b = 6, 9
        self.matrix_column_dict = collections.defaultdict(list)
        self.matrix_rows_dict = collections.defaultdict(list)
        self.calc_ssa_ssb = []
        self.total_mean = 0
        self.a = self.matrix_rows_dict.values()
        self.b = self.matrix_column_dict.values()

    def run(self):
        self.open_file()

    def writer(self, data):
        indep_list = data.fieldnames.copy()
        indep_list.remove(self.dep_var)
        for val in data:
            for i in indep_list[1:]:
                self.matrix_column_dict['column' + ' ' + val[i]] += [Decimal(val[self.dep_var])]  # TODO а дальше спец методы со словарями (объединение, zip и тд), чтобы вытянуть нужную строку и опознать среднее для нее
                self.matrix_rows_dict['rows' + ' ' + val[indep_list[0]]] += [Decimal(val[self.dep_var])]
                # TODO создать 3-й словарь с пересечение групп
                # >> > Matrix[(2, 3, 4)] = 88
                # >> > Matrix[(7, 8, 9)] = 99
        self.one_dict = self.matrix_column_dict.copy()
        self.two_dict = self.matrix_rows_dict.copy()  # TODO: теперь их перемножить, скомпехешнев

        for i in (self.one_dict, self.two_dict):
            self.calculate(i)
        print(self.matrix_column_dict)
        print(self.matrix_rows_dict)
        print(self.one_dict)
        print(self.two_dict)

        for i, j in self.matrix_column_dict.items():
            print(j)
            self.matrix_column_dict[i] = list(self.group(j, len(self.matrix_column_dict)))

        self.calc_ssw_new(values={**self.matrix_column_dict, **self.matrix_rows_dict},
                          mean_group_xi=self.one_dict,
                          mean_group_xj=self.two_dict,
                          mean_total=self.total_mean)

    # def representation(self, my_dict, too_my_dict):
    #     for a in group:
    #         [self.new_dict[i].append([mean]) for i in my_dict.keys() if i.endswith(a)]

    def calc_ssb(self, subtree, mean_gr):
        self.ssb = 0
        df = len(subtree)
        for i in subtree.values():
            self.ssb += (i - mean_gr) ** 2
        if df == 2:
            return f'ssb (sum Sq) - {3 * self.repeat * self.ssb}'
        else:
            return f'ssb (sum Sq) - {2 * self.repeat * self.ssb}'

    def group(self, iterable, count):
        """ Группировка элементов последовательности по count элементов """

        return zip(*[iter(iterable)] * count)

    def calc_ssw_new(self, values, mean_group_xi, mean_group_xj, mean_total):
        cprint(values, color='blue')
        cprint(mean_group_xi, color='blue')
        cprint(mean_group_xj, color='blue')

            # if i in mean_group_xi:
            #     self.ssw += [n - mean_group_xi[i]



if __name__ == '__main__':
    # my_statistic = Anova(file_for_analyze='genetherapy.csv')
    # my_statistic.run()
    my = MultiAnova('test_sample.csv', 'expr', 3)
    my.run()

# TODO после освоения расчетов статистических перенести работу на статистические пакеты
