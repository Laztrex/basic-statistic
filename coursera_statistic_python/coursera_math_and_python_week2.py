import numpy as np
import re
import scipy.spatial
import scipy.linalg

from collections import defaultdict
from math import sin, cos, exp
from matplotlib import pylab as plt
from scipy import interpolate


class Comparison:

    def __init__(self, file):
        self.data_text = file
        self.words_dict = defaultdict(int)
        self.sentences_list = []
        self.numbers_sentence = 0
        self.index = 0

    def run(self):
        self.open_file()
        self.cos_operation()

    def open_file(self):
        with open(self.data_text, mode='r', encoding='utf8') as data:
            for line in data.readlines():
                self.numbers_sentence += 1
                self.diff(line.strip().lower())

    def write_file(self):
        pass

    def diff(self, sentence):
        filter_list = list(filter(None, re.split('[^a-z]', sentence)))
        self.sentences_list.append(filter_list)
        for word in filter_list:
            if word not in self.words_dict:
                self.words_dict[word] += self.index
                self.index += 1

    def cos_operation(self):
        m = np.zeros((self.numbers_sentence, len(self.words_dict)))  # создаем массив размерность строки*слова
        print(m.shape)
        print(self.sentences_list)

        for i in range(self.numbers_sentence):
            for word in self.sentences_list[i]:
                curr_word = self.words_dict[word]
                m[i][curr_word] += 1

        distances = list()
        for i in range(self.numbers_sentence):
            distance = scipy.spatial.distance.cosine(m[0, :], m[i, :])  # считаем косинусную дистанцию
            distances.append((i, distance))
        sort = sorted(distances, key=lambda tup: tup[1])  # сортируем
        print(sort[1], sort[2])


class Approximate:

    def __init__(self):
        self.arrays = [[1, 15], [1, 8, 15], [1, 4, 8, 15], ]

    def f(self, x):
        return sin(x / 5.0) * exp(x / 10.0) + 5 * exp(-x / 2.0)

    def calculate(self):
        fx2 = np.vectorize(self.f)
        for i in self.arrays:
            p = np.array(i)
            b = fx2(p)
            a = np.zeros((len(p), len(p)))
            for j in range(0, len(i)):
                a[j, :] = np.array([p[j] ** n for n in range(0, len(i))])
            s1 = scipy.linalg.solve(a, b)
            print(s1)
            self.plot(p, b)

    def plot(self, x, y):
        if len(x) == 4:
            f = interpolate.interp1d(x, y, kind='quadratic')
            ynew = f(x)
            plt.plot(x, y, 'o', x, ynew, '-')
            plt.show()
            return
        plt.plot(x, y)
        plt.show()


if __name__ == '__main__':
    # comparison = Comparison('sentences.txt')
    # comparison.run()
    aprx = Approximate()
    aprx.calculate()
