import numpy as np
import re
import scipy.spatial
import scipy.linalg

from collections import defaultdict
from math import sin, cos, exp
from matplotlib import pylab as plt
from scipy import interpolate

# Дан набор предложений, скопированных с Википедии. Каждое из них имеет "кошачью тему" в одном из трех смыслов:
# кошки (животные) UNIX-утилита cat для вывода содержимого файлов версии операционной системы OS X, названные в честь
# семейства кошачьих.
# Задача — найти два предложения, которые ближе всего по смыслу к расположенному в самой
# первой строке. В качестве меры близости по смыслу используется косинусное расстояние.


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
        massive = np.zeros((self.numbers_sentence, len(self.words_dict)))
        print(massive.shape)
        print(self.sentences_list)

        for i in range(self.numbers_sentence):
            for word in self.sentences_list[i]:
                curr_word = self.words_dict[word]
                massive[i][curr_word] += 1

        distances = list()
        for i in range(self.numbers_sentence):
            distance = scipy.spatial.distance.cosine(massive[0, :], massive[i, :])
            distances.append((i, distance))
        sort = sorted(distances, key=lambda tup: tup[1])
        print(sort[1], sort[2])


class Approximate:

    def __init__(self):
        self.arrays = [[1, 15], [1, 8, 15], [1, 4, 10, 15], ]

    def func_given(self, x):
        return sin(x / 5.0) * exp(x / 10.0) + 5 * exp(-x / 2.0)

    def calculate(self):
        fx2 = np.vectorize(self.func_given)
        for i in self.arrays:
            dots = np.array(i)
            value = fx2(dots)
            koef = np.zeros((len(dots), len(dots)))
            for j in range(0, len(i)):
                koef[j, :] = np.array([dots[j] ** n for n in range(0, len(i))])
            solution_system = scipy.linalg.solve(koef, value)
            print(solution_system)
            self.plot(dots, value)

    def plot(self, x, y):
        if len(x) == 4:
            f = interpolate.interp1d(x, y, kind='quadratic')
            y_new = f(x)
            plt.plot(x, y, 'o', x, y_new, '-')
            plt.show()
            return
        plt.plot(x, y)
        plt.show()


if __name__ == '__main__':
    # comparison = Comparison('files/sentences.txt')
    # comparison.run()
    aprx = Approximate()
    aprx.calculate()
