from scipy.optimize import minimize, differential_evolution
import numpy as np
from math import sin, exp
from matplotlib import pyplot as plt


class MinSmoothFunc:
    def __init__(self, interval):
        self.array = interval

    def func_given(self, x):
        return sin(x / 5.0) * exp(x / 10.0) + 5 * exp(-x / 2.0)

    def h_func(self, x):
        return int(self.func_given(x))

    def default_minimize(self):
        x = np.arange(*self.array, 0.01)
        result = map(self.func_given, x)
        plt.plot(x, list(result))
        plt.show()
        min1, min2 = self._min_sum(func=self.func_given, x0s=[np.array(2), np.array(30)], method='BFGS')
        self.write_ans(f"{round(float(min1.fun), 2)} {round(float(min2.fun), 2)}", '1')

    def _min_sum(self, x0s, method, func):
        for x0 in x0s:
            yield minimize(fun=func, x0=x0, method=method)

    def global_optimize(self):
        scope = [(self.array[0], self.array[1])]
        answer = differential_evolution(self.func_given, scope)
        self.write_ans(f"{round(float(answer.fun), 2)}", '2')

    def rough_optimization(self):
        x = np.arange(*self.array, 0.01)
        result1 = map(self.func_given, x)
        result2 = map(self.h_func, x)
        plt.plot(x, list(result1))
        plt.plot(x, list(result2))
        plt.show()
        # ans1 = minimize(self.h_func, np.array(30), method='BFGS')
        temp, ans1 = self._min_sum(x0s=[np.array(2), np.array(30)], method='BFGS', func=self.h_func)
        ans2 = differential_evolution(self.h_func, bounds=[(1, 30)])

        self.write_ans(f"{round(float(ans1.fun), 2)} {round(float(ans2.fun), 2)}", '3')

    def write_ans(self, data, task):
        with open(f'out-{task}.txt', 'w') as file:
            file.writelines(data)


if __name__ == '__main__':
    task1 = MinSmoothFunc([1, 30])
    task1.default_minimize()
    task1.global_optimize()
    task1.rough_optimization()
