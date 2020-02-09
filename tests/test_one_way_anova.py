import unittest
from unittest.mock import Mock
from one_way_anova import Anova


class GlobalCaesarTest(unittest.TestCase):

    def setUp(self):
        self.one_way = Anova('simple.csv')
        print(f'Вызван {self.shortDescription()}', flush=True)

    def tearDown(self):
        print(f'Результаты будут прологированы, но потом :)')

    def test_simple(self):
        pass