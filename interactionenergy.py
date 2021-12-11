import math
from filereader import FileReader
from vector import Vector


class InteractionEnergy(object):
    """
    Класс, в котором реализован способ рассчета энергии взаимодействия атомов
    """

    def __init__(self, filepath, D0, r0, S, Betta, gamma, c, d, h, two_mu, R, D):
        """
        Конструктор класса, в который можно передать разные параметры, используемы
        для рассчета энергии взаимодействия атомов кремния и углерода.
        """
        self.filepath = filepath
        self.D0 = D0  # eV
        self.r0 = r0  # A
        self.S = S
        self.Betta = Betta  # A^(-1)
        self.gamma = gamma
        self.c = c
        self.d = d
        self.h = h
        self.two_mu = two_mu  # A^(-1)
        self.R = R  # A
        self.D = D  # A

        filreader = FileReader(filepath)
        self.atom_radiuses, _, self.atoms = filreader.file_reader_mod()

    def r(self, i, j):
        """
        Вызывает метод класса Point distance для подсчета модуля радиус-вектора между двумя точками
        Имеющими индексы i и j
        """
        vector = Vector(self.atoms[i], self.atoms[j])
        rij = vector.lenght()
        return rij

    def Vr(self, r):
        """
        Функция для подсчета энергии отталкивания
        """
        Vrepulsion = self.D0 * math.exp((-1) * self.Betta * math.sqrt(2 * self.S) * (r - self.r0)) / (self.S - 1)
        return Vrepulsion

    def Va(self, r):
        """
        Функция для подсчета энергии притяжения
        """
        Vattraction = self.S * self.D0 * math.exp((-1) * self.Betta * math.sqrt(2 / self.S) * (r - self.r0)) / (self.S - 1)
        return Vattraction

    def fc(self, r):
        """
        Функция расстояния между атомами, имеющая смысл зарезания вазимодейсвия на далеких расстояних
        В небольшой области сглаживает как синус
        """
        if r < self.R - self.D:
            return 1
        elif math.fabs(self.R - r) <= self.D:
            return 0.5 - 0.5 * math.sin(math.pi * (r - self.R) / (2 * self.D))
        else:
            return 0

    def cos_teta(self, i, j, k):
        """
        Косинус угла тета, угла между двумя векторами, идущими из атома i к атомам j и k
        """
        vector1 = Vector(self.atoms[i], self.atoms[j])
        vector2 = Vector(self.atoms[i], self.atoms[k])
        return vector1.cosine(vector2)

    def g(self, cos_teta):
        """

        """
        return self.gamma * (1 + self.c**2 / self.d**2 - self.c**2 / (self.d**2 + (self.h + cos_teta)**2))

    def X(self, i, j):
        """
        Угловая функция
        """
        X = 0
        for k in range(0, len(self.atoms)):
            if k != i and k != j:
                rik = self.r(i, k)
                rij = self.r(i, j)
                Xij = self.fc(rik) * math.exp(self.two_mu * (rij - rik)) * self.g(self.cos_teta(i, j, k))

                X = X + Xij
                Xij = 0

        return X

    def b(self, i, j):
        """
        Порядок связей определяется выражением
        """
        return 1 / math.sqrt(1 + self.X(i, j))

    def calculate_energy(self):
        E = 0
        for i in range(0, len(self.atoms)):
            for j in range(0, len(self.atoms)):
                if i > j:
                    rij = self.r(i, j)
                    mean_b = (self.b(i, j) + self.b(j, i)) / 2
                    Eij = self.fc(rij) * (self.Vr(rij) - mean_b * self.Va(rij))

                    E = E + Eij
                    Eij = 0
        return E


def main():
    PATH = "/Users/magomednikolaev/Documents/Вычислительные_и_инф_технологии/Interaction_energy/SiMolecule.xyz"

    IE = InteractionEnergy(filepath=PATH,
                           D0=3.24,
                           r0=2.222,
                           S=1.57,
                           Betta=1.4760,
                           gamma=0.09253,
                           c=1.13681,
                           d=0.63397,
                           h=0.335,
                           two_mu=0.0,
                           R=2.90,
                           D=0.15)


if __name__ == '__main__':
    main()
