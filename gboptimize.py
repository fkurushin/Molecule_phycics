import math
import numpy as np
from tqdm import tqdm
from point import Point
from vector import Vector
from filereader import FileReader
from interactionenergy import InteractionEnergy


class GBOptimize(InteractionEnergy):
    """Класс, оптимизирующий структуру молекулы методом градиентоного спуска"""

    def __init__(self, filepath, D0, r0, S, Betta, gamma, c, d, h, two_mu, R, D, delta, precision):
        """
        Конструктор класса, в который можно передать разные параметры, используемы
        для рассчета энергии взаимодействия атомов кремния и углерода.
        """
        InteractionEnergy.__init__(self, filepath, D0, r0, S, Betta, gamma, c, d, h, two_mu, R, D)
        self.delta = delta
        self.precision = precision
        self.forces = list()

    def add_new(self, atoms):
        self.atoms = atoms

    def derivative(self, axis, atom_i):
        """
        atom_i - номер атома, который будем толкать
        """
        Uq = self.calculate_energy()
        self.atoms[atom_i].move(axis, self.delta)
        Uq_plus_dq = self.calculate_energy()
        self.atoms[atom_i].move(axis, -self.delta)
        return (Uq_plus_dq - Uq) / self.delta

    def delta_sign(self, axis, atom_i):
        """
        Выбирает направление наискорейшего спуска
        """
        if self.derivative(axis, atom_i) > 0:
            return -1  # Возвращаем -1 чтобы шли в сторону уменьшения функции
        else:
            return 1

    def find_min(self, ):
        """
        """
        axises = ['x', 'y', 'z']
        print("\nCalculating optimal displacment...")
        for i in tqdm(range(len(self.atoms))):
            for axis in axises:
                sign = self.delta_sign(axis, i)
                while self.derivative(axis, i) > self.precision:
                    self.atoms[i].move(axis, sign * self.delta)

    def print_atoms(self, ):
        """
        Метод для распечатки атомов
        """
        for i in range(len(self.atoms)):
            self.atoms[i].print_point()

    def grad(self, atom_i):
        """
        Возвращает градиент функции
        """
        return self.derivative('x', atom_i), self.derivative('y', atom_i), self.derivative('z', atom_i)

    def find_forces(self, ):
        """
        """
        forces = list()
        for i in range(len(self.atoms)):
            fx, fy, fz = self.grad(i)
            force = Point(fx, fy, fz)
            forces.append(force)

        self.forces = forces


def main():
    PATH = "/Users/fedorkurusin/Documents/informatics/GBoost/molecules/si.xyz"

    GBO = GBOptimize(filepath=PATH,
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
                     D=0.15,
                     delta=0.0001,
                     precision=0.001)
    print("До:")
    GBO.print_atoms()
    print(GBO.calculate_energy())

    GBO.find_min()

    print("\nПосле:")
    GBO.print_atoms()
    print(GBO.calculate_energy())


if __name__ == '__main__':
    main()
