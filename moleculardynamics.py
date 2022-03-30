import os
import math
import copy
import random
import numpy as np
import numpy.linalg as LA
from tqdm import tqdm
from multiprocessing import Pool
# from multiprocessing import Process
from point import Point
from vector import Vector
from filereader import FileReader
import matplotlib.pyplot as plt

from datetime import datetime
import time


class MolecularDynamics(object):
    """docstring for MolecularDynamics."""

    def __init__(self, filepath, D0, r0, S, Betta, gamma, c, d, h, two_mu, R, D, delta, precision):
        super(MolecularDynamics, self).__init__()

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
        _, self.atoms, _ = filreader.file_reader()

        self.delta = delta
        self.precision = precision
        self.mass = 2912 * len(self.atoms)
        self.forces = np.zeros((len(self.atoms), 3))
        # self.derivatives = np.zeros((len(self.atoms), 3))
        self.signs = np.zeros((len(self.atoms), 3))
        # self.velocities = list()

        global lenght
        lenght = len(self.atoms)

        global axises
        # axises = {0: 'x', 1: 'y', 2: 'z'}
        axises = np.array([0, 1, 2])

    """
    INTERACTIONENERGY
    """
    def r(self, i, j):
        """
        подсчета модуля радиус-вектора между двумя точками
        Имеющими индексы i и j
        """
        return np.sqrt(np.sum(np.square(self.atoms[i] - self.atoms[j])))

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
        vector1 = self.atoms[i] - self.atoms[j]
        vector2 = self.atoms[i] - self.atoms[k]
        return (vector1 @ vector2) / (LA.norm(vector1) * LA.norm(vector2))

    def g(self, cos_teta):
        """

        """
        return self.gamma * (1 + self.c**2 / self.d**2 - self.c**2 / (self.d**2 + (self.h + cos_teta)**2))

    def X(self, i, j):
        """
        Угловая функция
        """
        X = 0
        for k in range(lenght):
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

    def calculate_energy(self, ):
        E = 0
        for i in range(lenght):
            for j in range(lenght):
                if i > j:
                    rij = self.r(i, j)
                    mean_b = (self.b(i, j) + self.b(j, i)) / 2
                    Eij = self.fc(rij) * (self.Vr(rij) - mean_b * self.Va(rij))
                    # print(f'{Eij:.2f} \t {rij:.2f} \t {self.fc(rij):.2f} \t {self.Vr(rij):.2f} \t {mean_b:.2f} \t {self.Va(rij):.2f}')
                    E = E + Eij
                    Eij = 0
        return E
    """
    END
    """

    """
    MOLECULEVOLUME
    """
    def add_new(self, atoms, radiuses):
        '''
        PseudoКонструктор класса
        '''
        self.atoms = atoms
        self.atom_radiuses = radiuses

    def find_parallelepiped(self):
        '''
        Находит координаты параллелепипеда
        '''
        list_x = list()
        list_y = list()
        list_z = list()
        list_x_rad = list()
        list_y_rad = list()
        list_z_rad = list()

        for atom in self.atoms:
            list_x.append(atom.x)
            list_y.append(atom.y)
            list_z.append(atom.z)

            for x in list_x:
                for rad in self.atom_radiuses.values():
                    p = x + rad
                    p1 = x - rad
                    list_x_rad.extend((p, p1))

            for y in list_y:
                for rad in self.atom_radiuses.values():
                    p = y + rad
                    p1 = y - rad
                    list_y_rad.extend((p, p1))

            for z in list_z:
                for rad in self.atom_radiuses.values():
                    p = z + rad
                    p1 = z - rad
                    list_z_rad.extend((p, p1))

        minmax = dict()
        minmax['x'] = [min(list_x_rad), max(list_x_rad)]
        minmax['y'] = [min(list_y_rad), max(list_y_rad)]
        minmax['z'] = [min(list_z_rad), max(list_z_rad)]
        return minmax

    def parallelepiped_volume(self, dict):
        '''
        Считает объем параллелепипеда
        '''
        v = (dict['x'][1] - dict['x'][0]) * (dict['y'][1] - dict['y'][0]) * (dict['z'][1] - dict['z'][0])
        return v

    def random_point_in_paral(self, dict):
        '''
        Создает рандомную точку в параллелепипеде
        '''
        x = random.uniform(dict['x'][0], dict['x'][1])
        y = random.uniform(dict['y'][0], dict['y'][1])
        z = random.uniform(dict['z'][0], dict['z'][1])
        return Point(x, y, z)

    def calculate_volume(self, precision, n_experiments, alpha, verbose=0):
        '''
        Расчитывает объем молекулы методом Монте Карло при заданной точности, считает погрешность,
        Можно регулировать количество измеренией
        '''
        volumes = list()
        parall_cords = self.find_parallelepiped()
        if verbose == 1:
            print("\nCalculating volume...")
            dis = False
        else:
            dis = True
        for i in tqdm(range(n_experiments), disable=dis):

            n_accuracy = 0
            ratio = 0.0
            n_accuracy = 0.0
            for i in range(precision):

                atom_counter = 0
                p_random = self.random_point_in_paral(parall_cords)
                for atom in self.atoms:

                    if p_random.is_in_atom(atom) is True:
                        atom_counter += 1

                if atom_counter == 0:
                    n_accuracy += 0
                elif atom_counter == 1:
                    n_accuracy += 1
                else:
                    n_accuracy += 1
                    '''
                    Я считаю, что законмментированный вариант правильный, но тот, что выше
                    совпадает с овтетом
                    '''
                    # n_accuracy += 1 / atom_counter

            n_total = precision
            ratio = n_accuracy / n_total
            volume = self.parallelepiped_volume(parall_cords) * ratio
            volumes.append(volume)

        volumes = np.array(volumes)
        # volumes, _, volume_std = stats.bayes_mvs(volumes, alpha)

        return round(volumes.mean(), 2), round(volumes.std(), 2)

    """
    END
    """

    """
    RANDNORMAL
    """
    def get_rand_num_arr(self, mu=0, sigma=1, number_el=int(10e5)):
        return np.random.normal(mu, sigma, number_el)

    def show_momentum(self, x, mu, sigma):
        tmp1 = (x - mu)
        tmp2 = tmp1**2
        tmp3 = tmp1**3
        tmp4 = tmp1**4
        print(0, " : ", tmp1.mean())
        print(sigma**2, " : ", tmp2.mean())
        print(0, " : ", tmp3.mean())
        print(3 * sigma**4, " : ", tmp4.mean(), "\n")
    """
    END
    """

    """
    GBOPTIMIZE
    """
    def derivative(self, axis, atom_i):
        """
        atom_i - номер атома, который будем толкать
        """
        Uq = self.calculate_energy()
        self.atoms[atom_i][axis] += self.delta

        Uq_plus_dq = self.calculate_energy()
        self.atoms[atom_i][axis] += -self.delta

        return (Uq_plus_dq - Uq) / self.delta

    def move_atoms(self, ):
        for i in range(lenght):
            self.atoms[i] += np.where(self.forces[i] > 0, -1, 1) * self.delta * abs(self.forces[i])

    def move_atoms_special(self, blacklist):
        for i in range(lenght):
            if i not in blacklist:
                self.atoms[i] += np.where(self.forces[i] > 0, -1, 1) * self.delta * abs(self.forces[i])

    def relaxate(self, verbose=0):
        """
        Определение оптимальной структуры молекулы
        методом градиентного спуска
        """
        self.find_forces()
        max_force = abs(self.forces).max()

        if verbose != 0:
            print("\nCalculating optimal displacment...")
            while max_force >= self.precision:

                self.find_forces()
                max_force = abs(self.forces).max()
                print(f"energy : {self.calculate_energy():.4f} \t max_force : {max_force:.4f}")
                self.move_atoms()

        else:
            while max_force >= self.precision:

                self.find_forces()
                max_force = abs(self.forces).max()
                self.move_atoms()

    def print_atoms(self, ):
        """
        Метод для распечатки атомов
        """
        print("Atoms:")
        for i in range(lenght):
            self.atoms[i].print_point()

    def atoms_to_file(self, filename):
        with open(filename, "w") as file:
            file.write(str(lenght) + '\n')
            file.write('Si' + '\t' + str(self.atom_radiuses['Si']) + '\n')
            for atom in self.atoms:
                file.write('Si' + '\t' + str(atom.x) + '\t' + str(atom.y) + '\t' + str(atom.z) + '\n')
        file.close()

    def grad(self, atom_i):
        """
        Возвращает градиент функции
        """
        return self.derivative(0, atom_i), self.derivative(1, atom_i), self.derivative(2, atom_i)
    """
    END
    """

    """
    VERLEMETHOD
    """
    def find_velocities(self, T):
        sigma = np.sqrt(8.6 * 10e-5 * 100 * T / self.mass)
        tmp = list()
        for i in range(lenght):
            velocitiy = Point(self.get_rand_num_arr(0, sigma, 1),
                              self.get_rand_num_arr(0, sigma, 1),
                              self.get_rand_num_arr(0, sigma, 1))
            tmp.append(velocitiy)

        self.velocities = tmp

    def print_velocities(self, ):
        """
        Метод для распечатки атомов
        """
        print("Velocities:")
        for i in range(len(self.velocities)):
            self.velocities[i].print_point()

    def find_forces(self, ):
        """
        """
        for i in range(lenght):
            # print(f'fx: {fx:.2f}, fy: {fy:.2f}, fz: {fz:.2f}')
            self.forces[i] = np.array([self.derivative(0, i), self.derivative(1, i), self.derivative(2, i)])

    def find_forces_special(self, blacklist):
        """
        """
        for i in range(lenght):
            if i not in blacklist:
                self.forces[i] = np.array([self.derivative(0, i), self.derivative(1, i), self.derivative(2, i)])

    def find_max_abs_force(self, ):
        """
        """
        force = abs(self.forces)
        return force.max()

    def find_max_abs_force_special(self, blacklist):
        """
        """
        force = np.delete(abs(self.forces), blacklist, axis=0)
        return force.max()

    def find_volume(self, T, dtime=1, num_iters=1000, error=True):
        """
        Метод для придачи движения атомам
        """
        print(f"Calculating volume(T={T})")
        self.find_velocities(T)
        frequency = num_iters / 5
        volumes = np.zeros(int(num_iters / frequency))
        self.find_forces()
        vforces = self.forces
        for i in tqdm(range(num_iters)):
            for atom, force, velocity in zip(self.atoms, self.forces, self.velocities):
                atom.x += velocity.x * dtime + 0.5 * force.x * dtime**2 / self.mass
                atom.y += velocity.y * dtime + 0.5 * force.y * dtime**2 / self.mass
                atom.z += velocity.z * dtime + 0.5 * force.z * dtime**2 / self.mass
            vforces = self.forces
            self.find_forces()
            for vf, force, velocity in zip(vforces, self.forces, self.velocities):
                velocity.x += 0.5 * dtime * (vf.x + force.x) / self.mass
                velocity.y += 0.5 * dtime * (vf.y + force.y) / self.mass
                velocity.z += 0.5 * dtime * (vf.z + force.z) / self.mass

            if i % frequency == 0:
                volumes[int(i / frequency)], _ = self.calculate_volume(precision=1000,
                                                                       n_experiments=1000,
                                                                       alpha=0.95,
                                                                       verbose=0)
        if error is True:
            return volumes.mean(), (volumes.max() - volumes.min()) / 2
        else:
            return volumes.mean()
    """
    END
    """

    """
    Изучение прочности
    """
    # def findlx(self, ):
    #     lst = list()
    #     for atom in self.atoms:
    #         lst.append(atom.x)
    #     return max(lst) - min(lst)
    #
    # def findly(self, ):
    #     lst = list()
    #     for atom in self.atoms:
    #         lst.append(atom.y)
    #     return max(lst) - min(lst)
    #
    # def findlz(self, ):
    #     lst = list()
    #     for atom in self.atoms:
    #         lst.append(atom.z)
    #     return max(lst) - min(lst)

    def relaxate_special(self, blacklist):
        """
        """
        max_force = self.find_max_abs_force_special(blacklist)

        while max_force >= self.precision:

            self.find_forces()
            max_force = self.find_max_abs_force_special(blacklist)
            # print(f"energy : {self.calculate_energy():.4f} \t max_force : {max_force:.4f}")
            self.move_atoms_special(blacklist)

    def pull(self, epsilon):
        for idx in range(lenght):
            if idx != 2 and idx != 4:  # 3 5
                self.atoms[idx][0] *= (1 + epsilon)


def tensile_one_eps(molecule):
    molecule.relaxate_special([2, 4, 8, 9])  # 3 5 9 10

    f1x = molecule.derivative(0, 2)  # 3
    f2x = molecule.derivative(0, 4)  # 5
    f3x = molecule.derivative(0, 8)  # 9
    f4x = molecule.derivative(0, 9)  # 10

    sigma = (abs(f1x) + abs(f2x) + abs(f3x) + abs(f4x)) / 40
    return sigma


def tensile_multiprocessing(molecule, epsilons):  # sigmas
    molecule.relaxate()
    molecule0 = molecule
    sigmas = list()
    molecules = list()

    for epsilon in epsilons:
        epsilon = epsilon * 3 / 100
        molecule.pull(epsilon)
        molecules.append(molecule)
        molecule = molecule0

    with Pool(10) as p:
        sigmas = p.map(tensile_one_eps, molecules)

    return (epsilons * 3 / 100), sigmas

    """
    Изучение прочности
    """

    """
    Изучение термосопротивления
    """
    def thermal_test(self, T):
        """
        Определить термическую устойчивость своей структуры, то есть Тмакс при которой, она может прожить 10пикосекунд, с точностью то 100 градусов К
        1) Отрелаксировать
        2) Нагреть до T хз и подождать 10 пикосекунд
        3) Отрелаксировать
        4) Оптимальная энергия
        5) Если Е0=Е1
        5) Значит не распалась
        6) Значит Т выше! (0 оценка - Тплав кремния, при которой точно распадется и понижать Т, делением поплоам например при 300К не распадется)
        7) Ищем Т макс
        """
        self.find_velocities(T)
        self.find_forces()
        dtime = 1
        f = np.zeros([10, 3])
        volumes = np.zeros(50)
        for i in tqdm(range(1000)):
            for atom, force, velocity in zip(self.atoms, self.forces, self.velocities):
                atom.x += velocity.x * dtime + 0.5 * force.x * dtime**2 / self.mass
                atom.y += velocity.y * dtime + 0.5 * force.y * dtime**2 / self.mass
                atom.z += velocity.z * dtime + 0.5 * force.z * dtime**2 / self.mass

            vforces = self.forces
            self.find_forces()

            for vf, force, velocity in zip(vforces, self.forces, self.velocities):
                velocity.x += 0.5 * dtime * (vf.x + force.x) / self.mass
                velocity.y += 0.5 * dtime * (vf.y + force.y) / self.mass
                velocity.z += 0.5 * dtime * (vf.z + force.z) / self.mass


def checkEquality(mol_base, mol):
    mol.relaxate()
    if np.abs(mol_base.calculate_energy() - mol.calculate_energy()) < 0.001:
        return "Molecule is ok"
    return "Molecule has been desintegrated"


def thermalResistance(mol, T):
    mol_base = mol
    mol.thermal_test(T)
    print(checkEquality(mol_base, mol))
    mol.print_atoms()


def find_period(mol, dtime=1, num_iters=1000):
    """
    Метод для придачи движения атомам
    """
    print(f"Calculating period")
    atoms0 = np.array([[0.2, 0.0, 0.0], [0.0, 0.0, 0.0]], np.float32)
    t_prev = 0
    mol.find_velocities(300)
    mol.find_forces()
    vforces = mol.forces
    for t in tqdm(range(num_iters)):
        for atom, force, velocity in zip(mol.atoms, mol.forces, mol.velocities):
            atom.x += velocity.x * dtime + 0.5 * force.x * dtime**2 / 2912
            atom.y += velocity.y * dtime + 0.5 * force.y * dtime**2 / 2912
            atom.z += velocity.z * dtime + 0.5 * force.z * dtime**2 / 2912

        vforces = mol.forces
        mol.find_forces()
        for vf, force, velocity in zip(vforces, mol.forces, mol.velocities):
            velocity.x += 0.5 * dtime * (vf.x + force.x) / 2912
            velocity.y += 0.5 * dtime * (vf.y + force.y) / 2912
            velocity.z += 0.5 * dtime * (vf.z + force.z) / 2912

        if np.abs(0.2 - mol.atoms[0].x) < 0.1 and np.abs(0.0 - mol.atoms[1].x) < 0.1:
            print(t, t - t_prev)
            t_prev = t
        # if t < 2:
        #     print("0")
        #     mol.print_atoms()

    """
    Изучение термосопротивления
    """

    """
    Изучение теплопроводности
    """
    # Вариант 4


def thermal_conductivity(molecule):

    molecule.relaxate(verbose=1)
    T = 300  # example
    dT = [20, 40, 60, 80, 100]
    time = 100
    molecule.find_velocities(2 * T)  # vx vy=0? vz=0?
    molecule.find_forces()
    dtime = 1

    for i in tqdm(range(100)):

        vforces = copy.deepcopy(molecule.forces)
        for atom, force, velocity in zip(molecule.atoms, molecule.forces, molecule.velocities):
            atom.x += velocity.x * dtime + 0.5 * force.x * dtime**2 / molecule.mass
            atom.y += velocity.y * dtime + 0.5 * force.y * dtime**2 / molecule.mass
            atom.z += velocity.z * dtime + 0.5 * force.z * dtime**2 / molecule.mass

        molecule.find_forces()
        for vf, force, velocity in zip(molecule, molecule.forces, molecule.velocities):
            velocity.x += 0.5 * dtime * (vf.x + force.x) / molecule.mass
            velocity.y += 0.5 * dtime * (vf.y + force.y) / molecule.mass
            velocity.z += 0.5 * dtime * (vf.z + force.z) / molecule.mass

        if (i % 10) == 0:
            """
            Повышаем и понижаем температуру
            """
            pass

    """
    Изучение теплопроводности
    """


def main():
    start_time = datetime.now()
    # /mnt/pool/rhic/1/fkurushin/informatics/Molecule_physics/molecules
    PATH = "/Users/fedorkurusin/Documents/informatics/Molecule_physics/molecules/dataex8.xyz"
    MD = MolecularDynamics(filepath=PATH,
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
                           delta=1e-3,
                           precision=1e-3)

    x, y = tensile_multiprocessing(MD, [i for i in range(10)])
    print(x)
    print(y)
    # plt.plot(x, y)
    # plt.show()

    # f0 = np.array([MD.derivative(i, 0) for i in range(3)])
    # sign = np.where(f0 > 0, -1, 1)
    # print(sign)

    # MD.find_forces()
    # print(MD.forces)
    # MD.relaxate_special([2, 4, 8, 9])
    # print(MD.forces)
    print(f"execution time : {datetime.now() - start_time}")


if __name__ == '__main__':
    main()
