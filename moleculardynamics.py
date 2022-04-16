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
from functools import partial
from vector import Vector
from filereader import FileReader
import matplotlib.pyplot as plt
import scipy.integrate as integrate

from datetime import datetime
import time


class MolecularDynamics(object):
    """docstring for MolecularDynamics."""

    def __init__(self, filepath, D0, r0, S, Betta, gamma, c, d, h, two_mu, R, D, delta, precision, mass):
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

        filereader = FileReader(filepath)
        _, self.atoms, _ = filereader.file_reader()

        self.delta = delta
        self.precision = precision
        self.m = mass
        self.mass = self.m * len(self.atoms)
        self.forces = np.zeros((len(self.atoms), 3))
        self.signs = np.zeros((len(self.atoms), 3))
        self.velocities = np.zeros((len(self.atoms), 3))

        self.lenght = len(self.atoms)

        global axises
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
        for k in range(self.lenght):
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
        for i in range(self.lenght):
            for j in range(self.lenght):
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
    def find_parallelepiped(self):
        '''
        Находит координаты параллелепипеда НАДО ПЕРЕПИСАТЬ
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
        Считает объем параллелепипеда НАДО ПЕРЕПИСАТЬ
        '''
        v = (dict['x'][1] - dict['x'][0]) * (dict['y'][1] - dict['y'][0]) * (dict['z'][1] - dict['z'][0])
        return v

    def random_point_in_paral(self, dict):
        '''
        Создает рандомную точку в параллелепипеде НАДО ПЕРЕПИСАТЬ
        '''
        x = random.uniform(dict['x'][0], dict['x'][1])
        y = random.uniform(dict['y'][0], dict['y'][1])
        z = random.uniform(dict['z'][0], dict['z'][1])
        return Point(x, y, z)

    def calculate_volume(self, precision, n_experiments, alpha, verbose=0):
        '''
        Расчитывает объем молекулы методом Монте Карло при заданной точности, считает погрешность,
        Можно регулировать количество измеренией НАДО ПЕРЕПИСАТЬ
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
        for i in range(self.lenght):
            self.atoms[i] += np.where(self.forces[i] > 0, -1, 1) * self.delta * abs(self.forces[i])

    def move_atoms_special(self, blacklist):
        for i in range(self.lenght):
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
        for i in range(self.lenght):
            self.atoms[i].print_point()

    def atoms_to_file(self, filename):
        with open(filename, "w") as file:
            file.write(str(self.lenght) + '\n')
            file.write('Si' + '\t' + str(self.atom_radiuses['Si']) + '\n')
            for atom in self.atoms:
                file.write('Si' + '\t' + str(atom.x) + '\t' + str(atom.y) + '\t' + str(atom.z) + '\n')
        file.close()

    def grad(self, atom_i):
        """
        Возвращает градиент функции
        """
        return self.derivative(0, atom_i), self.derivative(1, atom_i), self.derivative(2, atom_i)

    # Изменить одну из координат на шаг step для atom
    def move_i_coord(self, step, atom, coord):
        new_coords = self.atoms[atom].copy()
        new_coords[coord] = new_coords[coord] + step
        self.atoms[atom] = new_coords

    """
    END
    """

    """
    VERLEMETHOD
    """
    def find_velocities(self, T):
        # k = 0.00008617
        # sigma = np.sqrt(k * T / 2700)
        # for i in range(self.lenght):
        #     self.velocities[i] = np.random.normal(0, sigma, 3)

        for i in range(self.lenght):
            for j in range(3):
                phi = random.uniform(0, 2 * np.pi)
                r = random.uniform(0.0000001, 1)

                ro = np.sqrt(-2 * np.log(r))
                x = ro * np.cos(6.28 * phi)
                k = 0.00008617
                mu = 0
                sigma = np.sqrt(2 * k * T / self.m)
                self.velocities[i][j] = mu + sigma * x

    def find_forces(self, ):
        """
        """
        for i in range(self.lenght):
            # print(f'fx: {fx:.2f}, fy: {fy:.2f}, fz: {fz:.2f}')
            self.forces[i] = np.array([self.derivative(0, i), self.derivative(1, i), self.derivative(2, i)])

    def find_forces_special(self, blacklist):
        """
        """
        for i in range(self.lenght):
            if i not in blacklist:
                self.forces[i] = np.array([self.derivative(0, i), self.derivative(1, i), self.derivative(2, i)])

    def find_max_abs_force_special(self, blacklist):
        """
        """
        return np.delete(abs(self.forces), blacklist, axis=0).max()

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
    def relaxate_special(self, blacklist):
        """
        """
        max_force = np.delete(abs(self.forces), blacklist, axis=0).max()

        while max_force >= self.precision:

            self.find_forces()
            max_force = np.delete(abs(self.forces), blacklist, axis=0).max()
            print(f"energy : {self.calculate_energy():.4f} \t max_force : {max_force:.4f}")
            self.move_atoms_special(blacklist)

    def pull(self, epsilon):
        for idx in range(self.lenght):
            if idx != 2 and idx != 4:  # 3 5
                self.atoms[idx][0] *= (1 + epsilon)

    def kinetic_energy(self, ):
        return self.m * np.sum(self.velocities**2) / 2


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
    molecule0 = copy.deepcopy(molecule)
    sigmas = list()
    molecules = list()
    eps_new = list()

    for epsilon in epsilons:
        epsilon = epsilon * 4 / 100
        molecule.pull(epsilon)
        molecules.append(molecule)
        molecule = copy.deepcopy(molecule0)
        eps_new.append(epsilon)

    with Pool() as p:
        sigmas = p.map(tensile_one_eps, molecules)

    return eps_new, sigmas

    """
    Изучение прочности
    """


def thermal_test(molecule, num_steps=100):
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

    T = np.array([500 * i for i in range(1, 10)])
    energy = np.zeros(num_steps)
    energy_mean = np.zeros(len(T))

    for i in tqdm(range(len(T))):

        molecule.find_velocities(T[i])
        dtime = 1
        for j in range(num_steps):
            molecule.atoms += molecule.velocities * dtime + 0.5 * molecule.forces * dtime**2 / molecule.mass
            vf = molecule.forces
            molecule.find_forces()
            molecule.velocities += 0.5 * dtime * (vf + molecule.forces) / molecule.mass
            energy[j] = molecule.calculate_energy()

        energy_mean[i] = energy.mean()

    return T, energy_mean


def thermal_conductivity_one_processing(dT, time=100, T=300, relax=False):
    """
    1) 14 атомов 1D цепочка 2А
    2) Отрелаксировать только вдоль оси x
    3) Придать всем атомам случайные скорости temp =  2T
    4) Запустить алгоритм Верле
    5) Каждые 10 шагов умножать скорости первых 3 атомов на K > 1 K = m/2(v1**2 + V2**2 + v3**2), а Kдолжна = 3 * 3/2k(T+T/2)
    Отсюда к-цент => k = sqrt(kдолж / k). k для последних трех атомов Kдолжна = 3 * 3/2k(T-T/2)
    6) time = 100 шагов dtime = 1 T = 300K dT = 20 40 60 80 100
    7) График J от dT, где J = kдолжн - k
    """
    PATH = "/Users/fedorkurusin/Documents/informatics/Molecule_physics/molecules/si14_1d_chain.xyz"
    molecule = MolecularDynamics(filepath=PATH,
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

    if relax is True:
        molecule.relaxate()
    molecule.find_velocities(2 * T)
    k = 8.62 * 10e-5
    dtime = 1

    for t in range(time):
        molecule.atoms += molecule.velocities * dtime + 0.5 * molecule.forces * dtime**2 / molecule.mass
        vf = molecule.forces
        molecule.find_forces()
        molecule.velocities += 0.5 * dtime * (vf + molecule.forces) / molecule.mass

        if t % 10 == 0:
            coefficinet = np.sqrt((3 * 3 / 2 * k * (T + dT / 2)) / (molecule.m / 2 * (molecule.velocities[0]**2 + molecule.velocities[1]**2 + molecule.velocities[2]**2)))
            molecule.velocities[0] *= coefficinet
            molecule.velocities[1] *= coefficinet
            molecule.velocities[2] *= coefficinet

            coefficinet = np.sqrt((3 * 3 / 2 * k * (T - dT / 2)) / (molecule.m / 2 * (molecule.velocities[-1]**2 + molecule.velocities[-2]**2 + molecule.velocities[-3]**2)))
            molecule.velocities[-1] *= coefficinet
            molecule.velocities[-2] *= coefficinet
            molecule.velocities[-3] *= coefficinet

        return 3 * 3 / 2 * k * (T + dT / 2) - molecule.m / 2 * (molecule.velocities[0][0]**2 + molecule.velocities[1][0]**2 + molecule.velocities[2][0]**2)


def thermal_conductivity_multiprocessing():
    """
    Используя этот метод, я не могу инициализировать класс в main функции
    """
    with Pool(5) as p:
        J = p.map(thermal_conductivity_one_processing, [20, 40, 60, 80, 100])

    return [20, 40, 60, 80, 100], J


def thermal_conductivity(molecule, dT_list, time=100, T=300, relax=False):
    """
    1) 14 атомов 1D цепочка 2А
    2) Отрелаксировать только вдоль оси x
    3) Придать всем атомам случайные скорости temp =  2T
    4) Запустить алгоритм Верле
    5) Каждые 10 шагов умножать скорости первых 3 атомов на K > 1 K = m/2(v1**2 + V2**2 + v3**2), а Kдолжна = 3 * 3/2k(T+T/2)
    Отсюда к-цент => k = sqrt(kдолж / k). k для последних трех атомов Kдолжна = 3 * 3/2k(T-T/2)
    6) time = 100 шагов dtime = 1 T = 300K dT = 20 40 60 80 100
    7) График J от dT, где J = kдолжн - k
    """
    J = list()
    if relax is True:
        molecule.relaxate()
    molecule.find_velocities(2 * T)
    k = 8.62 * 10e-5
    dtime = 1

    for dT in tqdm(dT_list):
        for t in range(time):
            molecule.atoms += molecule.velocities * dtime + 0.5 * molecule.forces * dtime**2 / molecule.mass
            vf = molecule.forces
            molecule.find_forces()
            molecule.velocities += 0.5 * dtime * (vf + molecule.forces) / molecule.mass

            if t % 10 == 0:
                coefficinet = np.sqrt((3 * 3 / 2 * k * (T + dT / 2)) / (molecule.m / 2 * (molecule.velocities[0]**2 + molecule.velocities[1]**2 + molecule.velocities[2]**2)))
                molecule.velocities[0] *= coefficinet
                molecule.velocities[1] *= coefficinet
                molecule.velocities[2] *= coefficinet

                coefficinet = np.sqrt((3 * 3 / 2 * k * (T - dT / 2)) / (molecule.m / 2 * (molecule.velocities[-1]**2 + molecule.velocities[-2]**2 + molecule.velocities[-3]**2)))
                molecule.velocities[-1] *= coefficinet
                molecule.velocities[-2] *= coefficinet
                molecule.velocities[-3] *= coefficinet

        J.append(3 * 3 / 2 * k * (T + dT / 2) - molecule.m / 2 * (molecule.velocities[0][0]**2 + molecule.velocities[1][0]**2 + molecule.velocities[2][0]**2))

    return J

    """
    Изучение термосопротивления
    """

    """
    Калориметрическая кривая
    """


def calorimetry_curve(molecule, T):
    time = 10
    dtime = 1
    k = 1 / 11602

    molecule.find_velocities(T)
    kinetic_energy = molecule.kinetic_energy()
    temperature = np.zeros(time)
    for t in range(time):
        molecule.atoms += molecule.velocities * dtime + 0.5 * molecule.forces * dtime**2 / molecule.m
        vf = molecule.forces
        molecule.find_forces()
        molecule.velocities += 0.5 * dtime * (vf + molecule.forces) / molecule.m
        temperature[t] = 2 / 3 * 1 / (molecule.lenght * k) * molecule.kinetic_energy()

    return kinetic_energy, np.mean(temperature)


def calorimetry_curve_multiprocessing(molecule, T=[0, 200, 400, 600, 800, 1000], relax=False):
    if relax is True:
        molecule.relaxate()

    with Pool(6) as p:
        Q_temperature = p.map(partial(calorimetry_curve, molecule), T)

    return Q_temperature

    """
    Калориметрическая кривая
    """

    """
    Спектр молекулы
    """


def second_deriv_E(mol, i, j):
    atom_i, coord_i = (i // 3, i % 3)
    atom_j, coord_j = (j // 3, j % 3)
    mol_move = copy.copy(mol)
    mol_move.move_i_coord(mol.delta, atom_i, coord_i)
    mol_move.move_i_coord(-mol.delta, atom_j, coord_j)
    E_ul = mol_move.calculate_energy()
    mol_move.move_i_coord(-2 * mol.delta, atom_i, coord_i)
    E_dl = mol_move.calculate_energy()
    mol_move.move_i_coord(2 * mol.delta, atom_j, coord_j)
    E_dr = mol_move.calculate_energy()
    mol_move.move_i_coord(2 * mol.delta, atom_i, coord_i)
    E_ur = mol_move.calculate_energy()
    second_deriv = ((E_ur - E_ul) - (E_dr - E_dl)) / (4 * mol.delta ** 2)
    return second_deriv


def second_deriv_E_matrix(mol):
    sd_matrix = np.array([[second_deriv_E(mol, i, j) for i in range(3 * mol.lenght)] for j in range(3 * mol.lenght)])
    return sd_matrix


def spectrum_vals(omegas, min_val, max_val, dots):
    sigma = 10
    sum_of_gausses = 0
    spectrum_dots = np.empty(dots)
    omega_dots = np.arange(min_val, max_val, (max_val - min_val) / dots)
    i = 0
    for omega in omega_dots:
        for omega_k in omegas:
            sum_of_gausses += np.exp(-(omega - omega_k) ** 2 / (2 * sigma ** 2))
        spectrum_dots[i] = sum_of_gausses
        sum_of_gausses = 0
        i += 1
    return (spectrum_dots, omega_dots)


def spectrum_plot(mol, dots):
    sd_matrix = second_deriv_E_matrix(mol)
    eigen_values = np.linalg.eig(sd_matrix)[0]
    print(eigen_values)
    omegas = np.array([(1 / mol.m) * abs(eigen_values[i]) ** 0.5 * (10 ** 15) / (3 * 10 ** 9) for i in range(eigen_values.size) if np.abs(eigen_values[i]) > 0.01])
    min_val = omegas.min() * 0
    max_val = omegas.max() * 1.1
    spectrum_dots, omega_dots = spectrum_vals(omegas, min_val, max_val, dots)
    fig, ax = plt.subplots()
    ax.plot(omega_dots, spectrum_dots)
    fig.set_figwidth(12)
    fig.set_figheight(7)
    ax.set_title('Спектр')
    plt.show()

    """
    Спектр молекулы
    """

    """
    Автокорреляция
    """


def verle_vel(molecule, T):
    # Алгоритм Верле, надо взять димер C2 ~ 300 одномерный,
    # надо взять еще формулу и там по 6 атомам сложить
    time = 1000
    dtime = 1
    dim = 3
    v = np.zeros((time, molecule.lenght, dim))
    molecule.find_velocities(T)
    for t in range(time):
        molecule.atoms += molecule.velocities * dtime + 0.5 * molecule.forces * dtime**2 / molecule.m
        vf = molecule.forces
        molecule.find_forces()
        molecule.velocities += 0.5 * dtime * (vf + molecule.forces) / molecule.m
        v[t] = molecule.velocities

    return v.reshape(molecule.lenght, dim, time)


def autocorrelation(molecule, tau):
    vel_lenght_dim = verle_vel(molecule)
    coefficinets = np.zeros(molecule.lenght, dim)
    time = 1000
    for i in range(len(vel_lenght_dim)):
        for j in range(len(vel_dim)):
            time = len(vel_lenght_dim[i][j])
            vel1 = vel_lenght_dim[i][j][0:time - tau]
            vel2 = vel_lenght_dim[i][j][tau:]
            vel_mean = vel.mean()
            vel_var = np.array([i**2 for i in vel - vel_mean]).sum()
            auto_corr = 0
            for i in range(time - tau):
                temp = (vel1[i] - vel_mean) * (vel2[i] - vel_mean) / vel_var
                auto_corr = auto_corr + temp
            coefficinets[i][j] = auto_corr

    return auto_corr


def P(w):
    return integrate.quad(lambda t: np.cos(w * t), 0, np.inf)[0]


def spectrum_plot_autocorr(molecule, T):
    v = verle_vel(molecule, T)
    v0tv00 = v[0][0] * v[0][0][0]
    v1tv10 = v[1][0] * v[1][0][0]
    mean_v = np.mean(v0tv00 + v1tv10)
    omega = np.array([w for w in range(1000)])
    # min_val = omegas.min() * 0
    # max_val = omegas.max() * 1.1
    p = np.array([mean_v * P(w) for w in range(1000)])

    # spectrum_dots, omega_dots = spectrum_vals(omegas, min_val, max_val, dots)

    plt.plot(omega, p)
    # plt.plot(p)
    plt.show()


def spectrum_plot1(dots):
    eigen_values = np.array([0.0232854, 0.0116322, 0.0122042, 0.0101042, 0.0093638, 0.0085536, 0.0098193, 0.0006366, 0.0035970, 0.0023074, 0.0015176, 0.0012631, 0.0003403, -0.0002513, -0.0002415, 0.0001060, 0.0000742, 0.0000596, -0.0000000, -0.0000000, 0.0000000, 0.0000000, -0.0000000, 0.0000000])
    omegas = np.array([abs(eigen_values[i]) ** 0.5 * (10 ** 15) / (3 * 10 ** 10) for i in range(eigen_values.size) if np.abs(eigen_values[i]) > 0.00001])
    min_val = omegas.min() * 0
    max_val = omegas.max() * 1.1
    spectrum_dots, omega_dots = spectrum_vals(omegas, min_val, max_val, dots)
    # print("lol")
    print()
    print(omegas)
    # print(spectrum_dots)
    fig, ax = plt.subplots()
    ax.plot(omega_dots, spectrum_dots)
    fig.set_figwidth(12)
    fig.set_figheight(7)
    ax.set_title('Спектр')
    plt.show()

    """
    Автокорреляция
    """

    """
    Енергия активации
    """


def two_mols_research(mol57, mol66):
    mol = copy.deepcopy(mol57)
    energy = np.arange(0, 1, 0.01)
    alpha = np.arange(0, 1, 0.01)
    for i, a in enumerate(alpha):
        mol.atoms = a * mol57.atoms + (1 - a) * mol66.atoms
        energy[i] = mol.calculate_energy()
    Ea = energy.max() - mol57.calculate_energy()
    Et = mol66.calculate_energy() - mol57.calculate_energy()
    return energy, alpha, Ea, Et

    """
    Енергия активации
    """


def main():
    start_time = datetime.now()
    PATH = "/Users/fedorkurusin/Documents/informatics/Molecule_physics/molecules/C2.xyz"
    molecule = MolecularDynamics(filepath=PATH,
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
                                 precision=1e-2,
                                 mass=300)

    spectrum_plot1(100000)
    print(f"execution time : {datetime.now() - start_time}")


if __name__ == '__main__':
    main()
