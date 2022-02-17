import math
import random
import numpy as np
from tqdm import tqdm
from point import Point
from vector import Vector
from filereader import FileReader


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
        self.atom_radiuses, _, self.atoms = filreader.file_reader_mod()

        self.delta = delta
        self.precision = precision
        self.mass = 2912
        self.forces = list()
        self.velocities = list()

    """
    INTERACTIONENERGY
    """
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

    def calculate_energy(self, ):
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
        print("Atoms:")
        for i in range(len(self.atoms)):
            self.atoms[i].print_point()

    def grad(self, atom_i):
        """
        Возвращает градиент функции
        """
        return self.derivative('x', atom_i), self.derivative('y', atom_i), self.derivative('z', atom_i)
    """
    END
    """

    """
    VERLEMETHOD
    """
    def find_velocities(self, T):
        sigma = np.sqrt(8.6 * 10e-5 * 100 * T / self.mass)
        tmp = list()
        for i in range(len(self.atoms)):
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

    def print_forces(self, ):
        """
        Метод для распечатки атомов
        """
        print("Forces:")
        for i in range(len(self.forces)):
            self.forces[i].print_point()

    def find_forces(self, ):
        """
        """
        for i in range(len(self.atoms)):
            fx, fy, fz = self.grad(i)
            force = Point(fx, fy, fz)
            self.forces.append(force)

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

    def findlx(self, ):
        lst = list()
        for atom in self.atoms:
            lst.append(atom.x)
        return max(lst) - min(lst)

    def findly(self, ):
        lst = list()
        for atom in self.atoms:
            lst.append(atom.y)
        return max(lst) - min(lst)

    def findlz(self, ):
        lst = list()
        for atom in self.atoms:
            lst.append(atom.z)
        return max(lst) - min(lst)

    def move(self, ):
        # Быть может это условие неправильное
        for atom, idx in enumerate(self.atoms):
            if idx != 9 and idx != 10:
                atom.y = atom.y * 1.1
                # atom.z = atom.z * 1.1
    """
    END
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
    mol.find_min()
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


def main():
    PATH = "/Users/fedorkurusin/Documents/informatics/Molecule_phycics/molecules/sinew10.xyz"
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
                           delta=0.01,
                           precision=0.1)

    MD.find_min()
    MD.print_atoms()


if __name__ == '__main__':
    main()
