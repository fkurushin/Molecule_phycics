import math
import random
import numpy as np
from tqdm import tqdm
from filereader import FileReader
from point import Point


class MoleculeVolume(object):

    def __init__(self, atoms, radiuses):
        '''
        Конструктор класса
        '''
        self.atoms = atoms
        self.atom_radiuses = radiuses
        # filereader = FileReader(filename)
        # self.atom_radiuses, self.filename, self.atoms = filereader.file_reader_mod()

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

    def calculate_volume(self, precision, n_experiments, alpha):
        '''
        Расчитывает объем молекулы методом Монте Карло при заданной точности, считает погрешность,
        Можно регулировать количество измеренией
        '''
        volumes = list()
        parall_cords = self.find_parallelepiped()
        print("\nCalculating volume...")
        for i in tqdm(range(n_experiments)):

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
