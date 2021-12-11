import os
import regex as re
from collections import defaultdict
from point import Point


class FileReader(object):

    def __init__(self, filename):
        '''
        Инифиализирует класс файлридер
        '''
        self.filename = filename
        if os.path.exists(self.filename):
            self.file = open(os.path.abspath(self.filename), 'r')
            # print("\n", self.filename)

    def file_reader(self):
        '''
        Метод читающий файл формата .xyz
        '''
        atoms = list()
        xyz = open(self.filename)
        n_atoms = int(xyz.readline())
        title = xyz.readline()

        for line in xyz:
            name, x, y, z = line.split()
            atom = Point(float(x), float(y), float(z))
            atoms.append(atom)

        xyz.close()

        return self.filename, tuple(atoms), title

    def file_reader_mod(self):
        '''
        Метод читающий файл формата .xyz, который был дан преподавателем
        '''
        atoms = list()
        xyz = open(self.filename)
        n_atoms = int(xyz.readline())
        title = xyz.readline()
        radiuses = re.findall(r'\d*\.\d*', title)
        names = re.findall(r'[A-Za-z]+', title)
        atom_radiuses = dict()

        for name, radius in zip(names, radiuses):
            atom_radiuses.update({name: float(radius)})

        for line in xyz:
            name, x, y, z = line.split()
            atom = Point(float(x), float(y), float(z), float(atom_radiuses[name]), name)
            atoms.append(atom)

        xyz.close()

        return atom_radiuses, self.filename, tuple(atoms)
