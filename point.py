import random
import warnings
import numpy as np


class Point(object):
    '''
    Класс, создающий точку или атом. Для того, чтобы создать атом нужно присвоить ему имя название,
    как в периодической таблице химических элементов Менделеева и(если нужно) радиус атома
    '''

    def __init__(self, x=0, y=0, z=0, radius=np.NaN, name=np.NaN):
        '''
        Конструктор, по умолчанию, создающий точку(К сожалению
        полиморфизма констуркторов в python нет)
        '''

        self.x = x
        self.y = y
        self.z = z
        self.R = radius
        self.name = name

    def add(self, p):
        '''
        Складывает две точки
        '''
        x = self.x + p.x
        y = self.y + p.y
        z = self.z + p.z

    def min(self, p):
        '''
        Вычитает две точки
        '''
        self.x = self.x - p.x
        self.y = self.y - p.y
        self.z = self.z - p.z

    def print_point(self):
        '''
        Выводит на экран коодинаты точки или атома
        '''
        print("(%.10f, %.10f, %.10f)" % (self.x, self.y, self.z))

    def circle_equation(self, p):
        '''
        Подстановка в уравнение сферы, проверяемой точки
        '''
        return (self.x - p.x)**2 + (self.y - p.y)**2 + (self.z - p.z)**2

    def random_point_nearby(self, atom, precision=2):
        '''
        Создает рандомную точку в кубе описанном около сферы. То есть не все эти точки нам нужны
        Передаем информацию об атоме:
        '''
        self.x = round(random.uniform(atom.x - atom.R, atom.x + atom.R), precision)
        self.y = round(random.uniform(atom.y - atom.R, atom.y + atom.R), precision)
        self.z = round(random.uniform(atom.z - atom.R, atom.z + atom.R), precision)

    def is_in_atom(self, p):
        '''
        Метод проверяющий попадает ли заданная точка в окружность или нет.
        '''
        if self.circle_equation(p) <= p.R**2:
            return True
        else:
            return False

    def is_it_crossed(self, list_of_p):
        '''
        Принимает список, соседних атомов и узнает попадает ли точка в объемы этих атомов
        '''
        for p in list_of_p:
            if self.circle_equation(p) <= p.R**2:
                return True
                break
        return False

    def move(self, axis, delta):
        if axis == 'x':
            self.x = self.x + delta
        elif axis == 'y':
            self.y = self.y + delta
        elif axis == 'z':
            self.z = self.z + delta
        else:
            warnings.warn("This axis doesnt exist")


def main():
    p = Point()
    p.move(axis='x', delta=0.1)
    p.circle_equation(Point(0, 0, 1))


if __name__ == '__main__':
    main()
