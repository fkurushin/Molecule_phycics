"""
4. При помощи метода Верле заставить систему двигаться, и найти объем V(T) как среднее значение объема за все время моделирования. T=100, 200 и 300 К; dtime=1 фс:
задать начальные координаты x,y,z и скорости vx,vy,vz (распределение Максвелла)
https://ru.wikipedia.org/wiki/Распределение_Максвелла..
вычислить силы
do time=1, 10000
do i=1,nat
x(i)=x(i)+vx(i)*dtime+0.5d0* (i)*dtime*dtime/massa
y(i)=y(i)+vy(i)*dtime+0.5d0*fy(i)*dtime*dtime/massa
z(i)=z(i)+vz(i)*dtime+0.5d0*fz(i)*dtime*dtime/massa
vfx(i)=fx(i); vfy(i)=fy(i); vfz(i)=fz(i)
enddo
вычислить силы

do i=1,nat
vx(i)=vx(i)+0.5d0*dtime*(vfx(i)+fx(i))/massa
vy(i)=vy(i)+0.5d0*dtime*(vfy(i)+fy(i))/massa
vz(i)=vz(i)+0.5d0*dtime*(vfz(i)+fz(i))/massa
endd
enddo

Upot=0
Ucin=2Ucin(T)

sigma = sqrt(kT / m)
T 0 100 200 300
m = 2912
mu = 0
v0 - из распределения
r0 - отрелаксированная конфигурация

распределение максвелла гаусса P(x)=1/sqrt(2*pi*sigma**2) * exp(-x**2/(2sigma**2))

Построить график V(T). Если зависимость линейна, найти коэффициент теплового расширения.
"""

import numpy as np
from point import Point
from gboptimize import GBOptimize
from randnormal import get_number


class VerleMethod(GBOptimize):
    """docstring for VerleMethod."""

    def __init__(self, atoms):
        FR = FileReader(path)
        self.atoms = atoms
        self.velocity = list()

    def fill_velocity(self, T, mass=2912):
        sigma = np.sqrt(8.6 * 10e-5 * 100 * T / mass)
        tmp = list()
        for i in range(len(self.atoms)):
            velocitiy = Point(get_number(mu=0, sigma=sigma),
                              get_number(mu=0, sigma=sigma),
                              get_number(mu=0, sigma=sigma))
            tmp.append(velocitiy)

        self.velocity = tmp

    def print_atoms(self, ):
        """
        Метод для распечатки атомов
        """
        for atom in range(len(self.atoms)):
            self.atoms[i].print_point()

    def print_velocity(self, ):
        """
        Метод для распечатки атомов
        """
        for i in range(len(self.velocity)):
            self.velocity[i].print_point()

    def move(self, T, dtime=1 * 10e-15):
        """
        Метод для придачи движения атомам
        """
        gb = gboptimize()
        self.fill_velocity(T)
        for i in range(len(self.atoms)):
            atom.x = atom.x + vx * dtime + 0.5 * fx * dtime**2 / mass
            atom.y = atom.y + vy * dtime + 0.5 * fy * dtime**2 / mass
            atom.z = atom.z + vz * dtime + 0.5 * fz * dtime**2 / mass

        # vfx(i) = fx(i)
        # vfy(i) = fy(i)
        # vfz(i) = fz(i)


def main():
    PATH = "/Users/fedorkurusin/Documents/informatics/Molecule_phycics/molecules/check.xyz"
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
                     delta=0.01,
                     precision=0.1)

    VR = VerleMethod(GBO.atoms)
    VR.fill_velocity(T=100)
    VR.print_velocity()


if __name__ == '__main__':
    main()
