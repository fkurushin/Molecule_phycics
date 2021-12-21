"""
4. При помощи метода Верле заставить систему двигаться, и найти объем V(T) как среднее значение объема за все время моделирования. T=100, 200 и 300 К; dtime=1 фс:
задать начальные координаты x,y,z и скорости vx,vy,vz (распределение Максвелла)
https://ru.wikipedia.org/wiki/Распределение_Максвелла..
вычислить силы
do time=1, 10000
do i=1,nat
x(i)=x(i)+vx(i)*dtime+0.5d0*fx(i)*dtime*dtime/massa
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
from tqdm import tqdm
from point import Point
from gboptimize import GBOptimize
from randnormal import get_number
from filereader import FileReader
from moleculevolume import MoleculeVolume


class VerleMethod(GBOptimize):
    """docstring for VerleMethod."""

    def __init__(self, ):
        self.forces = list()
        self.mass = 2912
        self.velocities = list()

    def find_velocities(self, T):
        sigma = np.sqrt(8.6 * 10e-5 * 100 * T / self.mass)
        tmp = list()
        for i in range(len(self.atoms)):
            velocitiy = Point(get_number(mu=0, sigma=sigma),
                              get_number(mu=0, sigma=sigma),
                              get_number(mu=0, sigma=sigma))
            tmp.append(velocitiy)

        self.velocities = tmp

    def print_atoms(self, ):
        """
        Метод для распечатки атомов
        """
        print("Atoms:")
        for atom in range(len(self.atoms)):
            self.atoms[i].print_point()

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
        forces = list()
        for i in range(len(self.atoms)):
            fx, fy, fz = self.grad(i)
            force = Point(fx, fy, fz)
            forces.append(force)

        self.forces = forces

    def find_volume(self, T, dtime=1):
        """
        Метод для придачи движения атомам
        """
        MV = MoleculeVolume(self.atoms, self.radiuses)
        self.find_velocities(T)
        volumes = list()
        self.forces = self.find_forces()
        print(type(self.forces))
        vforces = self.forces
        print(vforces[0])
        print("Calculate volume(T)")
        for i in tqdm(range(100)):
            for atom, force, velocity in zip(self.atoms, self.forces, self.velocities):
                atom.x += velocity.x * dtime + 0.5 * force.x * dtime**2 / self.mass
                atom.y += velocity.y * dtime + 0.5 * force.y * dtime**2 / self.mass
                atom.z += velocity.z * dtime + 0.5 * force.z * dtime**2 / self.mass
            vforces = self.forces.copy()
            self.forces = self.find_forces()
            for vf, force, velocity in zip(vforces, self.forces, self.velocities):
                velocity.x += 0.5 * dtime * (vf.x + force.x) / self.mass
                velocity.y += 0.5 * dtime * (vf.y + force.y) / self.mass
                velocity.z += 0.5 * dtime * (vf.z + force.z) / self.mass

            if i % 200 == 0:
                MV.add_new(self.atoms, self.radiuses)
                V, _ = MV.calculate_volume(precision=1000,
                                           n_experiments=1000,
                                           alpha=0.95, verbose=0)
                volumes.append(V)


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

    GBO.find_min()
    print(GBO.atoms)
    # VM = VerleMethod()
    # V1 = VM.find_volume(100)
    # V2 = VM.find_volume(200)
    # V3 = VM.find_volume(300)
    print(V1, V2, V3)


if __name__ == '__main__':
    main()
