from gboptimize import GBOptimize
from moleculevolume import MoleculeVolume
"""
Удалить файл atomicenergy
"""
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

print("\nBefore:")
GBO.print_atoms()
MV = MoleculeVolume(GBO.atoms, GBO.atom_radiuses)  # Сделать так, чтобы не вызывать этот конструктор 100 раз
print("V : ", MV.calculate_volume(precision=1000,
                                  n_experiments=1000,
                                  alpha=0.95))
print("E :", GBO.calculate_energy())

GBO.relaxate()  # Ищем оптимальную структуру

print("\nAfter:")

GBO.print_atoms()
MV.add_new(GBO.atoms, GBO.atom_radiuses)
print("V : ", MV.calculate_volume(precision=1000,
                                  n_experiments=1000,
                                  alpha=0.95))
print("E :", GBO.calculate_energy())
