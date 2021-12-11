import math
from point import Point


class Vector(object):
    """docstring for ."""

    def __init__(self, p0, p1):
        """
        x0 - начало
        x1 - цонец
        """
        self.x = p1.x - p0.x
        self.y = p1.y - p0.y
        self.z = p1.z - p0.z

    def lenght(self, ):
        return math.sqrt((self.x)**2 + (self.y)**2 + (self.z)**2)

    def cosine(self, v):
        up = self.x * v.x + self.y * v.y + self.z * v.z
        down = self.lenght() * v.lenght()
        try:
            cos = up / down
            return cos
        except ZeroDivisionError as exc:
            raise ZeroDivisionError('Так как один из вектров - нулевой вектор, то невозможно найти угол между векторами!') from exc


def main():
    p0 = Point(0, 0, 0)
    p1 = Point(1, 1, 0.5)
    p2 = Point(2, 3, 2)
    vector1 = Vector(p0, p1)
    vector2 = Vector(p0, p2)
    len1 = vector1.lenght()
    len2 = vector2.lenght()
    print(vector2.cosine(vector1))
    return 0


if __name__ == '__main__':
    main()
