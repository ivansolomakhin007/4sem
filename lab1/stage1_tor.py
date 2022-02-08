import gmsh
import sys
from math import pi, cos, sin

gmsh.initialize()
gmsh.option.setNumber("Mesh.SaveAll", 1)
gmsh.option.setNumber("Mesh.Algorithm", 6)

gmsh.model.add("tor")

lc = 5

R1 = 2  # радиус основания
R2 = 1
h = 5  # расстояния от центра тора до центра окружности сечения тора


def tor(R, h, n, shift, orientation, tag):
    """Строит тороидальную поверхность
    R - радиус сечения тора
    n - количество углов сечения
    h - расстояние от центра тора до центра сечения
    tag - тег поверхности, принимает значение 0 или 1
    shift - сдвиг тегом точек/кривых для удобства пользования функцией и исключения пересечения тегов"""
    for i in range(n):
        for j in range(n):
            alpha = 0.5 * (2 * pi / n) + j * (2 * pi / n)
            beta = 0.5 * (2 * pi / n) + i * (2 * pi / n)
            gmsh.model.geo.addPoint(h * cos(beta) - R * cos(alpha) * cos(beta),
                                    h * sin(beta) - R * cos(alpha) * sin(beta), R * sin(alpha), lc, i * n + j + shift)
    for i in range(n):
        for j in range(n):
            gmsh.model.geo.addLine(i * n + j + shift, i * n + (j + 1) % n + shift, 4 * (i * n + j) + 1 + shift)
            gmsh.model.geo.addLine(i * n + (j + 1) % n + shift, ((i + 1) % n) * n + (j + 1) % n + shift, 4 * (i * n + j) + 2 + shift)
            gmsh.model.geo.addLine(((i + 1) % n) * n + (j + 1) % n + shift, ((i + 1) % n) * n + j % n + shift,
                                   4 * (i * n + j) + 3 + shift)
            gmsh.model.geo.addLine(((i + 1) % n) * n + j % n + shift, i * n + j + shift, 4 * (i * n + j) + 4 + shift)
            gmsh.model.geo.addCurveLoop(
                [4 * (i * n + j) + 1 + shift, 4 * (i * n + j) + 2 + shift, 4 * (i * n + j) + 3 + shift, 4 * (i * n + j) + 4 + shift], i * n + j + shift + 1)
            gmsh.model.geo.addPlaneSurface([i * n + j + shift + 1], i * n + j + shift)
    l = gmsh.model.geo.addSurfaceLoop([i * n + j + shift for i in range(n) for j in range(n)])
    gmsh.model.geo.addVolume([orientation*l])


tor(R1, h, 20, 0, 1, 1)
tor(R2, h, 20, 10000, -1, -1)


gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

gmsh.write("tor.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
