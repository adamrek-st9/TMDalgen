#small_functions.py
#plik z mniejszymi, pobocznymi funkcjami, ktora sa wykorzystane w innych funkcjach i skryptach

import numpy as np
from ase import Atoms
from ase.data import covalent_radii
import random


#funkcja obliczajaca promien atomu
def get_r(atomic_number):
    tmp_r = 0.9 * covalent_radii[atomic_number]
    return tmp_r

#funkcja obliczajaca odleglosc pomiedzy dwoma atomami
def calc_distance(atom1, atom2):
    tmp_d = np.linalg.norm(atom1.position - atom2.position)
    return tmp_d

#funkcja zwracajaca atomy z obiektu parent, ktorych nie mam w part
#potrzebna w crossover
def second_part(parent, part):
    tmp_second_part = Atoms()
    tmp_second_part.set_cell(parent.cell)
    for parent_atom in parent:
        found = False
        for part_atom in part:
            if np.allclose(parent_atom.position, part_atom.position, atol=1e-4):
                found = True
                break
        if not found:
            tmp_second_part.append(parent_atom)
    return tmp_second_part

#funkcja zwracajace losowe wartosci x,y,z w odpowiednich
#przedzialach zaleznych od rozmiarow komorki a,b,c
def get_rand_xyz(cell):
    array = cell.cellpar()
    a = array[0]  # wymiar w kierunku osi x
    b = array[1]
    c = array[2]  # wymiar w kierunku osi z
    b_y = b * np.sqrt(3) / 2  # skladowa wektora b w kierunku y

    y0 = 0  #zakres dozwolonych wartosci y
    y1 = b_y
    tmp_y = random.uniform(y0, y1)
    x_prim = tmp_y * np.sqrt(3) / 3
    x0 = 0 - x_prim  #zakres dozwolonych wartosci x
    x1 = a - x_prim
    tmp_x = random.uniform(x0, x1)
    tmp_z = c / 2 #wartosc z

    return np.array([tmp_x, tmp_y, tmp_z])

#funkcja zwracajaca pozycje obrazow atomu
def get_image_positions(cell, atom):
    array = cell.cellpar()
    a = array[0]  # wymiar w kierunku osi x
    b = array[1]
    b_y = b * np.sqrt(3) / 2  # skladowa wektora b w kierunku y
    b_x = b / 2  # skladowa wektora b w kierunku x

    tmp_image_positions = [atom.position,
                           atom.position + np.array([a, 0., 0.]),
                           atom.position + np.array([2 * a, 0., 0.]),
                           atom.position + np.array([-b_x, b_y, 0.]),
                           atom.position + np.array([a - b_x, b_y, 0.]),
                           atom.position + np.array([2 * a - b_x, b_y, 0.]),
                           atom.position + np.array([-2 * b_x, 2 * b_y, 0.]),
                           atom.position + np.array([a - 2 * b_x, 2 * b_y, 0.]),
                           atom.position + np.array([2 * a - 2 * b_x, 2 * b_y, 0.])]
    return tmp_image_positions
