#crossover.py
#funkcja operatora krzyzowania, tworzaca obiekt klasy Atoms bedacy wynikiem krzyzowania (wymiany zadanej ilosci atomow)
#pomiedzy strukturami rodzicow
#crossover(rodzic1, rodzic2, liczba_wymienianych_atomow, nazwa_pliku_ze_struktura(np. MoS2), rozmiar_struktury(np. 3x3))

from ase import Atoms
import math
import random
from itertools import combinations
from functions.small_functions import get_r, second_part
from functions.prep_struct import prep_struct


def crossover(atoms1, atoms2, n, struct_filename, struct_size):
    #obiekty rodzicow - takie same komorki jak atoms1 i atoms2
    #zostana do nich dodane atomy znajdujace sie miedzy warstwami w atoms1 i atoms2
    parent1 = Atoms()
    parent1.set_cell(atoms1.cell)
    parent2 = Atoms()
    parent2.set_cell(atoms2.cell)

    #wyjsciowe obiekt potomka
    child = prep_struct(struct_filename, struct_size)

    #wymiary komorki elementarnej rodzicow
    cell = atoms1.get_cell()
    array = atoms1.cell.cellpar()
    c = array[2]  #wymiar w kierunku osi z
    c_half = c / 2 #polowa wysokosci komorki

    tol_r = 0.1  #tolerancja uzywana przy sprawdzaniu odleglosci miedzy dodawanymi atomami

    #dodanie atomow znajdujacych sie miedzy warstwami w obiektach atoms1 i atoms2 do obiektow rodzicow
    for atom in atoms1:
        if math.isclose(atom.position[2], c_half, abs_tol=1E-1):
            parent1.append(atom)
    for atom in atoms2:
        if math.isclose(atom.position[2], c_half, abs_tol=1E-1):
            parent2.append(atom)

    combinations1 = list(combinations(parent1, n))
    combinations2 = list(combinations(parent2, n))

    while True:
        tmp_child1 = Atoms()
        tmp_child1.set_cell(cell)
        tmp_child2 = Atoms()
        tmp_child2.set_cell(cell)
        parent1_a = Atoms(random.choice(combinations1))  # czesc a parent1 stworzona z atomow wybranych podczas kombinacji
        parent1_b = second_part(parent1, parent1_a)  # czesc b parent1 stworzona z pozostalych atomow
        parent2_a = Atoms(random.choice(combinations2))  # czesc a parent2 stworzona z atomow wybranych podczas kombinacji
        parent2_b = second_part(parent2, parent2_a)  # czesc b parent2 stworzona z pozostalych atomow

        for atom in parent1_a + parent2_b: #dodanie atomow miedzy warstwami do potomka nr1
            tmp_child1.append(atom)

        for atom in parent1_b + parent2_a: #dodanie atomow miedzy warstwami do potomka nr 2
            tmp_child2.append(atom)

        #sprawdzenie kolizji miedzy atomami i ich obrazami w tmp_child1
        image_structure1 = tmp_child1 * [3, 3, 1]
        collision_in_child1 = False
        for i in range(len(image_structure1)):
            for j in range(i + 1, len(image_structure1)):
                distance = image_structure1.get_distance(i, j)
                if distance < (2 * get_r(image_structure1.get_atomic_numbers()[0]) + tol_r):
                    collision_in_child1 = True
                    break

        #sprawdzenie kolizji miedzy atomami i ich obrazami w tmp_child2
        image_structure2= tmp_child2 * [3, 3, 1]
        collision_in_child2 = False
        for i in range(len(image_structure2)):
            for j in range(i + 1, len(image_structure2)):
                distance = image_structure2.get_distance(i, j)
                if distance < (2 * get_r(image_structure2.get_atomic_numbers()[0]) + tol_r):
                    collision_in_child2 = True
                    break
        #jesli nie ma kolizji w tmp_child1 i tmp_child2 - utworzenie wyjsciowych struktur potomkow
        if not collision_in_child1:
            for atom in tmp_child1:
                child.append(atom)
            break
        if not collision_in_child2:
            for atom in tmp_child2:
                child.append(atom)
            break

    return child
