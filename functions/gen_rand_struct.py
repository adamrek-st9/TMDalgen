#gen_rand_struct.py
#funkcja generujaca na podstawie pliku .xyz  dwuwarstwową strukture o zadanym rozmiarze z losowo umieszczonymi atomami wskazanego pierwiastka
#gen_rand_struct(nazwa_struktury, rozmiar(np. 3x3), symbol_atomu(np. Mo), ilosc_atomow)

from ase import Atoms, Atom
import numpy as np
from functions.small_functions import get_r, calc_distance, get_image_positions, get_rand_xyz
from functions.prep_struct import prep_struct


def gen_rand_struct(structure_file_name, size, atom_symbol, n_atoms):
    #wyjsciowa struktura do ktorej zostana dodane atomy
    structure = prep_struct(structure_file_name, size)

    #wymiary komorki elementarnej
    cell = structure.get_cell()
    array = cell.cellpar()
    a = array[0] #wymiar w kierunku osi x
    b = array[1]
    b_y = b * np.sqrt(3) / 2 #skladowa wektora b w kierunku y
    b_x = b / 2 #skladowa wektora b w kierunku x

    #tymczasowa struktura przechowujaca atomy ktore zostana dodane do koncowej struktury
    tmp_structure = Atoms()
    tmp_structure.set_cell(cell)

    #tymczasowa struktura 3x3 sluzaca do sprawdzenia czy kolejne dodawane atomy nie koliduja z już dodanymi (uwzgledniajac symetrie translacyjna)
    image_structure = Atoms()
    image_structure.set_cell(cell)
    image_structure = image_structure * [3, 3, 1]

    tol_r = 0.1 #tolerancja uzywana przy sprawdzaniu odleglosci miedzy dodawanymi atomami

    #dodanie atomow w zadanej ilosci
    for i in range(n_atoms):
        while True:
            #atom z losowymi wspolrzednymi, ktory chcemy dodac
            random_positions = get_rand_xyz(cell)
            random_atom = Atom(atom_symbol, random_positions)

            #atom z losowymi wspolrzednymi przesunietymi o wektor a i b
            #srodkowy atom w tymczasowej strukturze 3x3
            middle_positions = random_positions + np.array([a - b_x, b_y, 0.])
            middle_atom = Atom(atom_symbol, middle_positions)

            #pozycje na jakie zostana dodane atomy w tymczasowej strukturze 3x3 (symetria translacyjna)
            image_positions = get_image_positions(cell, random_atom)

            collision_with_atom = False #flaga informujaca o kolizji z warstwami
            for atom in structure: #sprawdzenie odleglosci pomiedzy dodawanym atomem a atomami z warstw
                d = calc_distance(random_atom, atom)
                if d < (get_r(random_atom.number) + get_r(atom.number) + tol_r):
                    collision_with_atom = True
                    break

            if not collision_with_atom:
                collision_with_image_atom = False #flaga informujaca o kolizji z obrazami atomow pomiedzy warstwami
                for atom in image_structure: #sprawdzenie czy dodawany atom nie koliduje z obrazami
                    d = calc_distance(middle_atom, atom)
                    if d < (get_r(middle_atom.number) + get_r(atom.number) + tol_r):
                        collision_with_image_atom = True
                        break

                if not collision_with_image_atom:
                    for position in image_positions: #dodanie atomow do struktury z obrazami na wszystkie mozliwe pozycje
                        image_structure.append(Atom(atom_symbol, position))
                    structure.append(random_atom) #dodanie atomu nie kolidujacego z warstwami i obrazami do struktury docelowej
                    break

    return structure
