#prep_structure.py
#funkcja tworzaca na podstawie pliku z struktura .xyz dwuwarstwową strukture o zadanym rozmiarze
#prep_structure(nazwa_pliku_z_struktura.xyz, rozmiar(np. 3x3))

from ase import Atoms
from ase import io


def prep_struct(structure_file_name, size):
    #wymiary struktury
    n = int(size[0])
    m = int(size[2])

    structure = Atoms(io.read(structure_file_name))
    #ustawienie nowego wymiaru komorki w kierunku osi z
    tmp_cell = structure.get_cell()
    tmp_cell[2] *= 1.4
    structure.set_cell(tmp_cell)
    #tymczasowa struktura potrzebna do dodawania atomow
    tmp_structure = structure.copy()

    #dodanie atomow drugiej warstwy
    for atom in tmp_structure:
        atom_to_add = atom
        atom_to_add.position[2] = tmp_structure.cell[2, 2] - atom.position[2]
        structure.append(atom_to_add)

    return structure * [n, m, 1]
