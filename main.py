import numpy as np
from mayavi import mlab
from ase.atoms import Atoms, Cell
from ase.io import read
from ase.data import covalent_radii as radii, atomic_numbers as numbers
from ase.data.colors import jmol_colors as colors

from ase.visualize import view

def draw_cell(cell: np.ndarray | Cell, figure=None):
    lines = []
    vertices = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1], \
                        [0,1,1], [1,0,1], [1,1,0], [1,1,1]])
    for i in np.arange(vertices.shape[0]):
        for j in np.arange(i+1, vertices.shape[0]):
            if np.linalg.norm(vertices[i,:]-vertices[j,:]) <= 1:
                line = vertices[[i,j],:] @ cell
                l = mlab.plot3d(line[:,0], line[:,1], line[:,2])
                lines.append(l)
    return lines

def draw_atoms(atoms: Atoms, figure) -> dict:
    points = {}
    for symbol in set(atoms.symbols):
        number = numbers[symbol]
        mask = atoms.symbols == symbol
        p = mlab.points3d(atoms.positions[mask,0], atoms.positions[mask,1], atoms.positions[mask,2],
                      scale_factor=radii[number]/radii[1]/2,
                      color=tuple(colors[number]),
                      resolution=20,
                      figure=figure)
        points[symbol] = p
    return points

def ashow_single(atoms: Atoms, figure=None):
    draw_atoms(atoms)
    draw_cell(atoms.cell, figure)

def ashow_traj(traj: list[Atoms], figure=None):
    symbols = traj[0].symbols
    points = draw_atoms(atoms=traj[0], figure=figure)
    @mlab.animate(delay=50, ui=False)
    def anim():
        for image in traj:
            for symbol in points.keys():
                mask = symbols == symbol
                points[symbol].mlab_source.trait_set(x=image.positions[mask,0], y=image.positions[mask,1], z=image.positions[mask,2])
            yield
    anim()

# atoms = read('test/c2db-6222.xyz')
# ashow_single(atoms)

atoms = read(r'C:\Users\llu22\OneDrive - Oulun yliopisto\Documents\phd\invited\Al2I6\defect\go4\vasprun.xml', index=':')
print(type(atoms))
ashow_traj(atoms)
mlab.show()