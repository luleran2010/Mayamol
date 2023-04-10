import numpy as np
from ase.atoms import Atoms
from ase.io import read
from ase.data import covalent_radii as radii, atomic_numbers as numbers
from ase.data.colors import jmol_colors as colors
from mayavi import mlab

def ashow(atoms: Atoms, figure=None):
    for symbol in set(atoms.symbols):
        number = numbers[symbol]
        mask = atoms.symbols == symbol
        mlab.points3d(atoms.positions[mask,0], atoms.positions[mask,1], atoms.positions[mask,2],
                      scale_factor=radii[number]/radii[1]/2,
                      color=tuple(colors[number]),
                      resolution=50,
                      figure=figure)
    vertices = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1], \
                         [0,1,1], [1,0,1], [1,1,0], [1,1,1]])
    for i in np.arange(vertices.shape[0]):
        for j in np.arange(i+1, vertices.shape[0]):
            if np.linalg.norm(vertices[i,:]-vertices[j,:]) <= 1:
                line = vertices[[i,j],:] @ atoms.cell
                mlab.plot3d(line[:,0], line[:,1], line[:,2])

    mlab.show()

atoms = read('test/c2db-6222.xyz')
ashow(atoms)