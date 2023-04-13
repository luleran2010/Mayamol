import numpy as np

from mayavi import mlab
from traits.api import HasTraits, Range, Instance, on_trait_change
from traitsui.api import View, Item, HGroup, Action, UIInfo, Handler
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene

from ase.atoms import Cell, Atoms
from ase.io import read
from ase.data import covalent_radii as radii, atomic_numbers as numbers
from ase.data.colors import jmol_colors as colors

def draw_cell(cell: np.ndarray | Cell, figure=None, mlab=mlab):
    lines = []
    vertices = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1], \
                        [0,1,1], [1,0,1], [1,1,0], [1,1,1]])
    for i in np.arange(vertices.shape[0]):
        for j in np.arange(i+1, vertices.shape[0]):
            if np.linalg.norm(vertices[i,:]-vertices[j,:]) <= 1:
                line = vertices[[i,j],:] @ cell
                l = mlab.plot3d(line[:,0], line[:,1], line[:,2], figure=figure)
                lines.append(l)
    return lines

def draw_atoms(atoms: Atoms, figure=None, mlab=mlab) -> dict:
    points = {}
    for symbol in set(atoms.symbols):
        number = numbers[symbol]
        mask = atoms.symbols == symbol
        p = mlab.points3d(atoms.positions[mask,0], atoms.positions[mask,1], atoms.positions[mask,2],
                      scale_factor=radii[number]/radii[1]/2,
                      color=tuple(colors[number]),
                      resolution=50,
                      figure=figure)
        points[symbol] = p
    return points

class MayamolHandler(Handler):
    def play(info: UIInfo):
        self.mod

class Mayamol(HasTraits):
    # image = Range(0, 10, 6)
    scene = Instance(MlabSceneModel, ())

    def __init__(self, traj: list[Atoms] | Atoms) -> None:
        super().__init__()
        image = Range(0, len(traj)-1, value=0)
        self.add_trait("image", image)
        self.traj = traj
        self.mlab = self.scene.mlab
        self.cell = draw_cell(traj[0].cell, mlab=self.mlab)
        self.points = draw_atoms(traj[0], mlab=self.mlab)

    @on_trait_change('image')
    def update_plot(self):
        image = self.traj[self.image]
        for symbol in self.points.keys():
            mask = image.symbols == symbol
            self.points[symbol].mlab_source.trait_set(x=image.positions[mask,0],
                                                      y=image.positions[mask,1],
                                                      z=image.positions[mask,2])

    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene), height=250, width=300, show_label=False),
                HGroup('_',
                       'image'
                    ),
                handler=MayamolHandler(),
                buttons=[Action(name='Play', action='play')]
                )

traj = read(r'C:\Users\llu22\OneDrive - Oulun yliopisto\Documents\phd\invited\Al2I6\defect\go4\vasprun.xml', index=':')
mayamol = Mayamol(traj)
mayamol.configure_traits()