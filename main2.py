import numpy as np
from mayavi import mlab
from traits.api import HasTraits, Range, Instance, Button, observe
from traitsui.api import View, Item, HGroup, RangeEditor
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from pyface.timer.timer import Timer
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
    pos = atoms.positions
    s = np.array([radii[Z] for Z in atoms.numbers])
    points = mlab.points3d(pos[:,0], pos[:,1], pos[:,2], s, scale_mode='vector')
    c = 
    return points

class Mayamol(HasTraits):
    step_low = 0
    step_high = 1
    step = Range('step_low', 'step_high', value=0)
    scene = Instance(MlabSceneModel, ())
    play = Button('play')

    def __init__(self, traj: list[Atoms] | Atoms) -> None:
        super().__init__()
        self.traj = traj
        self.step_high = len(self.traj) - 1
        self.mlab = self.scene.mlab
        self.cell = draw_cell(traj[0].cell, mlab=self.mlab)
        self.points = draw_atoms(traj[0], mlab=self.mlab)

    @observe('step')
    def update_plot(self, event):
        pos = self.traj[self.step].positions
        self.points.mlab_source.trait_set(x=pos[:,0], y=pos[:,1], z=pos[:,2])
            
    # @observe('play')
    # def update_plot(self, event):
    #     def anim():
    #         for i in range(len(self.traj)):
    #             self.trait_set(step=i)
    #             yield
    #     a = anim()
    #     t = Timer(50, next(a))


traj = read(r'C:\Users\llu22\OneDrive - Oulun yliopisto\Documents\phd\invited\Al2I6\defect\go4\vasprun.xml', index=':')
mayamol = Mayamol(traj)

view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene), height=250, width=300, show_label=False),
            HGroup('_',
                    'step',
                    'play'
                )
            )

mayamol.configure_traits(view=view)