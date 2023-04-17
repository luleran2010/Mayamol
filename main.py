import numpy as np
from ase.atoms import Atoms
from ase.io import read
from ase.data import covalent_radii as radii, atomic_numbers as numbers
from ase.data.colors import jmol_colors as colors
from mayavi import mlab

from traits.api import HasTraits, Int, Range, Instance, Button, Enum, observe
from traitsui.api import View, Item, HGroup, ButtonEditor, RangeEditor
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene

ncolors = len(colors)
element_colormap = np.concatenate((colors*255, np.repeat([[255]], ncolors, axis=0)), axis=1)

def set_element_colormap(glyph):
    lut = glyph.module_manager.scalar_lut_manager.lut
    lut.number_of_colors = len(colors)
    lut.table = element_colormap

def draw(atoms: Atoms, figure=None, mlab=mlab):
    pos = atoms.positions
    s = np.array([radii[Z] for Z in atoms.numbers])
    c = np.array([colors[Z] for Z in atoms.numbers])
    points = mlab.quiver3d(pos[:,0], pos[:,1], pos[:,2], s, s, s, scalars=atoms.numbers,
                           mode='sphere', resolution=20, scale_factor=1,
                           vmax=ncolors-1, vmin=0,
                           figure=figure)
    points.glyph.color_mode = 'color_by_scalar'
    points.glyph.glyph_source.glyph_source.center = [0, 0, 0]
    set_element_colormap(points)

    # draw the outlines of the cell
    cell = []
    vertices = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1], \
                         [0,1,1], [1,0,1], [1,1,0], [1,1,1]])
    for i in np.arange(vertices.shape[0]):
        for j in np.arange(i+1, vertices.shape[0]):
            if np.linalg.norm(vertices[i,:]-vertices[j,:]) <= 1:
                line = vertices[[i,j],:] @ atoms.cell                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                l = mlab.plot3d(line[:,0], line[:,1], line[:,2], figure=figure)
                cell.append(l)
    mlab.plot3d([1, 2, np.nan, 3, 4], [1, 2, np.nan, 3, 4], [1, 2, np.nan, 3, 4], figure=figure)
    return points, cell

def draw_traj(traj: list[Atoms], figure=None, mlab=mlab):
    pts, cell = draw(traj[0], figure, mlab)
    @mlab.animate(delay=50)
    def anim():
        index = 0
        while True:
            index = (index + 1) % len(traj)
            pos = traj[index].positions
            pts.mlab_source.trait_set(x=pos[:,0], y=pos[:,1], z = pos[:,2])
            yield
    a = anim()

class TrajectoryAnimator(HasTraits):
    scene = Instance(MlabSceneModel, ())
    timestep = Range(1, 200, value = 30)
    idx_min = Int(0)
    idx_max = Int()
    index = Range(low='idx_min', high='idx_max', value=0)
    play = Button('>')
    stop = Button('||')
    goto_first = Button('<<')
    goto_last = Button('>>')

    def __init__(self, traj: list[Atoms]) -> None:
        super().__init__()
        self.traj = traj
        self.pts, self.cell = draw(traj[0], mlab=self.scene.mlab)

        self.idx_max = len(self.traj) - 1

        self.animator = None

    @mlab.animate(delay=30, ui=False)
    def animate(self):
        while True:
            self.index = (self.index + self.timestep) % len(traj)
            yield

    @observe('index')
    def on_index_changed(self, event):
        pos = self.traj[self.index].positions
        self.pts.mlab_source.trait_set(x=pos[:,0], y=pos[:,1], z = pos[:,2])

    @observe('play')
    def on_play_clicked(self, event):
        if self.animator is None:
            self.animator = self.animate()
        else:
            self.animator.timer.Start(50)
    
    @observe('stop')
    def on_stop_clicked(self, event):
        self.animator.timer.Stop()

    @observe('goto_first')
    def on_goto_first_clicked(self, event):
        self.index = 0
    
    @observe('goto_last')
    def on_goto_last_clicked(self, event):
        self.index = self.idx_max

    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene), height=250, width=300, show_label=False),
                HGroup('timestep', Item('index', editor_args={'mode': 'slider'})),
                HGroup('goto_first', 'play', 'stop', 'goto_last'),
                title='Animator'
                )

def draw_traj_adv(traj: list[Atoms]):
    animator = TrajectoryAnimator(traj)
    animator.configure_traits()
    return animator

if __name__ == '__main__':
    # mlab.figure(1, bgcolor=(0, 0, 0))
    traj = read(r'C:\Users\llu22\OneDrive - Oulun yliopisto\Documents\phd\invited\As2I6\defect\md\XDATCAR', index=':')
    # draw_traj(traj[0::30])
    # mlab.show()
    draw_traj_adv(traj)