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

from core import draw

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
            self.index = (self.index + self.timestep) % len(self.traj)
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