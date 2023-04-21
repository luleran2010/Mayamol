import numpy as np
from ase.atoms import Atoms
from ase.io import read
from ase.data import covalent_radii as radii, atomic_numbers as numbers
from ase.data.colors import jmol_colors as colors
from mayavi import mlab
import argparse

from traits.api import HasTraits, Int, Range, Instance, Button, observe
from traitsui.api import View, Item, HGroup
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.tools.animator import Animator
from mayavi.core.ui.mayavi_scene import MayaviScene

from Mayamol.core import StructureVisualizer

def draw_traj(traj: list[Atoms], figure=None, mlab=mlab) -> Animator:
    '''
    Handy function to draw a trajectory in Mayavi.

    traj: the trajectory
    figure: mayavi figure to draw
    mlab: mlab module to draw
    '''
    # draw the first image of the trajectory
    sv = StructureVisualizer(traj[0], copy=True)
    sv.draw(figure, mlab)
    # create the animation function
    @mlab.animate(delay=30)
    def anim():
        index = 0
        while True:
            index = (index + 1) % len(traj)
            sv.update_positions(traj[index].positions, rebuild_pairs=True)
            sv.update_scene()
            yield
    a = anim() # animate
    return a

class TrajectoryAnimator(HasTraits):
    '''
    The Traits class that provides more advanced control over the animation of trajectory.
    '''
    scene = Instance(MlabSceneModel, ())                  # the Mayavi scene model that offers 3D visualization
    timestep = Range(1, 200, value = 30)                  # the number of image to skip between each timeframe
    idx_min = Int(0)                                      # the minimum image index
    idx_max = Int()                                       # the maximum image index
    index = Range(low='idx_min', high='idx_max', value=0) # image index
    play = Button('>')
    stop = Button('||')
    goto_first = Button('<<')
    goto_last = Button('>>')

    def __init__(self, traj: list[Atoms]) -> None:
        '''
        traj: the trajectory
        '''
        super().__init__()
        self.traj = traj
        self.idx_max = len(self.traj) - 1
        self.sv = StructureVisualizer(traj[0], copy=True)

        self.mlab = self.scene.mlab
        self.figure = self.scene.mlab.gcf()
        self.animator = None

    def draw(self, figure=None, mlab=None) -> None:
        '''
        Draw the first image of trajectory in Mayavi,

        figure: mayavi figure to draw
        mlab: mlab module to draw
        '''
        if figure is None:
            figure = self.figure
        if mlab is None:
            mlab = self.mlab
        self.sv.draw(mlab=self.mlab)

    @mlab.animate(delay=30, ui=False)
    def animate(self):
        '''
        The animation function,
        '''
        while True:
            # update the image index
            self.index = (self.index + self.timestep) % len(self.traj)
            yield

    @observe('index')
    def on_index_changed(self, event):
        if self.figure.scene is None:
            return
        self.sv.update_positions(self.traj[self.index].positions, rebuild_pairs=True)
        self.sv.update_scene()

    @observe('play')
    def on_play_clicked(self, event):
        if self.animator is None:
            self.animator = self.animate()
        else:
            self.animator.timer.Start(30)
    
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

def draw_traj_adv(traj: list[Atoms]) -> TrajectoryAnimator:
    '''
    Handy function to draw the trajectory with advanced controls
    
    traj: the trajectory
    '''
    animator = TrajectoryAnimator(traj)
    animator.draw()
    animator.configure_traits()
    return animator

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Animate the trajectory file')
    parser.add_argument('traj', help='the trajector file')
    args = parser.parse_args()
    traj = read(args.traj, index=':')
    a = draw_traj_adv(traj)