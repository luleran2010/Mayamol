import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import argparse

from ase.atoms import Atoms
from ase.io import read
import h5py

from traits.api import HasTraits, Float, Array, Instance, Button, observe
from traitsui.api import View, Item, HGroup, VGroup
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene

from matplotlib.backend_bases import PickEvent
from matplotlib.axes import Axes
import addcopyfighandler

from Mayamol.core import StructureVisualizer

class DisplacementVisualizer(HasTraits):
    '''
    Class to visualize phonon displacement in an oscillating movement.
    '''
    scene = Instance(MlabSceneModel, ())
    amplitude = Float(1)                             # the amplitude of the oscilation
    length = Float(1)                                # the scale of lengths of the arrows
    rep = Array(np.int32, (3,), np.array([1, 1, 1])) # repeatition of the cell
    play = Button('>')
    stop = Button('||')

    def __init__(self, atoms: Atoms, displacements: np.ndarray = None) -> None:
        '''
        atoms: the atomic structure
        displacements: the initial displacements, default to a zero array
        '''
        super().__init__()
        self.atoms = atoms
        self.atoms_rep = self.atoms.repeat(self.rep)
        self.disps = np.zeros((len(atoms), 3)) if displacements is None else displacements
        self.disps_rep = np.tile(self.disps, (np.prod(self.rep),1))

        self.osci_x = np.linspace(0, 2*np.pi, 20)

        self.sv = StructureVisualizer(atoms, copy=True)
        self.mlab = self.scene.mlab
        self.figure = self.mlab.gcf()
        self.animator = None

    @observe('rep')
    def on_rep_change(self, event=None):
        self.atoms_rep = self.atoms.repeat(self.rep)
        self.disps_rep = np.tile(self.disps, (np.prod(self.rep),1))
        self.mlab.clf(self.sv.figure)
        self.sv.set_structure(self.atoms_rep, copy=True)
        self.draw()

    def draw(self, figure=None, mlab=None) -> None:
        '''
        Draw the structure and displacements in the scene.
        '''
        if figure is None:
            figure = self.figure
        if mlab is None:
            mlab = self.mlab
        self.sv.draw(figure=figure, mlab=mlab)
        self.arrows = mlab.quiver3d(*self.atoms_rep.positions.T, *self.disps_rep.T,
                                    mode='2darrow', line_width=2, scale_factor=1,
                                    figure=self.figure)

    def update_displacements(self, disp: np.ndarray):
        '''
        Update the displacements and also its repeatitive counterpart.
        '''
        self.disps = disp
        self.disps_rep = np.tile(self.disps, (np.prod(self.rep),1))

    @mlab.animate(delay=30, ui=False)
    def animate(self):
        cur, n = 0, 20
        # compute the amplitude of the oscilation in each timeframe and also its derivative
        osci = np.sin(self.osci_x)
        deriv = np.cos(self.osci_x)
        while True:
            if self.figure.scene is None:
                return
            # disable the render to speed up
            self.figure.scene.disable_render = True
            cur = (cur + 1) % n
            positions = self.atoms_rep.positions+self.disps_rep*self.amplitude*osci[cur]
            vectors = self.disps_rep*self.amplitude*self.length*deriv[cur]
            # update the positions of atoms
            self.sv.update_positions(positions)
            self.sv.update_scene(disable_render=False)
            # update the arrows
            self.arrows.mlab_source.trait_set(**dict(zip(['x', 'y', 'z'], positions.T)),
                                              **dict(zip(['u', 'v', 'w'], vectors.T)))
            self.figure.scene.disable_render = False
            yield

    @observe('play')
    def on_play_clicked(self, event=None):
        if np.all(self.disps == 0):
            return
        if self.animator is None:
            self.animator = self.animate()
        else:
            self.animator.timer.Start(30)
    
    @observe('stop')
    def on_stop_clicked(self, event=None):
        if self.animator is not None:
            self.animator.timer.Stop()

    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene), height=250, width=300, show_label=False),
                HGroup('rep', VGroup('amplitude',
                                     'length',
                                     HGroup('play', 'stop'))
                       ),
                title='Animator')

def read_bands(filename: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Phonon Dispersion
    with h5py.File(filename) as h5file:
        distance = h5file['distance'][:]   # the K-path
        frequency = h5file['frequency'][:] # the frequencies
        label = h5file['label'][:]         # labels of special points
        eigenvector = h5file['eigenvector'][:]
        # nqpoint = h5file['nqpoint'][0]
        # path = h5file['path'][:]
        # segment_nqpoint = h5file['segment_nqpoint'][:]
    eigenvector = np.real(eigenvector)     # displacements are the real parts of eigen vectors
    return distance, frequency, label, eigenvector

def plot_bands(bands, ax: Axes=None) -> None:
    distance, frequency, label, _ = bands
    nbands = frequency.shape[-1]
    nseg = distance.shape[0]
    # add np.nan to break the lines
    dist = np.concatenate((distance, np.array([np.nan]*nseg).reshape((nseg,1))), axis=1).flatten()
    dist = dist.reshape((1,-1))
    dist = np.repeat(dist, nbands, axis=0)
    dist = np.concatenate((dist, np.array([np.nan]*nbands).reshape((nbands,1))), axis=1).flatten()

    # add np.nan to break the lines
    bk = np.array([np.nan]*nbands).reshape((1,nbands))
    freq = np.concatenate(sum([[frequency[i,:,:], bk] for i in np.arange(frequency.shape[0])], []), axis=0).T
    freq = np.concatenate((freq, np.array([np.nan]*(nbands)).reshape((nbands,1))), axis=1).flatten()

    # retrieve the positions and labels of special points
    ticks = np.append(distance[:,0], distance[-1,-1])
    tickl = [bytes.decode(i) for i in np.append(label[:,0], label[-1,-1])]
    for i in np.arange(0, len(tickl)):
        if tickl[i] == 'Gamma':
            tickl[i] = '$\Gamma$'

    if ax is None:
        figure, ax = plt.subplots()

    ax.plot(dist, freq, picker=True, pickradius=5) # the dispersion
    ax.set_xticks(ticks, tickl) # ticks at special points
    ax.set_ylabel('Frequency/THz')
    ax.set_xlim(np.min(distance), np.max(distance))
    min_freq, max_freq = np.min(frequency), np.max(frequency)
    ax.set_ylim(min_freq, max_freq)
    # vertical lines at special points
    ax.vlines(ticks[1:-1], min_freq, max_freq, colors='gray', linestyles='dashed')
    return ax

class PhononVisualizer:
    '''
    Class to show the phonon dispersion and the displacements.
    '''
    def __init__(self, atoms, bands) -> None:
        self.atoms = atoms
        self.bands = bands
        self.marker = None
        self.animator = DisplacementVisualizer(self.atoms)

    def onpick(self, event: PickEvent) -> None:
        '''
        Callback when picking a data point in the phonon dispersion.
        '''
        thisline = event.artist
        mousex = event.mouseevent.xdata
        mousey = event.mouseevent.ydata
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        ind = ind[np.argmin(np.sqrt((xdata[ind]-mousex)**2+(ydata[ind]-mousey)**2))]

        # compute the index of the picked frequency
        len_seg = self.bands[0].shape[1]+1
        len_band = len_seg*3+1
        iband = ind // len_band
        iseg = (ind - len_band*iband) // len_seg
        iqpoint = ind-len_band*iband-len_seg*iseg

        x = xdata[ind]
        y = ydata[ind]
        # mark the position in the phonon dispersion
        if self.marker is None:
            self.marker, = event.mouseevent.inaxes.plot([x], [y], 'r+', markersize=10, linewidth=5)
        else:
            self.marker.set_xdata([x])
            self.marker.set_ydata([y])
        event.mouseevent.canvas.draw()
        event.mouseevent.canvas.flush_events()

        # show it in the displacement visualizer
        self.update_eigenvalues(iband, iseg, iqpoint)

    def plot_bands(self, ax: Axes=None) -> None:
        ax = plot_bands(self.bands, ax)
        figure = ax.figure
        # hook the callback to pick event
        figure.canvas.mpl_connect('pick_event', self.onpick)

    def draw(self, figure=None, mlab=None) -> None:
        self.animator.draw(figure, mlab)
        
    def update_eigenvalues(self, iband: int, iseg: int, iqpoint: int) -> None:
        '''
        Update the displacement in the displacement visualizer.

        iband, iseg, iqpoint: index of the frequency in phonon disperison
        '''
        egv = self.bands[3][iseg,iqpoint,:,iband]
        egv = egv.reshape((len(self.atoms), 3))
        self.animator.update_displacements(egv)
        self.animator.on_play_clicked()
    
def draw_phonons(atoms, bands, mlab_figure=None, mlab=None, plt_ax=None):
    '''
    Handy function to show phonon dispersion and displacements
    '''
    pv = PhononVisualizer(atoms, bands)
    plt.ion()
    pv.plot_bands(ax=plt_ax)
    pv.draw(figure=mlab_figure, mlab=mlab)
    plt.show()
    pv.animator.configure_traits()
    return pv

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Phonon dispersion & displacement visualizer')
    parser.add_argument('structure', help='the structure file')
    parser.add_argument('bands', help='the Phonopy-generated bands.h5 file')
    args = parser.parse_args()
    atoms = read(args.structure)
    bands = read_bands(args.bands)
    draw_phonons(atoms, bands)