import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

from ase.atoms import Atoms
from ase.io import read
import h5py

from traits.api import HasTraits, Instance, Button, observe
from traitsui.api import View, Item, HGroup, ButtonEditor
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene

from core import draw

class DisplacementVisualizer(HasTraits):
    scene = Instance(MlabSceneModel, ())
    play = Button('>')
    stop = Button('||')

    def __init__(self, atoms) -> None:
        super().__init__()
        self.atoms = atoms
        self.pts = None
        self.cell = None
        self.displacements = None
        self.animator = None
    
    def draw(self, figure=None, mlab=mlab):
        self.pts, self.cell = draw(self.atoms, figure=figure, mlab=self.scene.mlab)

    def update_displacements(self, displacements):
        self.displacements = displacements
        # self.animator = None
    
    @mlab.animate(delay=30, ui=False)
    def animate(self):
        cur, n = 0, 20
        amp = np.sin(np.linspace(0, 2*np.pi, 20))
        while True:
            cur = (cur + 1) % n
            pos = self.atoms.positions+self.displacements*10*amp[cur]
            self.pts.mlab_source.trait_set(x=pos[:,0], y=pos[:,1], z=pos[:,2])
            yield

    @observe('play')
    def on_play_clicked(self, event=None):
        if self.animator is None:
            self.animator = self.animate()
        else:
            self.animator.timer.Start(30)
    
    @observe('stop')
    def on_stop_clicked(self, event=None):
        self.animator.timer.Stop()

    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene), height=250, width=300, show_label=False),
                HGroup('play', 'stop'),
                title='Animator')

class PhononVisualizer:
    def __init__(self, struct_file, band_file) -> None:
        self.atoms = read(struct_file)
        self.bands = self.read_bands(band_file)
        self.marker = None
        self.animator = DisplacementVisualizer(self.atoms)

    def read_bands(self, filename):
        # Phonon Dispersion
        with h5py.File(filename) as h5file:
            distance = h5file['distance'][:]
            frequency = h5file['frequency'][:]
            label = h5file['label'][:]
            eigenvector = h5file['eigenvector'][:]
            # nqpoint = h5file['nqpoint'][0]
            # path = h5file['path'][:]
            # segment_nqpoint = h5file['segment_nqpoint'][:]
        eigenvector = np.real(eigenvector)
        return distance, frequency, label, eigenvector

    def onpick(self, event):
        thisline = event.artist
        mousex = event.mouseevent.xdata
        mousey = event.mouseevent.ydata
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        ind = ind[np.argmin(np.sqrt((xdata[ind]-mousex)**2+(ydata[ind]-mousey)**2))]

        len_seg = self.bands[0].shape[1]+1
        len_band = len_seg*3+1
        iband = ind // len_band
        iseg = (ind - len_band*iband) // len_seg
        iqpoint = ind-len_band*iband-len_seg*iseg

        x = xdata[ind]
        y = ydata[ind]
        if self.marker is None:
            self.marker, = event.mouseevent.inaxes.plot([x], [y], 'r+', markersize=10, linewidth=5)
        else:
            self.marker.set_xdata([x])
            self.marker.set_ydata([y])
        event.mouseevent.canvas.draw()
        event.mouseevent.canvas.flush_events()

        self.update_eigenvalues(iband, iseg, iqpoint)

    def plot_bands(self, ax=None, plt=plt):
        distance, frequency, label, _ = self.bands
        nbands = frequency.shape[-1]
        dist = np.concatenate((distance, np.array([np.nan]*3).reshape((3,1))), axis=1).flatten()
        dist = dist.reshape((1,-1))
        dist = np.repeat(dist, nbands, axis=0)
        dist = np.concatenate((dist, np.array([np.nan]*nbands).reshape((nbands,1))), axis=1).flatten()

        bk = np.array([np.nan]*nbands).reshape((1,nbands))
        freq = np.concatenate((frequency[0,:,:], bk, frequency[1,:,:], bk, frequency[2,:,:], bk), axis=0).T
        freq = np.concatenate((freq, np.array([np.nan]*(nbands)).reshape((nbands,1))), axis=1).flatten()

        ticks = np.append(distance[:,0], distance[-1,-1])
        tickl = [bytes.decode(i) for i in np.append(label[:,0], label[-1,-1])]
        for i in np.arange(0, len(tickl)):
            if tickl[i] == 'Gamma':
                tickl[i] = '$\Gamma$'

        figure = None
        if ax is None:
            figure, ax = plt.subplots()
        else:
            figure = ax.gcf()

        ax.plot(dist, freq, picker=True, pickradius=5)
        ax.set_xticks(ticks, tickl)
        ax.set_ylabel('Energy/eV')
        ax.set_xlim(np.min(distance), np.max(distance))
        ax.set_ylim(np.min(frequency), np.max(frequency))
        figure.canvas.mpl_connect('pick_event', self.onpick)

    def draw(self, figure=None, mlab=mlab):
        self.animator.draw(figure, mlab)
        
    def update_eigenvalues(self, iband, iseg, iqpoint):
        egv = self.bands[3][iseg,iqpoint,:,iband]
        egv = egv.reshape((len(self.atoms), 3))
        self.animator.update_displacements(egv)
        self.animator.on_play_clicked()
    
if __name__ == '__main__':
    h5name = r'C:\Users\llu22\OneDrive - Oulun yliopisto\Documents\phd\invited\Al2I6\defect\ph4\band.hdf5'
    strname = r'C:\Users\llu22\OneDrive - Oulun yliopisto\Documents\phd\invited\Al2I6\defect\scf\CONTCAR'
    pv = PhononVisualizer(strname, h5name)
    plt.ion()
    pv.plot_bands()
    pv.draw()
    plt.show()
    pv.animator.configure_traits()