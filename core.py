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