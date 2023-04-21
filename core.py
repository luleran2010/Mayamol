import numpy as np
from ase.atoms import Atoms
from ase.io import read
from ase.neighborlist import build_neighbor_list
from ase.data import covalent_radii as radii, atomic_numbers as numbers
from ase.data.colors import jmol_colors as colors
from mayavi import mlab
import argparse

from traits.api import HasTraits, Int, Range, Instance, Button, Enum, observe
from traitsui.api import View, Item, HGroup, ButtonEditor, RangeEditor
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene

# colormap for elements
ncolors = len(colors)
element_colormap = np.concatenate((colors*255, np.repeat([[255]], ncolors, axis=0)), axis=1)

def set_element_colormap(glyph):
    '''
    Set the colormap of a glyph object.
    '''
    lut = glyph.module_manager.scalar_lut_manager.lut
    lut.number_of_colors = len(colors)
    lut.table = element_colormap

class StructureVisualizer:
    # preset paths to draw the lines that represents the cell
    vertices = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1], \
                         [0,1,1], [1,0,1], [1,1,0], [1,1,1]])
    idx_lines = np.array([0, 1, 6, 2, 0, 3, 4, 2, 6, 7, 5, 3, 4, 7, 5, 1], dtype=int)

    def __init__(self, atoms: Atoms, copy: bool=False) -> None:
        '''
        atoms: the atomic structure
        copy: whether to make a copy of atoms
        '''
        self.set_structure(atoms, copy)
        
        self.figure = None
        self.mlab = None

        self.points = None
        self.bonds = None
        self.cell = None

        self.update_bonds = None

    def build_pairs(self) -> None:
        '''
        Build a list of pairs that represent the bonds between atoms.
        '''
        # compute a list of neighbors
        nl = build_neighbor_list(self.atoms, self_interaction=False, bothways=True)
        # the indices of the first and second atom of bonds
        self.pfirst = np.array(sum([[i] * len(nl.nl.neighbors[i]) for i in np.arange(len(self.atoms))], []))
        self.psecond = np.concatenate(nl.nl.neighbors, axis=0)
        # which cell the second atoms are in
        self.poffsets = np.concatenate(nl.nl.displacements, axis=0)
        # whether the second atoms are in the current cell
        self.pwithin = np.array(np.linalg.norm(self.poffsets, axis=1) > 0, dtype=int).reshape((-1,1))
        # the portions of bond length proportional to the radii of the first atoms
        self.pratios = (1 / (1+radii[self.atoms.numbers[self.psecond]]/radii[self.atoms.numbers[self.pfirst]])).reshape((-1,1))

    def set_structure(self, atoms: Atoms, copy: bool = False):
        '''
        Set the structure saved in the object.
        
        atoms: the atomic structure
        copy: whether to make a copy of atoms
        '''
        self.atoms = atoms.copy() if copy else atoms
        self.scales = radii[self.atoms.numbers]
        self.build_pairs()
        self.cell_lines = self.vertices[self.idx_lines, :] @ self.atoms.cell

    def draw(self, figure=None, mlab=mlab) -> None:
        '''
        Draw the structure using Mayavi.
        
        figure: the Mayavi figure object
        mlab: the Mayavi module
        '''
        self.mlab = mlab
        self.figure = figure

        # draw the atoms
        self.points = self.mlab.quiver3d(*self.atoms.positions.T, self.scales, self.scales, self.scales,
                                         scalars=self.atoms.numbers,
                                         mode='sphere', resolution=20, scale_factor=1,
                                         vmax=ncolors-1, vmin=0,
                                         figure=self.figure)
        self.points.glyph.color_mode = 'color_by_scalar'
        self.points.glyph.glyph_source.glyph_source.center = [0,0,0]
        set_element_colormap(self.points)

        # draw the bonds
        delta = self.get_bond_vectors()
        self.bonds = self.mlab.quiver3d(*self.atoms.positions[self.pfirst,:].T, *delta.T,
                                   mode='2ddash', line_width=10, scale_factor=1, color=(1,1,1),
                                   figure=self.figure)
        self.update_bonds = self.bonds.mlab_source.trait_set
        
        # draw the cell
        self.cell = self.mlab.plot3d(*self.cell_lines.T,
                                color=(1,1,1),
                                figure=self.figure)
        
        if self.figure is None:
            self.figure = self.mlab.gcf()

    def update_positions(self, positions:np.ndarray = None, rebuild_pairs: bool=False) -> None:
        '''
        Update the positions of atoms.
        
        positions: the positions
        rebuild_pairs: whether to rebuild the list of bond pairs
        '''
        if positions is not None:
            self.atoms.positions = positions
        if rebuild_pairs:
            self.build_pairs()
            self.update_bonds = self.bonds.mlab_source.reset

    def update_scene(self, disable_render=True) -> None:
        '''
        Update teh Mayavi scene.

        disable_render: whether to disable the render before all the object in scene are updated.
            Can speed up with large dataset.
        '''
        delta = self.get_bond_vectors()
        if disable_render and self.figure.scene is not None:
            self.figure.scene.disable_render = True
        self.points.mlab_source.trait_set(**dict(zip(['x', 'y', 'z'], self.atoms.positions[self.pfirst,:].T)))
        self.update_bonds(**dict(zip(['x', 'y', 'z'], self.atoms.positions[self.pfirst,:].T)),
                          **dict(zip(['u', 'v', 'w'], delta.T)))
        if disable_render:
            self.figure.scene.disable_render = False

    def get_bond_vectors(self):
        '''
        Generate the bond lengths and cut the ones that reach out of the cell.
        '''
        # the positions of second atoms
        delta = self.atoms.positions[self.psecond,:]+self.poffsets@self.atoms.cell
        # scale the bond length proprotional to the radii of first atoms
        delta = (1-(1-self.pratios)*self.pwithin)*delta + \
            self.pwithin*self.atoms.positions[self.pfirst,:]*(1-self.pratios)
        # substract the positions of first atoms
        delta -= self.atoms.positions[self.pfirst,:]
        return delta

def draw(atoms: Atoms, copy=False, figure=None, mlab=mlab) -> StructureVisualizer:
    '''
    Handy function to draw the atoms.

    atoms: the atomic structure
    copy: whether to make a copy of atoms
    figure: Mayavi figure to draw
    mlab: mlab module used to draw
    '''
    sv = StructureVisualizer(atoms, copy)
    sv.draw(figure=figure, mlab=mlab)
    return sv

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Animate the trajectory file')
    parser.add_argument('atoms', help='the structure file')
    args = parser.parse_args()
    atoms = read(parser.atoms)
    sv = draw(atoms)
    mlab.show()