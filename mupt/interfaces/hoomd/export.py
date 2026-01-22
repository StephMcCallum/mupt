'''Writer which converts the MuPT molecular representation to a HOOMD frame'''
'''
This code was moved over from the DPD builder. I'm working on generalizing it as an export tool to HOOMD.
TODO: remove generators, get positions from mupt representation, determine hierarchy (CG or AA, lowest level hierarchy), add angle and dihedral consideration. Add box information from mupt repr.
'''

__author__ = 'Stephanie McCallum'

import freud
import gsd, gsd.hoomd 
import time
import hoomd 
from hoomd.write import DCD
from hoomd.write import GSD
from hoomd.trigger import Periodic

from typing import (
    Generator,
    Hashable,
    Iterable,
    Iterator,
    Optional,
    Sized,
    Union,
    Sequence,
)
from numbers import Number
from collections import defaultdict
from itertools import count

import numpy as np
from scipy.spatial.transform import RigidTransform, Rotation
from networkx import all_simple_paths

from .base import PlacementGenerator
from ..mutils.iteration import flexible_iterator, sliding_window

from ..geometry.arraytypes import Shape, Dims, N
from ..geometry.measure import normalized
from ..geometry.coordinates.directions import random_unit_vector
from ..geometry.coordinates.reference import origin
from ..geometry.transforms.rigid import rigid_vector_coalignment
from ..geometry.shapes import Sphere, Ellipsoid

from ..mupr.topology import TopologicalStructure
from ..mupr.connection import Connector, TraversalDirection
from ..mupr.primitives import Primitive, PrimitiveHandle

# optional helper methods (to declutter casework from main logic)
def get_termini_handles(chain : TopologicalStructure) -> tuple[Hashable, Hashable]:
    '''
    Find the terminal node(s) of what is assumed to be a linear (path) graph
    Returns the pair of node labels of the termini (a pair of the same value twice for single-node graphs)
    '''
    termini = tuple(chain.termini)
    if len(termini) == 2:
        return termini
    elif len(termini) == 1:
        return termini[0], termini[0]
    else:
        raise ValueError('Unbranched topology must have either 1 or 2 terminal nodes')

def mupt_to_hoomd_frame(primitive : Primitive) -> frame : gsd.hoomd.frame:
    '''
    Trying to use universe of chains to set monomer positions
    primitive passed in here should be a universe primitive that has chains to loop over 
    paths are lists of handles
    If we assume chains are looped over in the same way, we can map from handles to indices
    '''
    # Initialize HOOMD Frame (initial snapshot) and periodic box
    frame = gsd.hoomd.Frame()
    
    ## Pre-allocate space for particles
    frame.particles.types = ['A'] # TODO: introduce HMT's?
    frame.particles.N = primitive.topology.number_of_nodes() # TB: would be nice to set after iterating over children, but needed to size box
    frame.particles.typeid = np.zeros(frame.particles.N)
    frame.particles.position = np.zeros((frame.particles.N, 3)) # populate with random walks

    bonds : list[tuple[int, int]] = []
    bond_types : list[str] = ['a']
    
    hoomd_chains : dict[int, tuple[int]] = dict() # for preserving chain order for orientation calc
    handle_to_particle_idx : dict[PrimitiveHandle, int] = dict()
    reference_anchor_positions : dict[int, np.ndarray[Shape[2, 3], float]] = dict()
    effective_radii : dict[int, float] = dict()
    
    particle_indexer : Iterator[int] = count(0)
    for chain_idx, chain in enumerate(primitive.topology.chains):
        head_handle, tail_handle = termini = get_termini_handles(chain)
        path : list[PrimitiveHandle] = next(all_simple_paths(chain, source=head_handle, target=tail_handle)) # raise StopIteration if no path exists

        chain_indices : list[int] = []
        for bead_handle in path:
            is_terminal : bool = ((bead_handle == head_handle) or (bead_handle == tail_handle))
            # determine unique int idx for corresponding particle in HOOMD Frame
            particle_idx = next(particle_indexer)
            chain_indices.append(particle_idx)
            handle_to_particle_idx[bead_handle] = particle_idx

            # determine reference anchor points for effective radius scaling and orientation back-calculation post-simulation
            anchor_positions = np.zeros((2, 3), dtype=float)
            bead_prim : Primitive = primitive.fetch_child(bead_handle)
            for conn_handle, conn in bead_prim.connectors.items():
                traver_dir : TraversalDirection = next(att for att in conn.anchor.attachables if isinstance(att, TraversalDirection))
                traver_dir_idx : dict[TraversalDirection, int] = {
                    TraversalDirection.ANTERO: 0,
                    TraversalDirection.RETRO: 1,
                }
                anchor_positions[traver_dir_idx[traver_dir],:] = conn.anchor.position
                if is_terminal:
                    radial_vector = conn.anchor.position - bead_prim.shape.centroid
                    diametric_anchor_pos = bead_prim.shape.centroid - radial_vector
                    anchor_positions[traver_dir_idx[TraversalDirection.complement(traver_dir)],:] = diametric_anchor_pos 
            reference_anchor_positions[particle_idx] = anchor_positions

            r_eff : float = np.linalg.norm(np.subtract(*anchor_positions)) / 2.0 
            effective_radii[particle_idx] = r_eff
        hoomd_chains[chain_idx] = tuple(chain_indices)
        
        # assign positions to LJ particle counterparts in simulation
        frame.particles.position[handle_to_particle_idx[head_handle]] = np.random.uniform( # place head randomly within box bounds
            low=(-L/2),
            high=(L/2),
            size=3,
        )
        for prim_handle_outgoing, prim_handle_incoming in sliding_window(path, 2):
            idx_outgoing, idx_incoming = idx_pair = handle_to_particle_idx[prim_handle_outgoing], handle_to_particle_idx[prim_handle_incoming]
            bonds.append(idx_pair)
            
            delta = self.bond_length * random_unit_vector()
            frame.particles.position[idx_incoming] = frame.particles.position[idx_outgoing] + delta
    
    frame.bonds.group = bonds
    frame.bonds.N = len(bonds)
    frame.bonds.types = bond_types
''' How to get these values from mupr?
    frame.angles.group = 
    frame.angles.N = 
    frame.angles.types = 

    frame.dihedrals.group = 
    frame.dihedrals.N = 
    frame.dihedrals.types = 

    frame.configuration.box = 
'''