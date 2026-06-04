"""
MuPT to mBuild Exporter

Converts MuPT Primitive hierarchies to mBuild Compound objects.

"""

__author__ = "Stephanie McCallum"
__email__ = "stephaniemccallu@u.boisestate.edu"

from collections.abc import Iterator
from typing import Optional

import numpy as np

import mbuild as mb
from ...mupr.primitives import Primitive
from .strategies import MBuildExportStrategy, AllAtomMBuildExportStrategy, MBuildMolData


def primitive_to_mbuild(
    primitive: Primitive,
    strategy: Optional[MBuildExportStrategy] = None,
) -> mb.Compound:
    """Convert Primitive hierarchy to mBuild Compound (merged multi-segment).

    Returns all segments merged into one hierarchy:
    root → segment_1 → residue_1 → particles
         → segment_2 → residue_2 → particles
         ...

    Parameters
    ----------
    primitive : Primitive
        Root Primitive of the molecular system to export.
        Must have SAAMR roles assigned (UNIVERSE, SEGMENT, RESIDUE, PARTICLE).
    strategy : MBuildExportStrategy, optional
        Export strategy (all-atom, coarse-grained, etc.).
        Defaults to AllAtomMBuildExportStrategy.

    Returns
    -------
    mb.Compound
        mBuild Compound hierarchy with particles and bonds.
    """
    if strategy is None:
        strategy = AllAtomMBuildExportStrategy()

    strategy.validate(primitive)
    root = mb.Compound(name=primitive.label or "system")

    for segment_idx, compound in enumerate(
        primitive_to_mbuild_compounds(primitive, strategy)
    ):
        root.add(compound)

    return root


def primitive_to_mbuild_compounds(
    primitive: Primitive,
    strategy: Optional[MBuildExportStrategy] = None,
) -> Iterator[mb.Compound]:
    """Yield one mBuild Compound per SEGMENT-role node.

    Matches RDKit pattern for multi-chain systems.

    Parameters
    ----------
    primitive : Primitive
        Root Primitive of the molecular system to export.
    strategy : MBuildExportStrategy, optional
        Export strategy (all-atom, coarse-grained, etc.).
        Defaults to AllAtomMBuildExportStrategy.

    Yields
    ------
    mb.Compound
        One Compound per segment in the hierarchy.
    """
    if strategy is None:
        strategy = AllAtomMBuildExportStrategy()

    strategy.validate(primitive)

    for segment_idx, mol_data in enumerate(strategy.iter_mol_data(primitive)):
        compound = _build_mbuild_compound(mol_data, segment_idx)
        yield compound


def _build_mbuild_compound(
    data: MBuildMolData,
    segment_idx: int,
) -> mb.Compound:
    """Construct mBuild Compound from collected topology data.

    Parameters
    ----------
    data : MBuildMolData
        Collected topology data for one segment.
    segment_idx : int
        Index of this segment.

    Returns
    -------
    mb.Compound
        mBuild Compound with hierarchical structure preserved.
    """
    root = mb.Compound(name=data.segment.label or f"segment_{segment_idx}")

    # Add hierarchical structure: segment → residues → particles
    atom_id_to_particle: dict[int, mb.Particle] = {}

    for residue_record in data.residue_records:
        res_compound = mb.Compound(name=residue_record.residue.label)
        root.add(res_compound)

        # Add particles for atoms in this residue
        for atom in residue_record.particles:
            try:
                global_idx = data.atoms.index(atom)
            except ValueError:
                continue

            if global_idx < len(data.atom_elements):
                particle = mb.Particle(
                    name=data.atom_elements[global_idx],
                    pos=data.atom_positions[global_idx],
                    element=data.atom_elements[global_idx],
                )
                res_compound.add(particle)
                atom_id_to_particle[id(atom)] = particle

    # Add bonds
    for idx1, idx2 in data.bonds:
        if idx1 < len(data.atoms) and idx2 < len(data.atoms):
            p1 = atom_id_to_particle.get(id(data.atoms[idx1]))
            p2 = atom_id_to_particle.get(id(data.atoms[idx2]))
            if p1 is not None and p2 is not None:
                root.add_bond([p1, p2])


    return root
