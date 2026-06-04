"""Strategy implementations for MuPT -> mBuild export."""

__author__ = "Stephanie McCallum"
__email__ = "stephaniemccallu@u.boisestate.edu"

from abc import ABC, abstractmethod
from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ...mupr.primitives import Primitive
from ...roles import PrimitiveRole
from .._shared.topology import (
    build_saamr_role_topology_index,
    connector_reference_sort_key,
    iter_saamr_residue_records,
    resolve_to_atom_cached,
)


@dataclass
class MBuildMolData:
    """Container for one segment's mBuild-exportable topology data."""

    segment: Primitive
    atoms: list[Primitive] = field(default_factory=list)
    atom_positions: list[np.ndarray] = field(default_factory=list)
    atom_elements: list[str] = field(default_factory=list)
    atom_resids: list[int] = field(default_factory=list)
    bonds: list[tuple[int, int]] = field(default_factory=list)
    residue_records: list = field(default_factory=list)


class MBuildExportStrategy(ABC):
    """Abstract strategy for collecting mBuild-exportable topology data."""

    @abstractmethod
    def validate(self, root: Primitive) -> None:
        """Validate role assignment and hierarchy preconditions for export."""

    @abstractmethod
    def iter_mol_data(self, root: Primitive) -> Iterator[MBuildMolData]:
        """Yield one topology dataset per SEGMENT-role node."""

    @property
    @abstractmethod
    def label(self) -> str:
        """Human-readable name for this strategy."""


class AllAtomMBuildExportStrategy(MBuildExportStrategy):
    """All-atom mBuild export strategy using PARTICLE-role nodes."""

    def __init__(self, default_atom_position: Optional[np.ndarray] = None) -> None:
        if default_atom_position is None:
            self.default_atom_position = np.array([0.0, 0.0, 0.0], dtype=float)
        else:
            default_atom_position = np.asarray(default_atom_position, dtype=float)
            if default_atom_position.shape != (3,):
                raise ValueError("default_atom_position must be a 3-dimensional vector")
            self.default_atom_position = default_atom_position

    @property
    def label(self) -> str:
        """Human-readable strategy name."""
        return "All-atom"

    def validate(self, root: Primitive) -> None:
        """Validate role assignments needed for all-atom mBuild export."""
        build_saamr_role_topology_index(root)

    def iter_mol_data(self, root: Primitive) -> Iterator[MBuildMolData]:
        """Yield one mBuild topology dataset per SEGMENT-role node."""
        index = build_saamr_role_topology_index(root)
        endpoint_cache: dict[tuple[int, object, object], Primitive] = {}
        residue_records_by_segment = {id(segment): [] for segment in index.segments}
        for residue_record in iter_saamr_residue_records(index):
            residue_records_by_segment[id(residue_record.segment)].append(residue_record)

        for segment in index.segments:
            data = MBuildMolData(segment=segment)
            atom_id_to_local: dict[int, int] = {}

            for residue_record in residue_records_by_segment[id(segment)]:
                for atom in residue_record.particles:
                    atom_id_to_local[id(atom)] = len(data.atoms)
                    data.atoms.append(atom)
                    if atom.shape is not None:
                        data.atom_positions.append(atom.shape.centroid)
                    else:
                        data.atom_positions.append(self.default_atom_position[:])
                    data.atom_elements.append(atom.element.symbol)
                    data.atom_resids.append(residue_record.residue_idx)

            data.residue_records = residue_records_by_segment[id(segment)]

            bonds_set: set[tuple[int, int]] = set()
            for node in index.bond_nodes_by_segment[id(segment)]:
                for conn_ref_pair in node.internal_connections:
                    conn_ref1, conn_ref2 = sorted(
                        conn_ref_pair,
                        key=connector_reference_sort_key,
                    )
                    atom1 = resolve_to_atom_cached(node, conn_ref1, endpoint_cache)
                    atom2 = resolve_to_atom_cached(node, conn_ref2, endpoint_cache)

                    idx1 = atom_id_to_local.get(id(atom1))
                    idx2 = atom_id_to_local.get(id(atom2))

                    if idx1 is None or idx2 is None:
                        continue

                    bond_pair = tuple(sorted((idx1, idx2)))
                    if bond_pair not in bonds_set:
                        data.bonds.append(bond_pair)
                        bonds_set.add(bond_pair)

            # Sort bonds deterministically
            if data.bonds:
                data.bonds = sorted(data.bonds)

            yield data


class CoarseGrainedMBuildExportStrategy(MBuildExportStrategy):
    """Coarse-grained mBuild export strategy using custom roles."""

    def __init__(
        self,
        cg_roles: list[str],
        default_atom_position: Optional[np.ndarray] = None,
    ) -> None:
        """Initialize coarse-grained strategy.

        Parameters
        ----------
        cg_roles : list[str]
            List of role names to export as particles (e.g., ["residue"], ["bead"]).
            When "residue" is in cg_roles, each residue becomes one CG bead (direct mapping).
            For other roles, collects only from segment's direct children.
        default_atom_position : np.ndarray, optional
            Default 3D position for nodes without explicit coordinates.
        """
        self.cg_roles = set(cg_roles)
        if default_atom_position is None:
            self.default_atom_position = np.array([0.0, 0.0, 0.0], dtype=float)
        else:
            default_atom_position = np.asarray(default_atom_position, dtype=float)
            if default_atom_position.shape != (3,):
                raise ValueError("default_atom_position must be a 3-dimensional vector")
            self.default_atom_position = default_atom_position

    @property
    def label(self) -> str:
        """Human-readable strategy name."""
        return f"Coarse-grained ({', '.join(sorted(self.cg_roles))})"

    def validate(self, root: Primitive) -> None:
        """Validate role assignments for coarse-grained export."""
        build_saamr_role_topology_index(root)

    def iter_mol_data(self, root: Primitive) -> Iterator[MBuildMolData]:
        """Yield one mBuild topology dataset per SEGMENT-role node.
        
        Special case: if "residue" in cg_roles, each residue becomes one CG bead.
        For other custom roles, collects only from segment's direct children.
        """
        index = build_saamr_role_topology_index(root)
        endpoint_cache: dict[tuple[int, object, object], Primitive] = {}

        for segment in index.segments:
            data = MBuildMolData(segment=segment)
            node_id_to_local: dict[int, int] = {}

            # Special case: "residue" in cg_roles means each residue is one CG bead
            if "residue" in self.cg_roles:
                residue_records_by_segment = {id(segment): []}
                for residue_record in iter_saamr_residue_records(index):
                    if id(residue_record.segment) == id(segment):
                        residue_records_by_segment[id(segment)].append(residue_record)
                
                data.residue_records = residue_records_by_segment[id(segment)]

                # Each residue becomes one CG particle (direct mapping)
                for residue_record in data.residue_records:
                    cg_node = residue_record.residue
                    node_id_to_local[id(cg_node)] = len(data.atoms)
                    data.atoms.append(cg_node)

                    if cg_node.shape is not None:
                        data.atom_positions.append(cg_node.shape.centroid)
                    else:
                        data.atom_positions.append(self.default_atom_position[:])

                    element_name = cg_node.label or cg_node.role.name
                    data.atom_elements.append(element_name)
                    data.atom_resids.append(residue_record.residue_idx)
            else:
                # Collect matching roles only from segment's direct children
                for child in segment.children:
                    if child.role.name in self.cg_roles:
                        node_id_to_local[id(child)] = len(data.atoms)
                        data.atoms.append(child)

                        if child.shape is not None:
                            data.atom_positions.append(child.shape.centroid)
                        else:
                            data.atom_positions.append(self.default_atom_position[:])

                        element_name = child.label or child.role.name
                        data.atom_elements.append(element_name)
                        data.atom_resids.append(0)

            # Process bonds between CG nodes
            bonds_set: set[tuple[int, int]] = set()
            for node in index.bond_nodes_by_segment[id(segment)]:
                for conn_ref_pair in node.internal_connections:
                    conn_ref1, conn_ref2 = sorted(
                        conn_ref_pair,
                        key=connector_reference_sort_key,
                    )
                    atom1 = resolve_to_atom_cached(node, conn_ref1, endpoint_cache)
                    atom2 = resolve_to_atom_cached(node, conn_ref2, endpoint_cache)

                    # Map atoms back to CG nodes
                    idx1 = self._find_parent_cg_node(atom1, node_id_to_local)
                    idx2 = self._find_parent_cg_node(atom2, node_id_to_local)

                    if idx1 is None or idx2 is None or idx1 == idx2:
                        continue

                    bond_pair = tuple(sorted((idx1, idx2)))
                    if bond_pair not in bonds_set:
                        data.bonds.append(bond_pair)
                        bonds_set.add(bond_pair)

            # Sort bonds deterministically
            if data.bonds:
                data.bonds = sorted(data.bonds)

            yield data

    def _find_parent_cg_node(
        self,
        node: Primitive,
        node_id_to_local: dict[int, int],
    ) -> Optional[int]:
        """Find the local index of the closest ancestor CG node."""
        current = node
        while current is not None:
            if id(current) in node_id_to_local:
                return node_id_to_local[id(current)]
            current = current.parent
        return None
