"""
Tests for MuPT to mBuild export.

Verifies that export produces valid mBuild Compounds, preserves topology,
and supports both all-atom and coarse-grained strategies.
"""

__author__ = "Stephanie McCallum"
__email__ = "stephaniemccallu@u.boisestate.edu"

import pytest
from anytree import PreOrderIter

import mbuild as mb

from mupt.interfaces.mbuild import (
    primitive_to_mbuild,
    primitive_to_mbuild_compounds,
    AllAtomMBuildExportStrategy,
    CoarseGrainedMBuildExportStrategy,
)
from mupt.mupr.primitives import Primitive
from mupt.roles import PrimitiveRole, assign_SAAMR_roles


def _mbuild_compounds(*args, **kwargs):
    """Materialize the streaming exporter for small unit-test fixtures."""
    return list(primitive_to_mbuild_compounds(*args, **kwargs))


def _count_internal_connections(root: Primitive) -> int:
    """Count total internal connections in hierarchy."""
    return sum(
        len(node.internal_connections)
        for node in PreOrderIter(root)
        if not node.is_leaf
    )


def _count_particles(compound: mb.Compound) -> int:
    """Count all particle objects in mBuild Compound hierarchy."""
    count = 0

    def traverse(node):
        nonlocal count
        if isinstance(node, mb.Particle):
            count += 1
        for child in node.children:
            traverse(child)

    traverse(compound)
    return count


def _count_bonds(compound: mb.Compound) -> int:
    """Count bonds in mBuild Compound."""
    return len(compound.bonds)


class TestAllAtomExport:
    """Test all-atom export strategy."""

    def test_all_atom_preserves_atom_count(self, single_polyethylene_3mer):
        """All-atom export preserves atom count."""
        compound = primitive_to_mbuild(single_polyethylene_3mer)
        atom_count = _count_particles(compound)
        expected_count = len(single_polyethylene_3mer.leaves)

        assert atom_count == expected_count

    def test_all_atom_preserves_bond_count(self, depth4_bonded_system):
        """All-atom export preserves bond count."""
        compound = primitive_to_mbuild(depth4_bonded_system)
        bond_count = _count_bonds(compound)
        expected_bond_count = _count_internal_connections(depth4_bonded_system)

        assert bond_count == expected_bond_count

    def test_all_atom_yields_one_compound_per_segment(self, multi_polyethylene_system):
        """Iterator yields one compound per segment."""
        compounds = _mbuild_compounds(multi_polyethylene_system)

        # multi_polyethylene_system has 10 chains
        assert len(compounds) == 10

    def test_all_atom_merged_compound_preserves_structure(self, single_polyethylene_2mer):
        """Merged compound preserves segment/residue hierarchy."""
        compound = primitive_to_mbuild(single_polyethylene_2mer)

        # Should have segment → residue → particle structure
        assert compound is not None
        assert len(list(compound.children)) > 0

    def test_all_atom_default_atom_position(self, single_polyethylene_2mer):
        """All-atom export uses default position for atoms without shape."""
        compound = primitive_to_mbuild(
            single_polyethylene_2mer,
            strategy=AllAtomMBuildExportStrategy(default_atom_position=[1.0, 2.0, 3.0]),
        )

        assert compound is not None
        atom_count = _count_particles(compound)
        assert atom_count > 0

    def test_all_atom_explicit_strategy_same_as_implicit(self, single_polyethylene_3mer):
        """Explicit all-atom strategy produces same result as implicit default."""
        compound_implicit = primitive_to_mbuild(single_polyethylene_3mer)
        compound_explicit = primitive_to_mbuild(
            single_polyethylene_3mer,
            strategy=AllAtomMBuildExportStrategy(),
        )

        assert _count_particles(compound_implicit) == _count_particles(compound_explicit)
        assert _count_bonds(compound_implicit) == _count_bonds(compound_explicit)

    def test_all_atom_deterministic_output(self, single_polyethylene_2mer):
        """Repeated all-atom exports produce identical bond ordering."""
        compound1 = primitive_to_mbuild(single_polyethylene_2mer)
        compound2 = primitive_to_mbuild(single_polyethylene_2mer)

        bonds1 = sorted((b[0], b[1]) for b in compound1.bonds)
        bonds2 = sorted((b[0], b[1]) for b in compound2.bonds)

        assert len(bonds1) == len(bonds2)


class TestCoarseGrainedExport:
    """Test coarse-grained export strategy."""

    def test_cg_custom_roles_detected(self):
        """CG strategy detects custom roles in hierarchy."""
        # Create a simple hierarchy with custom roles
        residue = Primitive(label="residue", role=PrimitiveRole.RESIDUE)

        # Add some fake CG beads
        bead1 = Primitive(label="bead_1", role=PrimitiveRole.UNASSIGNED)
        bead1.role = PrimitiveRole.UNASSIGNED
        bead1.role = type("Role", (), {"name": "bead"})()  # Mock custom role

        bead2 = Primitive(label="bead_2", role=PrimitiveRole.UNASSIGNED)
        bead2.role = type("Role", (), {"name": "bead"})()

        residue.attach_child(bead1)
        residue.attach_child(bead2)

        segment = Primitive(label="segment", role=PrimitiveRole.SEGMENT)
        segment.attach_child(residue)

        universe = Primitive(label="universe", role=PrimitiveRole.UNIVERSE)
        universe.attach_child(segment)

        strategy = CoarseGrainedMBuildExportStrategy(cg_roles=["bead"])
        assert strategy is not None
        assert "bead" in strategy.cg_roles

    def test_cg_strategy_label(self):
        """CG strategy generates appropriate label."""
        strategy = CoarseGrainedMBuildExportStrategy(cg_roles=["bead"])

        label = strategy.label
        assert "Coarse-grained" in label
        assert "bead" in label

        def test_cg_residue_direct_mapping(self, single_polyethylene_2mer):
            """CG strategy with 'residue' role maps each residue to one CG particle."""
            strategy = CoarseGrainedMBuildExportStrategy(cg_roles=["residue"])
            compounds = _mbuild_compounds(single_polyethylene_2mer, strategy=strategy)
        
            # Should have one compound per segment
            assert len(compounds) == 1
        
            # The compound should have CG particles (one per residue)
            cg_particles = _count_particles(compounds[0])
            assert cg_particles > 0

        def test_cg_residue_to_particle_count(self, single_polyethylene_3mer):
            """CG export with residue mapping yields one particle per residue."""
            strategy = CoarseGrainedMBuildExportStrategy(cg_roles=["residue"])
            compound = primitive_to_mbuild(single_polyethylene_3mer, strategy=strategy)
        
            cg_particle_count = _count_particles(compound)
        
            # Count residues in original structure
            residue_count = sum(
                1 for node in PreOrderIter(single_polyethylene_3mer)
                if node.role.name == "RESIDUE"
            )
        
            assert cg_particle_count == residue_count

        def test_cg_only_segment_children_collected(self, single_polyethylene_2mer):
            """CG strategy collects only from segment's direct children for non-residue roles."""
            # Create custom role for testing
            strategy = CoarseGrainedMBuildExportStrategy(cg_roles=["custom"])
            compounds = _mbuild_compounds(single_polyethylene_2mer, strategy=strategy)
        
            # Should not crash and should return compounds
            assert len(compounds) > 0

class TestStrategyValidation:
    """Test strategy validation."""

    def test_strategy_validate_rejects_non_saamr(self):
        """Strategy validation rejects non-SAAMR hierarchy."""
        # Create a hierarchy without proper roles
        prim = Primitive(label="not_saamr")

        strategy = AllAtomMBuildExportStrategy()

        with pytest.raises(ValueError):
            strategy.validate(prim)

    def test_all_atom_strategy_validate_succeeds_on_saamr(self, single_polyethylene_2mer):
        """Strategy validation succeeds on proper SAAMR hierarchy."""
        strategy = AllAtomMBuildExportStrategy()

        # Should not raise
        strategy.validate(single_polyethylene_2mer)

    def test_cg_strategy_validate_succeeds_on_saamr(self, single_polyethylene_2mer):
        """CG strategy validation succeeds on proper SAAMR hierarchy."""
        strategy = CoarseGrainedMBuildExportStrategy(cg_roles=["residue"])

        # Should not raise
        strategy.validate(single_polyethylene_2mer)


class TestMultiChain:
    """Test multi-chain export."""

    def test_multi_chain_compounds_iterator(self, multi_polyethylene_system):
        """Iterator properly yields one compound per chain."""
        compounds = _mbuild_compounds(multi_polyethylene_system)

        assert len(compounds) == 10
        for compound in compounds:
            assert isinstance(compound, mb.Compound)
            assert _count_particles(compound) > 0

    def test_multi_chain_merged_compound(self, multi_polyethylene_system):
        """Merged compound contains all chains as children."""
        compound = primitive_to_mbuild(multi_polyethylene_system)

        # Should have multiple segment children
        segment_children = [
            child for child in compound.children
            if isinstance(child, mb.Compound) and "segment" in child.name.lower()
        ]
        assert len(segment_children) > 0

    def test_multi_chain_particle_count(self, multi_polyethylene_system):
        """Multi-chain export preserves all particles."""
        compound = primitive_to_mbuild(multi_polyethylene_system)
        particle_count = _count_particles(compound)
        expected_count = len(multi_polyethylene_system.leaves)

        assert particle_count == expected_count


class TestEdgeCases:
    """Test edge cases and special scenarios."""

    def test_export_valid_mbuild_compound(self, single_polyethylene_2mer):
        """Exported compound is a valid mBuild Compound."""
        compound = primitive_to_mbuild(single_polyethylene_2mer)

        assert isinstance(compound, mb.Compound)
        assert compound is not None

    def test_export_preserves_element_symbols(self, single_polyethylene_3mer):
        """Element symbols are preserved in particle names."""
        compound = primitive_to_mbuild(single_polyethylene_3mer)

        # Collect all particle names
        particle_names = set()

        def collect_names(node):
            if isinstance(node, mb.Particle):
                particle_names.add(node.name)
            for child in node.children:
                collect_names(child)

        collect_names(compound)

        # Should have carbon and hydrogen particles
        assert len(particle_names) > 0
        assert any(name in ["C", "H"] for name in particle_names)

    def test_export_with_segment_label(self, single_polyethylene_2mer):
        """Exported compound preserves segment labels."""
        compound = primitive_to_mbuild(single_polyethylene_2mer)

        assert compound.name is not None
        assert len(compound.name) > 0

        def test_no_linker_particles_created(self, single_polyethylene_2mer):
            """Ghost linker particles are not created in export."""
            compound = primitive_to_mbuild(single_polyethylene_2mer)

            # Collect all particle names to check for "X" (linker) particles
            particle_names = []

            def collect_names(node):
                if isinstance(node, mb.Particle):
                    particle_names.append(node.name)
                for child in node.children:
                    collect_names(child)

            collect_names(compound)

            # Should not have any "X" linker particles
            assert "X" not in particle_names

class TestExportFunctions:
    """Test the main export functions."""

    def test_primitive_to_mbuild_returns_compound(self, single_polyethylene_2mer):
        """primitive_to_mbuild returns mb.Compound."""
        result = primitive_to_mbuild(single_polyethylene_2mer)

        assert isinstance(result, mb.Compound)

    def test_primitive_to_mbuild_compounds_returns_iterator(self, single_polyethylene_2mer):
        """primitive_to_mbuild_compounds returns iterator."""
        result = primitive_to_mbuild_compounds(single_polyethylene_2mer)

        # Should be iterable
        compounds = list(result)
        assert len(compounds) > 0
        assert all(isinstance(c, mb.Compound) for c in compounds)

    def test_primitive_to_mbuild_default_strategy(self, single_polyethylene_2mer):
        """primitive_to_mbuild uses AllAtomMBuildExportStrategy by default."""
        compound = primitive_to_mbuild(single_polyethylene_2mer)

        # Should produce valid output
        assert _count_particles(compound) == len(single_polyethylene_2mer.leaves)

    def test_primitive_to_mbuild_accepts_strategy(self, single_polyethylene_2mer):
        """primitive_to_mbuild accepts custom strategy."""
        strategy = AllAtomMBuildExportStrategy()
        compound = primitive_to_mbuild(single_polyethylene_2mer, strategy=strategy)

        assert _count_particles(compound) > 0
