"""
Properties of Primitives used to assess compatibility with a particular task
E.g. checking atomicity, linearity, adherence to a "standard" hierarchy, etc.
"""

__author__ = "Joseph R. Laforet Jr."
__email__ = "jola3134@colorado.edu"

from .primitives import Primitive


def is_SAAMR_compliant(prim: Primitive) -> bool:
   """
   Check whether a Primitive hierarchy is a strict SAAMR tree:
   exactly universe -> molecule -> repeat-unit -> atom (depth 3).

   SAAMR = Standard All-Atom Molecular Representation

   This checks for *strict* SAAMR compliance (all leaves are atoms at
   depth exactly 3).  Note that MDAnalysis export does **not** require
   strict SAAMR compliance — any tree with the four canonical SAAMR
   roles (UNIVERSE, SEGMENT, RESIDUE, PARTICLE) assigned can be exported,
   regardless of depth.  Use :func:`~mupt.mupr.roles.assign_SAAMR_roles`
   to tag a strict SAAMR tree, or assign roles manually for non-strict
   hierarchies.
   """

   return all(leaf.is_atom and (leaf.depth == 3) for leaf in prim.leaves)
