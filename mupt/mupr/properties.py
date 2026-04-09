"""
Properties of Primitives used to assess compatibility with a particular task
E.g. checking atomicity, linearity, adherence to a "standard" hierarchy, etc.
"""

__author__ = "Joseph R. Laforet Jr."
__email__ = "jola3134@colorado.edu"

from .primitives import Primitive


def is_SAAMR_compliant(prim: Primitive) -> bool:
   """
   Check whether a Primitive hierarchy is organized
   as universe -> molecule -> repeat-unit -> atom
   
   SAAMR = Standard All-Atom Molecular Representation

   A root-only or atom-less tree returns False because a valid SAAMR
   hierarchy requires all leaves to be atoms at depth 3.  Note that
   ``prim.leaves`` always includes at least ``prim`` itself (anytree
   treats childless nodes as their own leaf), so the ``all()`` predicate
   is never evaluated over an empty iterable.
   """

   return all(leaf.is_atom and (leaf.depth == 3) for leaf in prim.leaves)
