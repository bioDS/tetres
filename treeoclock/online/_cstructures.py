from ctypes import Structure, POINTER, c_long
from treeoclock.trees._ctrees import TREE


class PRECOMP(Structure):
    _fields_ = [('tree_changes', POINTER(POINTER(TREE))), ('rank_changes', POINTER(POINTER(c_long))),
                ('rank_identities', POINTER(POINTER(c_long))), ('dist', c_long)]

    def __init_(self):
        self.tree_changes = None
        self.rank_changes = None
        self.rank_identities = None
        self.dist = 0


class UPDATE(Structure):
    _fields_ = [('swap', POINTER(c_long)), ('distance_diff', c_long),
                ('paths_join', c_long), ('child_moves', c_long)]


