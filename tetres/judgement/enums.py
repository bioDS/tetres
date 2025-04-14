from enum import Enum

# todo does this go here or should it be in the higher level, also implement enums for the summary package i.e. centroid algorithm
# all enums used as arguments in this subpackage

class MDS(Enum):
    TSNE = 'tsne'
    PYTHON = 'py'
    R = 'r'
