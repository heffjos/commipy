import numpy as np
import nibabel as ni

from collections import namedtuple

SurfaceMapInfo = namedtuple("SurfaceMapInfo", 
    ["vertices", "brain_structure", "num_vertices", "type"])
VolumeMapInfo = namedtuple("SurfaceMapInfo",
    ["voxels", "brain_structure", "dimensions", "affine", "meter_exp", "type"])

def create_geometry_map(applies_to_matrix_dimension, info):
    brain_models = []
    
    
def create_dscalar_from_template():
    pass
