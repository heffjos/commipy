import numpy as np
import nibabel as ni

from collections import namedtuple

from nibabel import cifti2 as ci

HCP_32K_FS_L = 32492
HCP_32K_FS_R = 32492

# ci.CIFTI_MODEL_TYPES
class ModelType():
    SURFACE = "CIFTI_MODEL_TYPE_SURFACE"
    VOXEL = "CIFTI_MODEL_TYPE_VOXELS"

# ci.CIFTI_BRAIN_STRUCTURES
class Structure():
    ACCUMBENS_LEFT = "CIFTI_STRUCTURE_ACCUMBENS_LEFT"
    ACCUMBENS_RIGHT = "CIFTI_STRUCTURE_ACCUMBENS_RIGHT"
    ALL_WHITE_MATTER = "CIFTI_STRUCTURE_ALL_WHITE_MATTER"
    ALL_GREY_MATTER = "CIFTI_STRUCTURE_ALL_GREY_MATTER"
    AMYGDALA_LEFT = "CIFTI_STRUCTURE_AMYGDALA_LEFT"
    AMYGDALA_RIGHT = "CIFTI_STRUCTURE_AMYGDALA_RIGHT"
    BRAIN_STEM = "CIFTI_STRUCTURE_BRAIN_STEM"
    CAUDATE_LEFT = "CIFTI_STRUCTURE_CAUDATE_LEFT"
    CAUDATE_RIGHT = "CIFTI_STRUCTURE_CAUDATE_RIGHT"
    CEREBELLAR_WHITE_MATTER_LEFT = "CIFTI_STRUCTURE_CEREBELLAR_WHITE_MATTER_LEFT"
    CEREBELLAR_WHITE_MATTER_RIGHT = "CIFTI_STRUCTURE_CEREBELKLAR_WHITE_MATTER_RIGHT"
    CEREBELLUM = "CIFTI_STRUCTURE_CEREBELLUM"
    CEREBELLUM_LEFT = "CIFTI_STRUCTURE_CEREBELLUM_LEFT"
    CEREBELLUM_RIGHT = "CIFTI_STRUCTURE_CEREBELLUM_RIGHT"
    CEREBRAL_WHITE_MATTER_LEFT = "CIFTI_STRUCTURE_CEREBRAL_WHITE_MATTER_LEFT"
    CEREBRAL_WHITE_MATTER_RIGHT = "CIFTI_STRUCTURE_CEREBRAL_WHITE_MATTER_RIGHT"
    CORTEX = "CIFTI_STRUCTURE_CORTEX"
    CORTEX_LEFT = "CIFTI_STRUCTURE_CORTEX_LEFT"
    CORTEX_RIGHT = "CIFTI_STRUCTURE_CORTEX_RIGHT"
    DIENCEPHALON_VENTRAL_LEFT = "CIFTI_STRUCTURE_DIENCEPHALON_VENTRAL_LEFT"
    DIENCEPHALON_VENTRAL_RIGHT = "CIFTI_STRUCTURE_DIENCEPHALON_VENTRAL_RIGHT"
    HIPPOCAMPUS_LEFT = "CIFTI_STRUCTURE_HIPPOCAMPUS_LEFT"
    HIPPOCAMPUS_RIGHT = "CIFTI_STRUCTURE_HIPPOCAMPUS_RIGHT"
    OTHER = "CIFTI_STRUCTURE_OTHER"
    OTHER_GREY_MATTER = "CIFTI_STRUCTURE_OTHER_GREY_MATTER"
    OTHER_WHITE_MATTER = "CIFTI_STRUCTURE_OTHER_WHITE_MATTER"
    PALLIDUM_LEFT = "CIFTI_STRUCTURE_PALLIDUM_LEFT"
    PALLIDUM_RIGHT = "CIFTI_STRUCTURE_PALLIDUM_RIGHT"
    PUTAMEN_LEFT = "CIFTI_STRUCTURE_PUTAMEN_LEFT"
    PUTAMEN_RIGHT = "CIFTI_STRUCTURE_PUTAMEN_RIGHT"
    THALAMUS_LEFT = "CIFTI_STRUCTURE_THALAMUS_LEFT"
    THALAMUS_RIGHT = "CIFTI_STRUCTURE_THALAMUS_RIGHT"

# ci.CIFTI_MAP_TYPES
class Map():
    BRAIN_MODELS = "CIFTI_INDEX_TYPE_BRAIN_MODELS"
    PARCELS = "CIFTI_INDEX_TYPE_PARCELS"
    SERIES = "CIFTI_INDEX_TYPE_SERIES"
    SCALARS = "CIFTI_INDEX_TYPE_SCALARS"
    LABELS = "CIFTI_INDEX_TYPE_LABELS"

MapInfo = namedtuple("MapInfo",
    ["brain_structure", "indices", "model_type", "surface_number_of_vertices"])
NamedMapInfo = namedtuple("NamedMapInfo",
    ["name", "meta"])

def create_brain_models(info):
    """Creates a list of brain_models from a list of MapInfo"""
    offset = 0
    brain_models = []
    structures = []
    for i in info:
        if i.brain_structure in structures:
            raise Exception
        else:
            structures.append(i.brain_structure)
            
        if i.model_type == ModelType.SURFACE:
            vertices = ci.Cifti2VertexIndices(i.indices)
            brain_models.append(ci.Cifti2BrainModel(
                index_offset=offset, index_count=len(vertices), 
                model_type=i.model_type, brain_structure=i.brain_structure,
                vertex_indices=vertices, n_surface_vertices=i.surface_number_of_vertices
            ))
            offset += len(vertices)
        else:
            voxels = ci.Cifti2VoxelIndicesIJK(i.indices)
            brain_models.append(ci.Cifti2BrainModel(
                index_offset=offset, index_count=len(voxels),
                model_type=i.model_type, brain_structure=i.brain_structure,
                voxel_indices_ijk=voxels
            ))
            offset += len(voxels)
    return brain_models

def create_volume(dims, affine, exp):
    """Creates a cifti volume"""
    return ci.Cifti2Volume(dims,
        ci.Cifti2TransformationMatrixVoxelIndicesIJKtoXYZ(exp, affine)
    )
        
def create_geometry_map(applies_to_matrix_dimension, brain_models):
    """Creates a geometry map"""
    return ci.Cifti2MatrixIndicesMap(applies_to_matrix_dimension,
        Map.BRAIN_MODELS, maps=brain_models)

def create_scalar_map(applies_to_matrix_dimension, info):
    """Creates a scalar map form a list of NamedMapInfo"""
    maps = [ci.Cifti2NamedMap(i.name, ci.Cifti2MetaData(i.meta))
            for i in info]
    return ci.Cifti2MatrixIndicesMap(applies_to_matrix_dimension,
        Map.SCALARS, maps=maps)

def create_dscalar(scalar_map, geometry_map, data):
    matrix = ci.Cifti2Matrix()
    matrix.append(scalar_map)
    matrix.append(geometry_map)
    hdr = ci.Cifti2Header(matrix)
    img = ci.Cifti2Image(data, hdr)
    img.nifti_header.set_intent("NIFTI_INTENT_CONNECTIVITY_DENSE_SCALARS")
    return img

    
def equivalent_brain_models(c1, c2):
    # bm1 = c1.matrix.header.
    pass

"""
check what happens in cifti-math when data matrix has same dimensions, but
different vertices. To do this, I need to create a dense scalar with 
different selected vertices

cifti-math throws an error stating the cifti files do not have same brainordinates
"""

"""
check to see if remapping labels with wb_command is easier
"""

"""
correlation ideas:
surface roi in fs_LR space
csv file listing participant resting functional cifti
output prefix
flag for writing individual connectivity maps
"""
