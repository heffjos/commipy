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
ParcelInfo = namedtuple("ParcelInfo",
    ["name", "surfaces", "voxel_ijk"]) 
SurfaceInfo = namedtuple("SurfaceInfo",
    ["brain_structure", "indices"])

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

def create_series_map(applies_to_matrix_dimension, info):
    """Creates a series map from a list of SeriesInfo"""
    # It is simple enougth to directly use Cifti2MatrixIndicesMap
    pass

def create_img(maps, data):
    matrix = ci.Cifti2Matrix()
    matrix.extend(maps)
    hdr = ci.Cifti2Header(matrix)
    img = ci.Cifti2Image(data, hdr)
    return img

def create_dscalar(maps, data):
    img = create_img(maps, data)
    img.nifti_header.set_intent("NIFTI_INTENT_CONNECTIVITY_DENSE_SCALARS")
    return img

def create_dscalar_from_template(template, data, map_names):
    """create a dscalar with data and map_name from a template dscalar cifti"""
    geometry_map = template.header.matrix.get_index_map(1)
    scalar_map = create_scalar_map((0,), [NamedMapInfo(x[0], x[1])
        for x in map_names])
    return create_dscalar([scalar_map, geometry_map], data)

def create_pscalar():
    img = create_img(maps, data)
    img.nifti_header.set_intent("NIFTI_INTENT_CONNECTIVITY_PARCELLATED_SCALAR")
    return img

def create_parcel(info):
    """Create a list of Cifti2Parcels from a list of ParcelInfo"""
    allowed_brain_structures = (
        Structure.CORTEX_LEFT,
        Structure.CORTEX_RIGHT,
        Structure.CEREBELLUM
    )
    surfaces = []
    volume = None

    has_surfaces = info.surfaces is not None 
    has_ijk = info.voxel_ijk is not None

    if not has_vertices and not has_ijk:
        raise Exception
    if has_vertices:
        for surface in info.surfaces:
            if info.brain_structure not in allowed_brain_structure:
                raise Exception
            surfaces.append(
                ci.Cifti2Vertices(surface.brain_structure, surface.indices))
    if has_indices:
        volume = ci.Cifti2VoxelIndicesIJK(info.voxel_ijk)

    return ci.Cifti2Parcel(info.name, volume, surfaces)

def create_parcel_map(applies_to_matrix_dimension, surfaces, volume, pinfo):
    """Creaters a pracel map

    PARAMETERS
    ---------
    applies_to_matrix_dimension : tuple
    surfaces                    : list of Cifti2Surface
    volume                      : Cifti2Volume
    pinfo                       : list of ParcelInfo
    """
    mapping = ci.Cifti2MatrixIndicesMap(applies_to_matrix_dimension,
        Map.PARCELS)
    for p in pinfo:
        mapping.append(create_parcel(p))
    mapping.extend(surfaces)
    mapping.volume = volume
    return mapping

def create_pscalar_from_dlabel(dlabel, data, name_infos):
    """Create a pscalar from a dlabel with data data"""
    pinfo = []
    d_cii = ci.load(dlabel)
    labels = d_cii.header.matrix.get_index_map(0)
    brain_models = d_cii.header.matrix.get_index_map(1)

    for element in brain_models:
        if isinstance(element, ci.Cifti2BrainModel):
            if element.model_type == ModelType.SURFACE:
                pinfo.append(ParcelInfo(
                    name=element.brain_structure,
                    surfaces=SurfaceInfo(element.brain_structure, element.vertex_indices),
                    voxel_ijk=None
                ))
            if element.model_type == ModelType.VOXEL:
                pinfo.append(ParcelInfo(
                    name=element.brain_structure,
                    surfaces=None,
                    voxel_ijk=voxel_indices_ijk
                ))

def create_dtseries(maps, data):
    img = create_img(maps, data)
    img.nifti_header.set_intent("NIFTI_INTENT_CONNECTIVITY_DENSE_SERIES")
    return img

def equal_voxel_indices(vi1, vi2):
    """Checks equality of two Cifti2VoxelIndicesIJK

    Probably a good idea to implement this as __eq__ in class
    """
    if len(vi1) != len(vi2):
        return False
    for r1, r2 in zip(vi1, vi2):
        for c1, c2 in zip(r1, r2):
            if c1 != c2:
                return False
    return True

def equal_volume_space(v1, v2):
    """Checks equality of Cifti2Volume

    Probably a good idea to implement this as __eq__ in class
    """
    t1 = v1.volume.transformation_matrix_voxel_indices_ijk_to_xyz
    t2 = v2.volume.transformation_matrix_voxel_indices_ijk_to_xyz
    return (v1.volume_dimensions != v2.volume_dimensions 
        and np.all(t1.matrix == t2.matrix) 
        and t1.meter_exponent == t2.meter_exponent)

def equal_vertex_indices(vi1, vi2):
    """Check equality of two Cifti2VertexIndices

    Probably a good idea to implement this as __eq__ in class
    """
    if len(vi1) != len(vi2):
        return False
    for i, j in zip(vi1, vi2):
        if i != j:
            return False
    return True
    
def approximately_equal_brain_models(bm1, bm2):
    """ Checks if two brain models are equivalent.

    Equivalent is the follwoing:
    same model type
    both have volume data
    same vertices/voxels
    same brain structure

    This should be integrated into the Cifti2BrainModel class eventually
    """
    if bm1.model_type != bm2.model_type:
        return (False, "Model types do not match: {}, {}".format(
            bm1.model_type, bm2.model_type))
    if bm1.brain_structure != bm2.brain_structure:
        return (False, "brain structures do not match: {}, {}".format(
            bm1.brain_structure, bm2.brain_structure))
    if (bm1.voxel_indices_ijk is None) != (bm2.voxel_indices_ijk is None):
        return (False, "one of the brain models has no volume data")
    if (bm1.voxel_indices_ijk 
        and not equal_voxel_indices(bm1.voxel_indices_ijk, bm2.voxel_indices_ijk)):
        return (False, "mappings have a different volume space")
    if (bm1.vertex_indices 
        and (bm1.surface_number_of_vertices != bm2.surface_number_of_vertices
        or not equal_vertex_indices(bm1.vertex_indices, bm2.vertex_indices))):
        return (False, "mappings include different brainordinates")
    return (True, "")

def approximately_equal_indices_map(m1, m2):
    if m1.indices_map_to_data_type != m2.indices_map_to_data_type:
        return (False, "unequal number of indices maps")
    if (m1.volume is None) != (m2.volume is None):
        return (False, "one of the indices maps has no volume data")
    
    if (m1.indices_map_to_data_type == Map.BRAIN_MODELS
        or m1.indices_map_to_data_type == Map.SERIES):
        if len(m1) != len(m2):
            return (False, "unequal number of brain structures.")

        sm1 = sorted(m1, key=lambda x : x.brain_structure)
        sm2 = sorted(m2, key=lambda x : x.brain_structure)
        for x, y in zip(sm1, sm2):
            are_equal, exp = approximately_equal_brain_models(x, y)
            if not are_equal:
                return (False, "Structures 1,2: {},{}\n".
                    format(x.brain_structure, y.brain_structure) + exp)

        if m1.volume and not equal_volume_space(m1.volume, m2.volume):
            return (False, "brain models in different volume space")
    elif m1.indices_map_to_data_type == Map.PARCELS:
        raise Exception("Parcels equality not implemented yet.")
    elif m1.indcies_map_to_data_type == Map.SCALARS:
        raise Exception("Scalar equality not implemented yet.")
    elif m1.indices_map_to_data_type == Map.LABELS:
        raise Exception("Labels equality not implemented yet.")
    else:
        raise Exception("Unknown map type: {}.".format(m1.indices_map_to_data_type))
    return (True, "")
    

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
