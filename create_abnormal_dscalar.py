import nibabel as ni
import numpy as np
import cifti_helpers as ch

from nibabel import cifti2 as ci

out_file = "./data/random.dscalar.nii"

cortex_left_indices = np.loadtxt("./data/CortexLeftIndices.32k_fs_LR.txt", np.integer)
cortex_right_indices = np.loadtxt("./data/CortexRightIndices.32k_fs_LR.txt", np.integer)
wrong_cortex_right_indices = np.random.choice(ch.HCP_32K_FS_R, 
    size=len(cortex_right_indices), replace=False)
wrong_cortex_right_indices.sort()

scalar_info = ch.NamedMapInfo(name="random", meta={})
left_cortex = ch.MapInfo(brain_structure=ch.Structure.CORTEX_LEFT,
    indices=cortex_left_indices, model_type=ch.ModelType.SURFACE,
    surface_number_of_vertices=ch.HCP_32K_FS_L)
right_cortex = ch.MapInfo(brain_structure=ch.Structure.CORTEX_RIGHT,
    indices=wrong_cortex_right_indices, model_type=ch.ModelType.SURFACE,
    surface_number_of_vertices=ch.HCP_32K_FS_R)

scalar_map = ch.create_scalar_map((0,), [ch.NamedMapInfo(name="wrong", meta={})])
brain_models = ch.create_brain_models([left_cortex, right_cortex])
geometry_map = ch.create_geometry_map((1,), brain_models)
data = np.random.uniform(size=len(cortex_left_indices) + len(wrong_cortex_right_indices))


     
