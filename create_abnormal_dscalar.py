import nibabel as ni
import numpy as np
import cifti_helpers as ch

from nibabel import cifti2 as ci

out_wrong_file = "./data/random.dscalar.nii"
out_correct_file = "./data/hcp.dscalar.nii"

voxels = [[0, 1, 2],
          [1, 2, 3],
          [2, 2, 3],
          [3, 2, 3]]

affine = [[2, 0, 0, -50],
          [0, 2, 0, -50],
          [0, 0, 2, -50],
          [0, 0, 0, 1]]

cortex_left_indices = np.loadtxt("./data/CortexLeftIndices.32k_fs_LR.txt", np.integer)
cortex_right_indices = np.loadtxt("./data/CortexRightIndices.32k_fs_LR.txt", np.integer)
wrong_cortex_right_indices = np.random.choice(ch.HCP_32K_FS_R, 
    size=len(cortex_right_indices), replace=False)
wrong_cortex_right_indices.sort()

right_cortex = ch.MapInfo(brain_structure=ch.Structure.CORTEX_RIGHT,
    indices=cortex_right_indices, model_type=ch.ModelType.SURFACE,
    surface_number_of_vertices=ch.HCP_32K_FS_R)
wrong_right_cortex = ch.MapInfo(brain_structure=ch.Structure.CORTEX_RIGHT,
    indices=wrong_cortex_right_indices, model_type=ch.ModelType.SURFACE,
    surface_number_of_vertices=ch.HCP_32K_FS_R)
left_cortex = ch.MapInfo(brain_structure=ch.Structure.CORTEX_LEFT,
    indices=cortex_left_indices, model_type=ch.ModelType.SURFACE,
    surface_number_of_vertices=ch.HCP_32K_FS_L)
right_thalamus = ch.MapInfo(brain_structure=ch.Structure.THALAMUS_LEFT,
    indices=voxels, model_type=ch.ModelType.VOXEL, 
    surface_number_of_vertices=None)
volume = ch.create_volume([50, 50, 50], affine, -3)

wrong_scalar_map = ch.create_scalar_map((0,), 
    [ch.NamedMapInfo(name="wrong", meta={})])
correct_scalar_map = ch.create_scalar_map((0,), 
    [ch.NamedMapInfo(name="right", meta={})])

wrong_brain_models = ch.create_brain_models([wrong_right_cortex, left_cortex, 
    right_thalamus])
wrong_brain_models.append(volume)
correct_brain_models = ch.create_brain_models([right_cortex, left_cortex, 
    right_thalamus])
correct_brain_models.append(volume)

wrong_geometry_map = ch.create_geometry_map((1,), wrong_brain_models)
correct_geometry_map = ch.create_geometry_map((1, ), correct_brain_models)

data = np.random.uniform(
    size=(1, len(cortex_left_indices) + len(wrong_cortex_right_indices))
)
correct_dscalar_img = ch.create_dscalar(
    correct_scalar_map, correct_geometry_map, data
)
wrong_dscalar_img = ch.create_dscalar(
    wrong_scalar_map, wrong_geometry_map, data
)

check0, explanation0 = ch.approximately_equal_brain_models(
    correct_brain_models[2], correct_brain_models[2])
check1, explanation1 = ch.approximately_equal_brain_models(
    correct_brain_models[0], correct_brain_models[0])
check2, explanation2 = ch.approximately_equal_brain_models(
    correct_brain_models[0], correct_brain_models[1])
check3, explanation3 = ch.approximately_equal_brain_models(
    correct_brain_models[0], wrong_brain_models[0])

# ci.save(correct_dscalar_img, out_correct_file)
# ci.save(wrong_dscalar_img, out_wrong_file)

# surface and volume structures cannot be repeated in cifti file
