
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
	fsi_bypass
	fsi_offline_partition
    fsi_interface_mass
    fsi_atheroma_only_thick
    fsi_atheroma_double_layer
    fsi_atheroma_only_thick_varying_young
    fsi_atheroma_double_layer_varying_young
    hdf5ToVTK
)
