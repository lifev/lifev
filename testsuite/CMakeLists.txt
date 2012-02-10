set (FOLDER_LIST test_MatrixEpetraStructuredFramework test_TemplateTest lifearray lifecore test_Parser test_VectorContainer lifefem test_FEFunction lifefilters lifemesh test_navierStokes test_ADRAssembler test_bdf test_darcy test_fsi test_hyperbolic test_LevelSetSolver test_structure test_partition test_timeAdvance_orderI test_timeAdvance_orderII test_monolithic)
set (TEST_LIST test_MatrixEpetraStructuredFramework test_TemplateTest test_meshentitycontainer test_vectorsmall test_getPot test_Parser test_VectorContainer test_currentFE test_quadRule test_interpolate test_FEFunction test_filters testExport testImport test_mesh test_navierStokes test_navierStokes_2d test_ADRAssembler test_ADRAssembler_2d test_bdf test_bdf_2d test_darcy test_fsi test_hyperbolic test_LevelSetSolver test_structure test_partition test_cylinder test_timeAdvance1 test_timeAdvance2 test_monolithic)

foreach (component ${FOLDER_LIST})
    add_subdirectory(${component})
endforeach()

# Go down into the data subdirectory too
add_subdirectory (data)

add_custom_target (all_tests
                   DEPENDS ${TEST_LIST})
