add_executable(test_igl_geom test_igl_geom.cc)
target_link_libraries(test_igl_geom
${Boost_LIBRARIES}
)

add_executable(test_lscm test_lscm.cc)
target_link_libraries(test_lscm
${Boost_LIBRARIES}
)

add_executable(test_harmonic test_harmonic.cc)
target_link_libraries(test_harmonic
${Boost_LIBRARIES}
)

add_executable(test_curvature test_curvature.cc)
target_link_libraries(test_curvature
${Boost_LIBRARIES}
)

add_executable(test_cotmatrix test_cotmatrix.cc)

add_executable(test_unit test_unit.cc)
target_link_libraries(test_unit
    riemann
    zjucad-ptree
    ${Boost_LIBRARIES}
    lapack
)

add_executable(test_dcp test_dcp.cc)
target_link_libraries(test_dcp
    riemann
    ${Boost_LIBRARIES}
)

add_executable(test_green_2d test_green_2d.cc)
target_link_libraries(test_green_2d
    riemann
    ${Boost_LIBRARIES}
)

add_executable(test_green_3d test_green_3d.cc)
target_link_libraries(test_green_3d
    riemann
    ${Boost_LIBRARIES}
)

add_executable(test_tri_def_grad test_tri_def_grad.cc)
target_link_libraries(test_tri_def_grad
    jtf-mesh
    ${Boost_LIBRARIES}
    lapack
)

add_executable(test_vec_field_deform test_vec_field_deform.cc)
target_link_libraries(test_vec_field_deform
    riemann
    ${Boost_LIBRARIES}
)

add_executable(test_arap test_arap.cc)
target_link_libraries(test_arap
    riemann
    ${Boost_LIBRARIES}
)

add_executable(test_ff_deform test_ff_deform.cc)
target_link_libraries(test_ff_deform
    riemann
    ${Boost_LIBRARIES}
)

add_executable(test_deform_transfer test_deform_transfer.cc)
target_link_libraries(test_deform_transfer
    riemann
    ${Boost_LIBRARIES}
)

find_package(Arpack REQUIRED)
add_executable(test_lap_eigen test_lap_eigen.cc)
target_link_libraries(test_lap_eigen
    riemann
    ${Boost_LIBRARIES}
    ${ARPACK_LIBRARIES}
)

add_executable(test_grad_edit test_grad_edit.cc)
target_link_libraries(test_grad_edit
    riemann
    ${Boost_LIBRARIES}
)

add_executable(test_cross_param test_cross_param.cc)
target_link_libraries(test_cross_param
    riemann
    ${Boost_LIBRARIES}
)

add_executable(test_vtk_io test_vtk_io.cc)
target_link_libraries(test_vtk_io
    riemann
    ${Boost_LIBRARIES}
)

add_executable(test_wave_construction test_wave_construction.cc)
target_link_libraries(test_wave_construction
    riemann
    ${Boost_LIBRARIES}
)

if(VIENNACL_FOUND)
    add_executable(test_viennacl test_viennacl.cc)
endif(VIENNACL_FOUND)

add_executable(test_shell_interp test_shell_interp.cc)
target_link_libraries(test_shell_interp
    riemann
    ${Boost_LIBRARIES}
)

add_executable(test_bounded_dist test_bounded_dist.cc)
target_link_libraries(test_bounded_dist
    riemann
    ${Boost_LIBRARIES}
)
