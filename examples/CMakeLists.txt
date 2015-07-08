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
    surfparam
    zjucad-ptree
    ${Boost_LIBRARIES}
)

add_executable(test_dcp test_dcp.cc)
target_link_libraries(test_dcp
    surfparam
    ${Boost_LIBRARIES}
)

add_executable(test_green_2d test_green_2d.cc)
target_link_libraries(test_green_2d
    surfparam
    ${Boost_LIBRARIES}
)

add_executable(test_green_3d test_green_3d.cc)
target_link_libraries(test_green_3d
    surfparam
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
    surfparam
    ${Boost_LIBRARIES}
)

add_executable(test_arap test_arap.cc)
target_link_libraries(test_arap
    surfparam
    ${Boost_LIBRARIES}
)

add_executable(test_ff_deform test_ff_deform.cc)
target_link_libraries(test_ff_deform
    surfparam
    ${Boost_LIBRARIES}
)

add_executable(test_deform_transfer test_deform_transfer.cc)
target_link_libraries(test_deform_transfer
    surfparam
    ${Boost_LIBRARIES}
)

find_package(Arpack REQUIRED)
add_executable(test_lap_eigen test_lap_eigen.cc)
target_link_libraries(test_lap_eigen
    surfparam
    ${Boost_LIBRARIES}
    ${ARPACK_LIBRARIES}
)
