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

add_executable(test_lap_eigen test_lap_eigen.cc)
target_link_libraries(test_lap_eigen
  riemann
  ${Boost_LIBRARIES}
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

add_executable(boundvert boundvert.cc)
target_link_libraries(boundvert
  jtf-mesh
  )

add_executable(test_spin_trans test_spin_trans.cc)
target_link_libraries(test_spin_trans
  riemann
  ${Boost_LIBRARIES}
  )

add_executable(test_local_edit test_local_edit.cc)
target_link_libraries(test_local_edit
  riemann
  ${Boost_LIBRARIES}
  )

add_executable(test_conformal_volume test_conformal_volume.cc)
target_link_libraries(test_conformal_volume
  riemann
  ${Boost_LIBRARIES}
  )

add_executable(test_vol_frame test_vol_frame.cc)
target_link_libraries(test_vol_frame
  riemann
  ${Boost_LIBRARIES}
  zjucad-ptree
  )

add_executable(test_polycube test_polycube.cc)
target_link_libraries(test_polycube
  riemann
  ${Boost_LIBRARIES}
  zjucad-ptree
  lapack
  )
