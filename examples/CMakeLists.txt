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
