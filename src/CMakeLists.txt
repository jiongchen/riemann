file(GLOB header *.h *.hpp)
file(GLOB source *.cc)

enable_language(Fortran)
set(SRC_DIR "{${CMAKE_CURRENT_SOURCE_DIR},$ENV{HOME}/.maxima}")

add_custom_command(OUTPUT quad_scalar_field.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b quad_scalar_field.mac
    DEPENDS quad_scalar_field.mac
)

add_custom_command(OUTPUT dt_energy.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b dt_energy.mac
    DEPENDS dt_energy.mac
)

add_library(riemann SHARED ${header} ${source} quad_scalar_field.f90 dt_energy.f90)
target_link_libraries(riemann
    jtf-mesh
    ${UMFPACK_LIBRARIES}
)
