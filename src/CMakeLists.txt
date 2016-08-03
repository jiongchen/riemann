file(GLOB header *.h *.hpp)
file(GLOB source *.cc *.cpp *.c)

enable_language(Fortran)
set(SRC_DIR "{${CMAKE_CURRENT_SOURCE_DIR}}")

add_custom_command(OUTPUT quad_scalar_field.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b quad_scalar_field.mac
    DEPENDS quad_scalar_field.mac
)

add_custom_command(OUTPUT dt_energy.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b dt_energy.mac
    DEPENDS dt_energy.mac
)

add_custom_command(OUTPUT wave_condition.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b wave_condition.mac
    DEPENDS wave_condition.mac
)

add_custom_command(OUTPUT shell.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b shell.mac
    DEPENDS shell.mac
)

add_custom_command(OUTPUT mips.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b mips.mac
    DEPENDS mips.mac
)

add_custom_command(OUTPUT tet_arap.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b tet_arap.mac
    DEPENDS tet_arap.mac
)

add_custom_command(OUTPUT aqp.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b aqp.mac
    DEPENDS aqp.mac
)

add_custom_command(OUTPUT poly_cubic.f90
  COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b poly_cubic.mac
  DEPENDS poly_cubic.mac
)

add_custom_command(OUTPUT l1_frame.f90
  COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b l1_frame.mac
  DEPENDS l1_frame.mac
)

add_custom_command(OUTPUT polycube.f90
  COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b polycube.mac
  DEPENDS polycube.mac
)

set(f90_source quad_scalar_field.f90 dt_energy.f90 wave_condition.f90
               shell.f90 mips.f90 tet_arap.f90 aqp.f90 cubic_sym_sh.f90
               poly_cubic.f90 l1_frame.f90 polycube.f90)

add_library(riemann SHARED ${header} ${source} ${f90_source})
target_link_libraries(riemann
    jtf-mesh
    ${PETSC_LIBRARIES}
    ${UMFPACK_LIBRARIES}
    ${ARPACK_LIBRARIES}
    alglib
    minpack
    lbfgs
)
