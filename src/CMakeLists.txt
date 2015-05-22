file(GLOB header *.h)
file(GLOB source *.cc)

enable_language(Fortran)
set(SRC_DIR "{${CMAKE_CURRENT_SOURCE_DIR},$ENV{HOME}/.maxima}")

add_custom_command(OUTPUT quad_scalar_field.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b quad_scalar_field.mac
    DEPENDS quad_scalar_field.mac
)

add_library(surfparam SHARED ${header} ${source} quad_scalar_field.f90)
target_link_libraries(surfparam
    jtf-mesh
)
