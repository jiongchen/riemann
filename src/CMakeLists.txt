file(GLOB header *.h)
file(GLOB source *.cc)

add_library(surfparam SHARED ${header} ${source})
target_link_libraries(surfparam
    jtf-mesh
)
