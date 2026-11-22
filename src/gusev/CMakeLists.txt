# gusev module: end-linked networks built by Gusev's Monte Carlo procedure, with the strands as
# Brownian bridges between their cross-linkers, as starting structures for LAMMPS. Builds the
# static library `gusev` and the command-line front end `gusev_build`.

add_library(gusev STATIC
  lib/packing.cpp
  lib/linking.cpp
  lib/config.cpp
  lib/run.cpp)
target_include_directories(gusev PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(gusev PUBLIC poly::common poly::chains poly::endlink PRIVATE poly_warnings)
add_library(poly::gusev ALIAS gusev)

add_executable(gusev_build apps/gusev_build.cpp)
target_link_libraries(gusev_build PRIVATE gusev poly_warnings)

install(TARGETS gusev_build RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
