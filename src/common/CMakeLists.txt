# Code shared by every module.
#   poly::common     header-only: random number generator, configuration-file reader, index set,
#                    geometry (vectors, FCC lattice), topology (molecules, bond lists) and output
#                    (files, PDB, PSF and LAMMPS data writers); see README.md
#   poly::chemistry  radical polymerization chemistry (unit states, reaction schemes, bookkeeping
#                    and per-sample tables), used by the reactive lattice modules
#   poly::chains     chain-length distributions, used by the modules that build chains
#   poly::endlink    end-linking systems of the packing modules (randomwalk, gusev): the system,
#                    its LAMMPS data, XYZ and link writers, the shared keys and run steps

add_library(poly_common INTERFACE)
target_include_directories(poly_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(poly::common ALIAS poly_common)

add_library(poly_chemistry STATIC
  lib/scheme.cpp
  lib/tables.cpp)
target_link_libraries(poly_chemistry PUBLIC poly::common PRIVATE poly_warnings)
add_library(poly::chemistry ALIAS poly_chemistry)

add_library(poly_chains STATIC
  lib/lengths.cpp)
target_link_libraries(poly_chains PUBLIC poly::common PRIVATE poly_warnings)
add_library(poly::chains ALIAS poly_chains)

add_library(poly_endlink STATIC
  lib/endlink_system.cpp
  lib/endlink_writers.cpp
  lib/endlink_config.cpp
  lib/endlink_driver.cpp)
target_link_libraries(poly_endlink PUBLIC poly::common poly::chains PRIVATE poly_warnings)
add_library(poly::endlink ALIAS poly_endlink)

install(DIRECTORY include/ DESTINATION include)
