cmake_minimum_required(VERSION 3.10)
project(GravitySimulator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find SFML
find_package(SFML 2.5 COMPONENTS graphics window system QUIET)
if(NOT SFML_FOUND)
    # Fallback: try pkg-config
    find_package(PkgConfig)
    if(PkgConfig_FOUND)
        pkg_check_modules(SFML sfml-graphics sfml-window sfml-system)
    endif()
endif()

# Source files
set(SOURCES
    src/main.cpp
    src/CelestialBody.cpp
    src/GravitySimulator.cpp
    src/Scenarios.cpp
)

add_executable(gravitysim ${SOURCES})

target_include_directories(gravitysim PRIVATE include)

if(SFML_FOUND)
    target_link_libraries(gravitysim sfml-graphics sfml-window sfml-system)
elseif(SFML_LIBRARIES)
    target_include_directories(gravitysim PRIVATE ${SFML_INCLUDE_DIRS})
    target_link_libraries(gravitysim ${SFML_LIBRARIES})
else()
    # Manual linking as last resort
    target_link_libraries(gravitysim sfml-graphics sfml-window sfml-system)
endif()