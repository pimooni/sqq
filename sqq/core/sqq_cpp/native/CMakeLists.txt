cmake_minimum_required(VERSION 3.20)

if(DEFINED SKBUILD_PROJECT_VERSION)
    set(SQQ_PROJECT_VERSION "${SKBUILD_PROJECT_VERSION}")
else()
    file(STRINGS
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../_version.py"
        SQQ_VERSION_LINE
        REGEX "^__version__[ \\t]*="
    )
    string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" SQQ_PROJECT_VERSION "${SQQ_VERSION_LINE}")
endif()

if(NOT SQQ_PROJECT_VERSION)
    message(FATAL_ERROR "Unable to determine the SQQ version")
endif()

project(sqq_cpp_core VERSION ${SQQ_PROJECT_VERSION} LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_sqq_cpp MODULE
    src/core.cpp
    src/python_bindings.cpp
)

target_include_directories(_sqq_cpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(_sqq_cpp PRIVATE cxx_std_17)
target_compile_definitions(_sqq_cpp PRIVATE SQQ_VERSION="${SQQ_PROJECT_VERSION}")
set_target_properties(_sqq_cpp PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
)

if(MSVC)
    target_compile_options(_sqq_cpp PRIVATE /O2 /EHsc /permissive- /bigobj /fp:precise)
else()
    target_compile_options(_sqq_cpp PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

install(TARGETS _sqq_cpp
    LIBRARY DESTINATION sqq/core/sqq_cpp
    RUNTIME DESTINATION sqq/core/sqq_cpp
)
