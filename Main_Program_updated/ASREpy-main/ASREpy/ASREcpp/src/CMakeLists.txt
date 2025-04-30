# Define minimum cmake version
cmake_minimum_required(VERSION 3.12)

# Define project name
project(ASRElib VERSION 1.0)

# Declare libraries
add_library(ASRElibTimoBeam SHARED ASRElibTimoBeam.cpp)
add_library(ASRElib2DFrame SHARED ASRElib2DFrame.cpp)
add_library(ASRElibTimoBeam_ps SHARED ASRElibTimoBeam_ps.cpp)

# Include the headers in inc
target_include_directories(ASRElibTimoBeam PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/inc")
target_include_directories(ASRElib2DFrame PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/inc")
target_include_directories(ASRElibTimoBeam_ps PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/inc")

# Create an INTERFACE target for the header-only library Eigen
add_library(Eigen INTERFACE)

# Specify the include directories for the header files of Eigen
target_include_directories(Eigen INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/Eigen")

# Link Eigen library to ASRElib
target_link_libraries(ASRElibTimoBeam PRIVATE Eigen)
target_link_libraries(ASRElib2DFrame PRIVATE Eigen)
target_link_libraries(ASRElibTimoBeam_ps PRIVATE Eigen)

# Get parent directory
get_filename_component(parent_dir "${CMAKE_CURRENT_SOURCE_DIR}" DIRECTORY)

# If Windows add a MACRO _MSVC
IF(MSVC AND "${CMAKE_VS_PLATFORM_NAME}" MATCHES "(x64)")
    message("   ************************************")
    message("   Compiling on x64 Windows using MSVC")
    message("   ************************************")
    target_compile_definitions(ASRElibTimoBeam PRIVATE _MSVC)
    target_compile_definitions(ASRElib2DFrame PRIVATE _MSVC)
    target_compile_definitions(ASRElibTimoBeam_ps PRIVATE _MSVC)
    # Define the dir to install the compiled dll
    install(TARGETS ASRElibTimoBeam DESTINATION "${parent_dir}/bin/win32")
    install(TARGETS ASRElib2DFrame DESTINATION "${parent_dir}/bin/win32")
    install(TARGETS ASRElibTimoBeam_ps DESTINATION "${parent_dir}/bin/win32")
ELSEIF(MSVC)
    message("   ************************************")
    message("   Compiling on ${CMAKE_VS_PLATFORM_NAME} Windows using MSVC")
    message("   ************************************")
    target_compile_definitions(ASRElibTimoBeam PRIVATE _MSVC)
    target_compile_definitions(ASRElib2DFrame PRIVATE _MSVC)
    target_compile_definitions(ASRElibTimoBeam_ps PRIVATE _MSVC)
    install(TARGETS ASRElibTimoBeam DESTINATION "${parent_dir}/bin/win32")
    install(TARGETS ASRElib2DFrame DESTINATION "${parent_dir}/bin/win32")
    install(TARGETS ASRElibTimoBeam_ps DESTINATION "${parent_dir}/bin/win32")
ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message("   ************************************")
    message("   Compiling using GNU")
    message("   ************************************")
    target_compile_definitions(ASRElibTimoBeam PRIVATE _GCC)
    target_compile_definitions(ASRElib2DFrame PRIVATE _GCC)
    target_compile_definitions(ASRElibTimoBeam_ps PRIVATE _GCC)
    IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message(STATUS "Running on macOS")
        install(TARGETS ASRElibTimoBeam DESTINATION "${parent_dir}/bin/macOS")
        install(TARGETS ASRElib2DFrame DESTINATION "${parent_dir}/bin/macOS")
        install(TARGETS ASRElibTimoBeam_ps DESTINATION "${parent_dir}/bin/macOS")
    ELSEIF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        message(STATUS "Running on Linux")
        install(TARGETS ASRElibTimoBeam DESTINATION "${parent_dir}/bin/unix")
        install(TARGETS ASRElib2DFrame DESTINATION "${parent_dir}/bin/unix")
        install(TARGETS ASRElibTimoBeam_ps DESTINATION "${parent_dir}/bin/unix")
    ENDIF()
ELSE(TRUE)
    message("   ************************************")
    message("   Compiling on non-MSVC nor GNU")
    message("   ************************************")
ENDIF(MSVC AND "${CMAKE_VS_PLATFORM_NAME}" MATCHES "(x64)")
