add_executable(exampleMuscleSplines
    mainMuscleSplines.cpp)
target_link_libraries(exampleMuscleSplines ${Simbody_LIBRARIES} osimTools)
set_target_properties(exampleMuscleSplines PROPERTIES FOLDER "Examples")

# The answers file uses `#pragma region` to facilitate code folding
# when using this example in a live demonstration, but these pragmas are only
# recognized by MSVC.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
        CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(exampleHopperDeviceAnswers PUBLIC
            -Wno-unknown-pragmas)
endif()
