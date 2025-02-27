project(ExampleMyCable)

# find_package(OpenSim REQUIRED PATHS "${OPENSIM_INSTALL_DIR}")

# target_link_libraries(${TARGET} ${OpenSim_LIBRARIES})

add_executable(exampleMyCable
    mainMyCable.cpp)
target_link_libraries(exampleMyCable ${Simbody_LIBRARIES} osimTools osimExampleComponents)
set_target_properties(exampleMyCable PROPERTIES FOLDER "Examples")
