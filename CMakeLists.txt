#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(myProgram)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include  
                     ${Geant4_INCLUDE_DIR} )
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(myProgram main.cc ${sources} ${headers})
target_link_libraries(myProgram ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Enable switching between diamond detector and reference silicon
#
option(USE_SILICON "If ON, the reference silicon is" OFF) # Disabled by default
if(USE_SILICON)
	add_definitions(-DUSING_SILICON)
endif(USE_SILICON)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build basic_radioprotection. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(radioprotection_SCRIPTS
   vis.mac plot_fd.py run.mac
  )

foreach(_script ${radioprotection_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Copy the appropriate GPS macro depending on the setup used
#
if(NOT USE_SILICON)
	configure_file(
		${PROJECT_SOURCE_DIR}/gps_diamond.mac
		${PROJECT_BINARY_DIR}/gps.mac
    		COPYONLY
    	)
	message("**Using diamond detector**")
else()
	configure_file(
		${PROJECT_SOURCE_DIR}/gps_silicon.mac
		${PROJECT_BINARY_DIR}/gps.mac
    		COPYONLY
    	)
	message("**Using silicon detector**")
endif()

#----------------------------------------------------------------------------
# Remove USE_SILICON from the cache, otherwise the next time you run cmake the
# previous value would still be used (comment the following line if you want to
# remember it instead)
unset(USE_SILICON CACHE)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS myProgram DESTINATION bin)


