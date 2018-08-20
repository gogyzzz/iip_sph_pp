cmake_minimum_required(VERSION 3.4)
project(iip_sph_pp)

# TARGET_NAME will be the name of executable.
set(TARGET_NAME ${PROJECT_NAME})

#### BUILD ENVIRONMENT OPTIONS ######
option(USE_BLAS    "Using BLAS"      ON)
option(USE_OPENMP  "Using OpenMP"    ON)

# Not supported, yet.
#option(USE_CUDA    "Using CUDA"      OFF)

#################################
#         CMAKE INDEX           #
#  1. Macros                    #
#  2. Sources and an executable #
#  3. Header Includes           #
#  4. Check options             #
#  5. Defines and Libraries     #
#  6. Compile optinos           #
#  7. WIP : CUDA                #
#################################

# setting options by 'set' is unstable. Often not working.
# recommed using 'target_<blah>_<blah>'
# which is more robust and modern.

#### 1. Macros ####

macro(display_elements
	DISPLAY_MESSAGE
	LIST_TO_DISPLAY)
    message("[ ${DISPLAY_MESSAGE} ]")
    foreach(_DISPLAY_ELEMENT ${LIST_TO_DISPLAY})
	message(STATUS ${_DISPLAY_ELEMENT})
    endforeach()
endmacro()
#### 2. sources and executable ####

set(MAIN_SRC 
	# Add your source files here.
	# including the source which contains main().
	# root directory is location of CMakeLists.txt 
	# EXAMPLE
	# test/test_matlab.c
	
	)

#library source files
set(C_SOURCE
    source/iip_matrix.c
    source/iip_wav.c
    source/iip_blas_lv1.c
    source/iip_blas_lv2.c
    source/iip_blas_lv3.c
    source/iip_math.c
    source/iip_memory.c
    source/iip_time.c
    source/iip_io.c
    source/iip_invert.c
    source/iip_test.c
    )

add_executable(${TARGET_NAME} ${MAIN_SRC} ${C_SOURCE})

#### 3. Header includes ####

target_include_directories(${TARGET_NAME}
	PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/header
	# Add additional header folders here.
	#${CMAKE_CURRENT_SOURCE_DIR} means root directory.
	)

#### 4. check options ####

if(USE_BLAS)
# MKL also can be found by find_package(BLAS)
# But, MKL requires specific compile & link option.
# just using "target_link_libraries" doesn't work for MKL.
  set(BLA_VENDOR "OpenBLAS")
    find_package(BLAS)
	if(BLAS_FOUND)
    message(STATUS "Found OpenBLAS : " ${BLAS_LIBRARIES})
	set(USE_OPEN ON)
	set(USE_MKL  OFF)
    else()
		message(STATUS "OpenBLAS not found, assume usnig MKL")
	set(USE_OPEN OFF)
	set(USE_MKL  ON)
    endif()
endif()

if(USE_OPENMP)
   find_package(OpenMP)
   if(NOT(OpenMP_FOUND))
	    message(SEND_ERROR "USE_OPENMP is set ON but OpenMP not found")
	    set(USE_OPENMP OFF)
   else()
	  message(STATUS "OpenMP found")
   endif()
endif()

message(STATUS "USE_CUDA   : " ${USE_CUDA})
message(STATUS "USE_OPEN  : " ${USE_OPEN})
message(STATUS "USE_MKL   : " ${USE_MKL})
message(STATUS "USE_OPENMP : " ${USE_OPENMP})

#### 5. defines and links  ####

#compile flags
set(CF)
#link flags
set(LF)
#compule definitions
set(CD)

# Addition flags and defines can be attached by
# list(APPEND $(FLAG or Defines))
# EXAMPLE
# list(APPEND CF /W3) -> add /W3 to compile option
# list(APPEND LF -pthread) -> add -pthread as link option
# list(APPEND CD my_define 100) ->add #define my_define 100 at pre-process time

#GET OS
#UNIX   TRUE for UNIX,APPLE
#WIN32  TRUE for Windows 32bit and 64bit
#APPLE  TRUE for APPLE
if(UNIX)
  list(APPEND CD
    OS_UNIX=1
    )
elseif(WIN32)
  list(APPEND CD  
    OS_WIN=1
    )
endif()

# Check whether PATH is set correctly or not.
# message(STATUS "ENV_PATH : " $ENV{PATH})

if(NOT MSVC)
  list(APPEND LF m)
endif(NOT MSVC)


if(USE_OPEN OR USE_MKL )

  list(APPEND CD
    USE_CBLAS=1
    )

 if(USE_OPEN)
	message(STATUS ${BLAS_LIBRARIES})
  list(APPEND LF
    openblas
    )
   
  list(APPEND CD
	    USE_OPEN=1
      )
  elseif(USE_MKL)

  list(APPEND CD  
	    USE_MKL=1
      )
		
  endif()
endif(USE_OPEN OR USE_MKL)

if(USE_OPENMP)

if(MSVC)
   # Visual Studio uses '/' for compile option
    list(APPEND CF
    /openmp
    )
 # MSVC works well with just compile option /openmp 

else(MSVC)
  list(APPEND CF 
   -fopenmp
   )
  list(APPEND LF  
   -fopenmp
   )
endif(MSVC)

list(APPEND CD    
     USE_OMP=1
    )
endif(USE_OPENMP)


##### 6. compile optinos ####

if(WIN32 AND MSVC)
  # No /O3 in VS2013 (higher, too)?

  # /O2 - fastest code in the majority of cases
  # /MP - compiled in mult-process : faster build time
  #list(APPEND  CMAKE_C_FLAGS " /O2 /MP ") -- This method is unstable.
  list(APPEND CF
     /MP
     # /O2 - not compatible with debug mod flag /RTC1
     # you have to set manually /O2 in release mod
     )
 else(WIN32 AND MSVC)
  # list(APPEND CMAKE_C_FLAGS "-g -O3 -Wall -Wextra ") : Wall is too much for now
  list(APPEND CF
      -g
      -O3
      -Wextra
      )
  endif(WIN32 AND MSVC)

target_compile_definitions(${TARGET_NAME} PRIVATE ${CD} )
target_compile_options    (${TARGET_NAME} PRIVATE ${CF} )
target_link_libraries     (${TARGET_NAME} PRIVATE ${LF} )

display_elements("MAIN SOURCE" "${MAIN_SRC}")
display_elements("LIBRARY SOURCE" "${C_SOURCE}")
display_elements("Compile Flags" "${CF}")
display_elements("Linker Flags" "${LF}")
display_elements("Compile Defines" "${CD}")

message(STATUS "CMAKE_GENERATOR : "  ${CMAKE_GENERATOR})

#### 7. misc : CUDA(Work in Progress) ####

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

if(USE_CUDA)
    find_package(CUDA)
    if(NOT(CUDA_FOUND))
		message(SEND_ERROR "USE_CUDA is set ON but CUDA not found")
		set(USE_CUDA OFF)
    else()
	if(NOT CUDA_ARCH)
		message(STATUS "cuda arch wasn't provided, running feature test")

    include(feature_test)
	else()
	    set(CUDA_GEN_CODE "-gencode arch=compute_${CUDA_ARCH},code=sm_${CUDA_ARCH}")
		endif()
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}; ${CUDA_GEN_CODE}")

	message(STATUS "CUDA found")
    endif()
endif()


if(USE_CUDA)
    set(CUDA_TARGET_NAME "${PROJECT_NAME}_cuda")
    set(CUDA_SRC 
	${PROJECT_SOURCE_DIR}/source/iip_matrix.cu
	${PROJECT_SOURCE_DIR}/source/iip_wav.cu
	${PROJECT_SOURCE_DIR}/source/iip_blas_lv1.cu
	${PROJECT_SOURCE_DIR}/source/iip_blas_lv2.cu
	${PROJECT_SOURCE_DIR}/source/iip_blas_lv3.cu)
	
    display_elements("Source files for CUDA_SRC" "${CUDA_SRC}")
    # compile .c files by nvcc
    #SET_SOURCE_FILES_PROPERTIES(${COMMON_SRC} PROPERTIES LANGUAGE CUDA)

    if(WIN32)
	enable_language(CXX)
	SET_SOURCE_FILES_PROPERTIES(${MAIN_SRC} PROPERTIES LANGUAGE CXX)	
    else(WIN32)
	SET_SOURCE_FILES_PROPERTIES(${MAIN_SRC} PROPERTIES LANGUAGE CUDA)
    endif(WIN32)

    cuda_add_executable(${CUDA_TARGET_NAME} ${MAIN_SRC} ${CUDA_SRC})		
    cuda_add_cublas_to_target(${CUDA_TARGET_NAME})
    target_include_directories(${CUDA_TARGET_NAME}
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/header)
    target_link_libraries(${CUDA_TARGET_NAME}
	${CUDA_LIBRARIES})
    target_compile_definitions(${CUDA_TARGET_NAME}
	PRIVATE
	USE_CUDA=1)
endif()

