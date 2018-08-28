# Development Environment.
# CMAKE version of Ubuntu 16.04 is 3.5.
cmake_minimum_required(VERSION 3.5)
# This will be Name of Executable.
project(iip_sph_pp)

#################################
#         CMAKE INDEX           #
#  1. Declarations              #
#  2. ** User Setting **        #
#  3. Macros                    #
#  4. Add executable            #
#  5. Header Includes           #
#  6. Availability of Options   #
#  7. Library Setting           #
#  8. Basic Compile options     #
#  9. WIP : CUDA                #
#################################

#### 1. Declarations ####
# set excutable name as project name.
set(TARGET_NAME ${PROJECT_NAME})

# Compile flags
set(CF)
# Link flags
set(LF)
# Compile definitions
set(CD)
# Header directory
set(HD)

#### 2. User Setting ####

## BUILD ENVIRONMENT OPTIONS ##
# option(USE_BLAS        "Using BLAS"       OFF)

## WIP ##
option(USE_OPENBLAS    "Using OpenBLAS"   OFF)
option(USE_MKL         "Using Intel MKL"  ON)

option(USE_OPENMP      "Using OpenMP"     OFF)

# Not supported, yet.
# option(USE_CUDA    "Using CUDA"      OFF)

set(MAIN_SRC 
	# Add your source files here.
	# including the source which contains main().
	# root directory is location of CMakeLists.txt 
	# EXAMPLE
  #	test/test_matlab.c	
  # test/test_wav.c
   test/test_fft.c
	)
	
list(APPEND HD
    # Add additional header folders here.
	# ${CMAKE_CURRENT_SOURCE_DIR} means root directory.
	# EXAMPLE
	# ${CMAKE_CURRENT_SOURCE_DIR}/my_new_header_folder
    ) 

# Additional link flags and defines can be attached.
list(APPEND CF
	# my compile option
	)
list(APPEND LF
	# my link option
  fftw3
	)
list(APPEND CD
	# my pre-process defines
	)	
	
#### 3. Macros ####

macro(display_elements
	DISPLAY_MESSAGE
	LIST_TO_DISPLAY)
    message("[ ${DISPLAY_MESSAGE} ]")
    foreach(_DISPLAY_ELEMENT ${LIST_TO_DISPLAY})
	message(STATUS ${_DISPLAY_ELEMENT})
    endforeach()
endmacro()

#### 4. sources and executable ####

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
    source/iip_fft.c
    )

add_executable(${TARGET_NAME} ${MAIN_SRC} ${C_SOURCE})

#### 5. Header includes ####

# Setting options by 'set(CMAKE_<blah>_<blah>)' is unstable. Often not working.
# Recommed using 'target_<blah>_<blah>'
# which is more robust and modern.


target_include_directories(${TARGET_NAME}
	PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/header
	$(HD)	
	)

#### 6. check options ####
# This section checks whether your option is available or not. 

#if(USE_BLAS)
if(NOT)
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

#### 7. defines and Libraries  ####
# This section sets pre-processing defines and flags specified by environment

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

# #include <math.h> for gcc 
if(NOT MSVC)
  # MKL need to set 'm' but need to be in order.
  if(NOT USE_MKL)
    list(APPEND LF m)
  endif(NOT USE_MKL)
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
 # SEE
 # https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
 
 # == Options ==
 # Intel(R) MKL 2018.0
 # Linux
 # None
 # GNU C/C++
 # Intel(R) 64
 # Dynamic
 # 32-bit integer 
 # OpenMP threading
 # Intel(R) (libmop5)

  list(APPEND LF
         -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
         # -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
    )

  list(APPEND CF
     -m64 -I${MKLROOT}/include
    )

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

target_link_libraries     (${TARGET_NAME} PRIVATE ${LF} )
target_compile_definitions(${TARGET_NAME} PRIVATE ${CD} )

##### 8. compile optinos ####
# This section adds basic compile options 

if(WIN32 AND MSVC)
  # No /O3 in VS2013 (higher, too)?

  #list(APPEND  CMAKE_C_FLAGS " /O2 /MP ") -- This method is unstable.
  list(APPEND CF
     # /O2 - fastest code in the majority of cases
     # /MP - compiled in mult-process : faster build time
     /MP
     # /O2 - not compatible with debug mode flag /RTC1
     # /O2 is set as default in release mode of visual Studio
     )
 else(WIN32 AND MSVC)
  # list(APPEND CMAKE_C_FLAGS "-g -O3 -Wall -Wextra ") : Wall is too much for now
  list(APPEND CF
      -g
      # -O3 : maximize performace
      # -O2 : maximize performance as long as not increasing size of executable
      # -O1 : minimize compile time and others are secondary
      # -Os : minimize size of executable
       -O3
      -Wextra
      )
  endif(WIN32 AND MSVC)

target_compile_options    (${TARGET_NAME} PRIVATE ${CF} )

message(STATUS "CMAKE_GENERATOR : "  ${CMAKE_GENERATOR})
display_elements("MAIN SOURCE" "${MAIN_SRC}")
display_elements("LIBRARY SOURCE" "${C_SOURCE}")
display_elements("Compile Flags" "${CF}")
display_elements("Linker Flags" "${LF}")
display_elements("Compile Defines" "${CD}")

#### 9. misc : CUDA(Work in Progress) ####

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

