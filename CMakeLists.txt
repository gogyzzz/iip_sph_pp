cmake_minimum_required(VERSION 3.5)
project(iip_sph_pp)

set(TARGET_NAME ${PROJECT_NAME})

#### BUILD ENVIRONMENT OPTIONS ######
option(USE_CUDA    "Using CUDA"      OFF) # Not supported, yet.
option(USE_BLAS    "Using BLAS"      OFF)
option(USE_OPENMP  "Using OpenMP"    ON)

#### set Your Source files ####

# example for adding code
# 
# set(MAIN_SRC
# my_main.c
# src/my_c_code.c
# src2/my_other_c_code.c
# )
#
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/my_header_folder)

set(MAIN_SRC test/test_matlab.c)
#set(MAIN_SRC test/chunk.c)


###########################################
#       CMAKE INDEX                       #
# 1. functions                            #
# 2. includes                             #
# 3. sources and executable               #
# 4. check options           + Order      #
# 5. defines and links        a. BLAS     #
# 6. compile optinos          b. OMP      #
# 7. misc : CUDA(WIP)         c. Native C #
###########################################

# set options by 'set' is unstable. Often not working.
# recommed using 'target_<blah>_<blah>'
# which is robust and modern.

#### 1. functions ####

macro(display_elements
	DISPLAY_MESSAGE
	LIST_TO_DISPLAY)
    message("[ ${DISPLAY_MESSAGE} ]")
    foreach(_DISPLAY_ELEMENT ${LIST_TO_DISPLAY})
	message(STATUS ${_DISPLAY_ELEMENT})
    endforeach()
endmacro()

#### 2. includes ####

# Below method is better. But, in this point, 'TARGET_NAME' isn't set.
# target_include_directories(${TARGET_NAME}
#	PUBLIC
#	${CMAKE_CURRENT_SOURCE_DIR}/header)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/header)


# MKL also can be found by find_package(BLAS)
# But, MKL requires specific compile & link option.
# just using "target_link_libraries" doesn't work for MKL.
# In case you want to use MKL in linux. I wrote MKL option list BLAS section below.
# If you modify this cmake with those options, you can use MKL in linux.

#### 3. sources and executable ####

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
  display_elements("MAIN SOURCE" "${MAIN_SRC}")
  display_elements("LIBRARY SOURCE" "${C_SOURCE}")

add_executable(${TARGET_NAME} ${MAIN_SRC} ${C_SOURCE})

#### 4. check options ####

if(USE_BLAS)
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
# Not going to use omp.h so, just linker option is required.   
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

#GET OS
#UNIX   TRUE for UNIX,APPLE
#WIN32  TRUE for Windows 32bit and 64bit
#APPLE  TRUE for APPLE
if(UNIX)
    target_compile_definitions(${TARGET_NAME}
    PRIVATE
    OS_UNIX=1
    )
elseif(WIN32)
    target_compile_definitions(${TARGET_NAME}
    PRIVATE
    OS_WIN=1
    )
endif()


message(STATUS "ENV_PATH : " $ENV{PATH})

if(NOT MSVC) 
  #include <math.h> 
  target_link_libraries(${TARGET_NAME}
  PRIVATE
	m
  )
endif(NOT MSVC)
	
if(USE_OPEN OR USE_MKL )

 if(USE_OPEN)
	message(STATUS ${BLAS_LIBRARIES})
  target_link_libraries(${TARGET_NAME}
	    openblas)
    target_compile_definitions(${TARGET_NAME}
	    PRIVATE
	    USE_OPEN=1
	    USE_CBLAS=1)
    elseif(USE_MKL)
  # MKL Options 
  # Linux gnu-compiler intel-64 <linking> 64bit <threading> limiomp5 
  # 
  #   1. dynamic 
  #       2. OMP
  #         + link :   -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
  #         + compile :  -DMKL_ILP64 -m64 -I${MKLROOT}/include 
  #       
  #       2. sequential
  #         + link :  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
  #         + compile :  -DMKL_ILP64 -m64 -I${MKLROOT}/include
  #
  #   1. static
  #       2. OMP
  #         + line :  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
  #         + compile :  -DMKL_ILP64 -m64 -I${MKLROOT}/include
  #
  #       2. sequential
  #         + link :  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
  #         + compile :  -DMKL_ILP64 -m64 -I${MKLROOT}/include
  #
  #   
  
	if(UNIX)
	#    set(MKL_PATH)
	#    if("${MKL_PATH}" STREQUAL "")
	#		message(FATAL_ERROR "ADD AN ARGUMENT '-DMKL_PATH=<ALL REQUIRED MKL libs>'")
	##    endif()
	#    message(STATUS "MKL_PATH : " ${MKL_PATH})
	#    target_link_libraries(${BLAS_TARGET_NAME}
	#		${MKL_PATH})
	endif()	
	
  target_compile_definitions(${TARGET_NAME}
	    PRIVATE
	    USE_MKL=1)
		
    endif()
endif()

if(USE_OPENMP)
   #    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
   #    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS})
	if(MSVC)
	    # Visual Studio uses '/' for compile option
     target_compile_options(${TARGET_NAME}
     PRIVATE
     /openmp
   )
   # MSVC works well with just compile option /openmp 
  
else(MSVC)
    target_compile_options(${TARGET_NAME}
     PRIVATE
     -fopenmp
     )
   target_link_libraries(${TARGET_NAME}
     PRIVATE
     -fopenmp  
     )
  endif(MSVC)
    target_compile_definitions(
      ${TARGET_NAME}
      PRIVATE
      USE_OMP=1
      )
    # add_definitions(-DUSE_OMP=1)
endif(USE_OPENMP)

##### 6. compile optinos ####

if(WIN32)
  #There is something wrong about the way CMAKE set /W3 flag
  string (REPLACE "/W3" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}") 
  string (REPLACE ";" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}") 
endif(WIN32)

if(WIN32 AND MSVC)
  # No /O3 in VS2013 (higher, too)?

  # /O2 - fastest code in the majority of cases
  # /MP - compiled in mult-process : faster build time
  #list(APPEND  CMAKE_C_FLAGS " /O2 /MP ") -- This method is unstable.
     target_compile_options(${TARGET_NAME}
     PRIVATE
     /MP
     # /O2 - not compatible with debug mod flag /RTC1
     # you had to set manuallu /O2 in release mod
     )
else()
  # list(APPEND CMAKE_C_FLAGS "-g -O3 -Wall -Wextra ") : Wall is too much for now
     target_compile_options(${TARGET_NAME}
     PRIVATE
      -g
      -O3
      -Wextra
      )
endif()



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


message(STATUS "CMAKE_GENERATOR : "  ${CMAKE_GENERATOR})


