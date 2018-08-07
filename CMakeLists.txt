cmake_minimum_required(VERSION 3.5)
project(iip_sph_pp)

###################################################
###################################################
###################################################

option(USE_CUDA    "Using CUDA"      OFF)
option(USE_BLAS    "Using BLAS"      ON)

option(USE_OPENMP  "Using OpenMP"    ON)

#설정을 다 하셨다면 이 옵션을 ON 해주세요
option(I_DID_WHAT_I_HAVE_TO_DO "SET THIS ON" ON)
if(NOT I_DID_WHAT_I_HAVE_TO_DO)
    message(FATAL_ERROR "*****************************************\nYOU NEED TO SET CONFIGURATION OPTIONS MANUALLY\n****************************************")
endif()
### MAIN FILE ###### MAIN FILE ###### MAIN FILE ###
### MAIN FILE ###### MAIN FILE ###### MAIN FILE ###
#set(MAIN_SRC test/test_matrix.c)
#set(MAIN_SRC test/test_matrix_2.c)
#set(MAIN_SRC test/test_wav.c)
#set(MAIN_SRC test/test_blas_lv1.c)
#set(MAIN_SRC test/test_blas_lv2.c)
#set(MAIN_SRC test/test_blas_lv3.c)
#set(MAIN_SRC test/test_memory.c)
#set(MAIN_SRC test/test_aABpbC.c)
#set(MAIN_SRC test/test_broad.c)
#set(MAIN_SRC test/test_broad_matmul.c)
#set(MAIN_SRC test/test_dia_tr.c)
#set(MAIN_SRC test/test_invert.c)
#set(MAIN_SRC test/test_cinvert.c)
#set(MAIN_SRC test/test_dim.c)
#set(MAIN_SRC test/test_matlab.c)
set(MAIN_SRC test/test_combine.c)
#set(MAIN_SRC test/test_permute.c)
#set(MAIN_SRC test/test_bin.c)
#set(MAIN_SRC test/test_invelem.c)
#set(MAIN_SRC test/test_invers.c)
#set(MAIN_SRC test/test_row_col.c)
#set(MAIN_SRC test/test_rand.c)
#set(MAIN_SRC test/test_math.c)
#set(MAIN_SRC test/test_trans.c)
#set(MAIN_SRC test/openBLAS_with_openMP.c)
#set(MAIN_SRC test/test_blas_lv2_big.c)

### MAIN FILE ###### MAIN FILE ###### MAIN FILE ###
### MAIN FILE ###### MAIN FILE ###### MAIN FILE ###


set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

if(WIN32 AND MSVC)
    set(BUILD_FLAGS
	"${CMAKE_C_FLAGS} /O3 /MP")
  #There is something wrong about the way CMAKE set /W3 flag
    string (REPLACE "/W3" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}") 
else()
    if(WIN32)
	set(BUILD_FLAGS
	    "${CMAKE_C_FLAGS} -g -O3 -Wall -Wextra")
    else()
	set(BUILD_FLAGS
            "${CMAKE_C_FLAGS} -g -O3 -Wall -Wextra")
    endif()
endif()


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
    find_package(OpenMP)
    if(NOT(OpenMP_FOUND))
		message(SEND_ERROR "USE_OPENMP is set ON but OpenMP not found")
		set(USE_OPENMP OFF)

    else()
	message(STATUS "OpenMP found")
	set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
	if(MSVC)
	    # Visual Studio uses '/' for compile option
	    # CMAKE set -openmp option for OpenMP_FLAGS
	    string(REPLACE "-openmp" " /openmp" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
	endif(MSVC)
    endif()
endif()


message(STATUS "USE_CUDA   : " ${USE_CUDA})
message(STATUS "USE_OPEN  : " ${USE_OPEN})
message(STATUS "USE_MKL   : " ${USE_MKL})
message(STATUS "USE_OPENMP : " ${USE_OPENMP})

#GET OS

#UNIX   TRUE for UNIX,APPLE

#WIN32  TRUE for Windows 32bit and 64bit

#APPLE  TRUE for APPLE

macro(display_elements
	DISPLAY_MESSAGE
	LIST_TO_DISPLAY)
    message("[ ${DISPLAY_MESSAGE} ]")
    foreach(_DISPLAY_ELEMENT ${LIST_TO_DISPLAY})
	message(STATUS ${_DISPLAY_ELEMENT})
    endforeach()
endmacro()

if(UNIX)
    add_definitions(-DOS_UNIX=1)
elseif(WIN32)
    add_definitions(-DOS_WIN=1)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/header)
include_directories()
message(STATUS "ENV_PATH : " $ENV{PATH})

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
display_elements("Source files for C build" "${C_SOURCE}")

if(USE_OPEN OR USE_MKL )
    set(BLAS_TARGET_NAME "${PROJECT_NAME}_blas")

    message(STATUS "CMAKE_C_FLAGS : " ${CMAKE_C_FLAGS})	
    message(STATUS "CMAKE_EXE_LINKER_FLAGS : " ${CMAKE_EXE_LINKER_FLAGS})

    string (REPLACE "/W3" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}") 
	add_executable(${BLAS_TARGET_NAME} ${MAIN_SRC} ${C_SOURCE})
    string (REPLACE ";" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}") 

    if(USE_OPEN)
	message(STATUS ${BLAS_LIBRARIES})
  if(NOT MSVC)  
    target_link_libraries(${BLAS_TARGET_NAME}
  	m
    )
  endif(NOT MSVC)
	target_link_libraries(${BLAS_TARGET_NAME}
	    openblas)
	target_include_directories(${BLAS_TARGET_NAME}
	    PUBLIC
	    ${CMAKE_CURRENT_SOURCE_DIR}/header)
	target_compile_definitions(${BLAS_TARGET_NAME}
	    PRIVATE
	    USE_OPEN=1
	    USE_CBLAS=1)
    elseif(USE_MKL)
  
  # 리눅스 gnu-compiler intel-64 <linking> 64bit <threading> limiomp5 
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
  list(APPEND CMAKE_C_FLAGS -DMKL_ILP64 -m64 -I${MKLROOT}/include )
  list(APPEND CMAKE_EXE_LINKER_FLAGS -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
 )
	endif()	
	
	target_compile_definitions(${BLAS_TARGET_NAME}
	    PRIVATE
	    USE_MKL=1)
		
    endif()
elseif(USE_OPENMP)
    set(OPENMP_TARGET_NAME "${PROJECT_NAME}_openmp")
    #    list(APPEND CMAKE_C_FLAGS )

    #This Line Causes Error, ';' no such file or directory
    #		list(APPEND CMAKE_EXE_LINKER_FLAGS ${OpenMP_EXE_LINKER_FLAGS})
    
    add_definitions(-DUSE_OMP=1)
    string (REPLACE "/W3" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}") 
    string (REPLACE ";" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}") 
  
    list(APPEND CMAKE_EXE_LINKER_FLAGS )

	add_executable(${OPENMP_TARGET_NAME} ${MAIN_SRC} ${C_SOURCE})	
    target_link_libraries(${OPENMP_TARGET_NAME} )

  #  set_target_properties(${OPENMP_TARGET_NAME}
  #PROPERTIES COMPILE_FLAGS "${BUILD_FLAGS}"
  #     CXX_STANDARD 14
  #      CXX_STANDARD_REQUIRED ON)
  if(NOT MSVC)  
  target_link_libraries(${OPENMP_TARGET_NAME}
	m)
  endif(NOT MSVC)
  target_include_directories(${OPENMP_TARGET_NAME}
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/header)



  #Fully native C
else()
	set(NATIVE_TARGET_NAME "${PROJECT_NAME}_native")
	add_executable(${NATIVE_TARGET_NAME} ${MAIN_SRC} ${C_SOURCE})	
  If(UNIX)
  target_link_libraries(${NATIVE_TARGET_NAME} m)
  endif(UNIX)
	target_include_directories(${NATIVE_TARGET_NAME}
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/header)



endif()

message(STATUS "CMAKE_GENERATOR : "  ${CMAKE_GENERATOR})
