
#CMAKE 3.9 버전 부터 CUDA관련 함수 및 기능 대폭 변경

#CMAKE 3.10.3 이용


cmake_minimum_required(VERSION 3.5)
project(iip_sph_pp)

###################################################
###################################################
###################################################

option(USE_CUDA    "Using CUDA"      OFF)
option(USE_OPEN    "Using BLAS"      OFF)
option(USE_OPENMP  "Using OpenMP"    ON)

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
#set(MAIN_SRC test/test_dim.c)
#set(MAIN_SRC test/test_permute.c)
set(MAIN_SRC test/test_bin.c)
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




message(STATUS "USE_CUDA   : " ${USE_CUDA})
message(STATUS "USE_OPEN  : " ${USE_OPEN})
message(STATUS "USE_MKL   : " ${USE_MKL})
message(STATUS "USE_OPENMP : " ${USE_OPENMP})

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

if(WIN32 AND MSVC)
    set(BUILD_FLAGS
	"${CMAKE_C_FLAGS} /O3 /MP")
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
			message(STATUS "found OpenBLAS")
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

#설정을 다 하셨다면 이 옵션을 ON 해주세요

option(I_DID_WHAT_I_HAVE_TO_DO "SET THIS ON" ON)

###################################################
###################################################
###################################################


if(NOT I_DID_WHAT_I_HAVE_TO_DO)
    message(FATAL_ERROR "*****************************************\nYOU NEED TO SET CONFIGURATION OPTIONS MANUALLY\n****************************************")
endif()
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

    )
display_elements("Source files for C build" "${C_SOURCE}")

if(USE_OPEN OR USE_MKL )
    set(BLAS_TARGET_NAME "${PROJECT_NAME}_blas")

    message(STATUS "CMAKE_C_FLAGS : " ${CMAKE_C_FLAGS})	
    message(STATUS "CMAKE_EXE_LINKER_FLAGS : " ${CMAKE_EXE_LINKER_FLAGS})

	add_executable(${BLAS_TARGET_NAME} ${MAIN_SRC} ${C_SOURCE})
    set_target_properties(${BLAS_TARGET_NAME}
	PROPERTIES COMPILE_FLAGS "${BUILD_FLAGS}"
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON)

    if(USE_OPEN)
	message(STATUS ${BLAS_LIBRARIES})
	target_link_libraries(${BLAS_TARGET_NAME}
	    openblas m)
	target_include_directories(${BLAS_TARGET_NAME}
	    PUBLIC
	    ${CMAKE_CURRENT_SOURCE_DIR}/header)
	target_compile_definitions(${BLAS_TARGET_NAME}
	    PRIVATE
	    USE_OPEN=1
	    USE_CBLAS=1)
    elseif(USE_MKL)
	
	#linux 라면 mkl 사용시 컴파일 옵션을 구질구질하게 걸어줘야하는데
	#Windows VS를 쓰면 mkl 플러그인이 있어서 쓰기 편하다
	#그렇다면 Window일땐 그냥 냅두는것이?

	#	if(UNIX)
	#   #MKL 은 원래도 수동으로 경로를 잡아서 cmake가 해줄순 없다
	#    set(MKL_PATH)
	#    if("${MKL_PATH}" STREQUAL "")
	#		message(FATAL_ERROR "ADD AN ARGUMENT '-DMKL_PATH=<ALL REQUIRED MKL libs>'")
	##    endif()
	#    message(STATUS "MKL_PATH : " ${MKL_PATH})
	#    target_link_libraries(${BLAS_TARGET_NAME}
	#		${MKL_PATH})
	#endif()	
	
	target_compile_definitions(${BLAS_TARGET_NAME}
	    PRIVATE
	    USE_MKL=1)
    endif()
elseif(USE_OPENMP)
    set(OPENMP_TARGET_NAME "${PROJECT_NAME}_openmp")
    list(APPEND CMAKE_C_FLAGS )

    #This Line Causes Error, ';' no such file or directory
    #		list(APPEND CMAKE_EXE_LINKER_FLAGS ${OpenMP_EXE_LINKER_FLAGS})
    
    add_definitions(-DUSE_OMP=1)
	add_executable(${OPENMP_TARGET_NAME} ${MAIN_SRC} ${C_SOURCE})	

    set_target_properties(${OPENMP_TARGET_NAME}
	PROPERTIES COMPILE_FLAGS "${BUILD_FLAGS}"
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON)
    target_link_libraries(${OPENMP_TARGET_NAME}
	m)
    target_include_directories(${OPENMP_TARGET_NAME}
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/header)
else()
	set(NATIVE_TARGET_NAME "${PROJECT_NAME}_native")
	add_executable(${NATIVE_TARGET_NAME} ${MAIN_SRC} ${C_SOURCE})	
	target_link_libraries(${NATIVE_TARGET_NAME} m)
	target_include_directories(${NATIVE_TARGET_NAME}
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/header)
endif()

message(STATUS "CMAKE_GENERATOR : "  ${CMAKE_GENERATOR})
