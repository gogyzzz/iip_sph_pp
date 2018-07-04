#CMAKE 3.9 버전 부터 CUDA관련 함수 및 기능 대폭 변경
#CMAKE 3.10.3 이용

cmake_minimum_required(VERSION 3.10.3)
set(PROJECT_NAME iip_sph_pp)

#include(CMakeDependentOption)

###################################################
###################################################
###################################################

#BLAS
option(USE_CUDA "Using CUDA"                 OFF)

option(BLAS_CU  " BLAS : cublas"             OFF)
option(BLAS_OPEN "BLAS : OpenBLAS"           OFF)
option(BLAS_MKL "BLAS : IntelMKL"            OFF)

#OpenMP
option(USE_OPENMP "Using OpenMP"             OFF)


#설정을 다 하셨다면 이 옵션을 ON 해주세요
option(I_DID_WHAT_I_HAVE_TO_DO "SET THIS ON" ON)

###################################################
###################################################
###################################################
message(STATUS "USE_CUDA   : " ${USE_CUDA})
message(STATUS "BLAS_OPEN  : " ${BLAS_OPEN})
message(STATUS "BLAS_MKL   : " ${BLAS_MKL})
message(STATUS "USE_OPENMP : " ${USE_OPENMP})

set(LINKLIBS)
set(COMMON_SRC)
set(C_SRC)
set(CUDA_SRC)


list(APPEND COMMON_SRC main.c)


if(NOT I_DID_WHAT_I_HAVE_TO_DO)
			message(FATAL_ERROR "*****************************************\nYOU NEED TO SET CONFIGURATION OPTIONS MANUALLY\n****************************************")
endif()


#Untested
if(USE_CUDA)
	project(${PROJECT_NAME} CUDA)
	list(APPEND CUDA_SRC iip_matrix.cu)
	
	add_definitions(-DUSE_CUDA=1)

	if(BLAS_CU)
	add_definitions(-DBLAS_CU=1)
	endif()

#.c 파일도 다 .cu로 컴파일
	SET_SOURCE_FILES_PROPERTIES(${COMMON_SRC} PROPERTIES LANGUAGE CUDA)
	add_executable(out ${COMMOC_SRC} ${CUDA_SRC})
#	CUDA_ADD_EXECUTABLE(out ${COMMON_SRC} ${CUDA_SRC})
#	CUDA_ADD_CUBLAS_TO_TARGET(out)
	#CUDA_ADD_CUFFT_TO_TARGET(out)

else()
	project(${PROJECT_NAME} C)
	list(APPEND C_SRC iip_matrix.c iip_blas_lv1.c)

	if(BLAS_OPEN)
		find_package(BLAS REQUIRED)
		if(NOT BLAS_FOUND)
			message(FATAL_ERROR "Can't find BLAS")
		endif()
		add_definitions(-DUSE_CBLAS=1)
		message(STATUS ${BLAS_LIBRARIES})
		list(APPEND LINKLIBS ${BLAS_LIBRARIES})
	
	elseif(BLAS_MKL)
		#MKL 은 원래도 수동으로 경로를 잡아서 cmake가 해줄순 없다
		set(MKL_PATH)
		if("${MKL_PATH}" STREQUAL "")
				message(FATAL_ERROR "ADD AN ARGUMENT '-DMKL_PATH=<ALL REQUIRED MKL libs>'")
		endif()
		add_definitions(-DUSE_CBLAS=1)
		message(STATUS "MKL_PATH : " ${MKL_PATH})
		list(APPEND LINKLIBS ${MKL_PATH})
	endif()


	if(USE_OPENMP)
		find_package(OpenMP REQUIRED)
		if(NOT OPENMP_FOUND)
			message(FATAL_ERROR "Can't find OpenMP")
		endif()
		list(APPEND CMAKE_C_FLAGS ${OpenMP_C_FLAGS})
		list(APPEND CMAKE_EXE_LINKER_FLAGS ${OpenMP_EXE_LINKER_FLAGS})
	endif()



	add_executable(out ${COMMON_SRC} ${C_SRC})
endif()

#generic function 
target_link_libraries(out ${LINKLIBS})
