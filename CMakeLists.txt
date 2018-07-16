#CMAKE 3.9 버전 부터 CUDA관련 함수 및 기능 대폭 변경

#CMAKE 3.10.3 이용



cmake_minimum_required(VERSION 3.10.3)

set(PROJECT_NAME iip_sph_pp)



#include(CMakeDependentOption)



###################################################

###################################################

###################################################



option(USE_CUDA "Using CUDA"                 OFF)



#option(BLAS_CU  " BLAS : cublas"             ON)

option(USE_OPEN "BLAS : OpenBLAS"           OFF)

option(USE_MKL "BLAS : IntelMKL"            OFF)

#OpenMP

option(USE_OPENMP "Using OpenMP"             ON)



#설정을 다 하셨다면 이 옵션을 ON 해주세요

option(I_DID_WHAT_I_HAVE_TO_DO "SET THIS ON" ON)



###################################################

###################################################

###################################################



set(LINKLIBS)

set(MAIN_SRC)

set(C_SRC)

set(CUDA_SRC)



### MAIN FILE ###### MAIN FILE ###### MAIN FILE ###

### MAIN FILE ###### MAIN FILE ###### MAIN FILE ###


#set(MAIN_SRC test/test_matrix.c)
#set(MAIN_SRC test/test_wav.c)
#set(MAIN_SRC test/test_blas_lv1.c)
#set(MAIN_SRC test/test_blas_lv2.c)
set(MAIN_SRC test/test_blas_lv3.c)

#set(MAIN_SRC test/test_blas_lv2_big.c)



### MAIN FILE ###### MAIN FILE ###### MAIN FILE ###

### MAIN FILE ###### MAIN FILE ###### MAIN FILE ###



if(NOT I_DID_WHAT_I_HAVE_TO_DO)

	message(FATAL_ERROR "*****************************************\nYOU NEED TO SET CONFIGURATION OPTIONS MANUALLY\n****************************************")

endif()





message(STATUS "USE_CUDA   : " ${USE_CUDA})

message(STATUS "USE_OPEN  : " ${USE_OPEN})

message(STATUS "USE_MKL   : " ${USE_MKL})

message(STATUS "USE_OPENMP : " ${USE_OPENMP})



#GET OS

#UNIX   TRUE for UNIX,APPLE

#WIN32  TRUE for Windows 32bit and 64bit

#APPLE  TRUE for APPLE

if(UNIX)

	add_definitions(-DOS_UNIX=1)

elseif(WIN32)

	add_definitions(-DOS_WIN=1)

endif()



include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test)

message(STATUS "ENV_PATH : " $ENV{PATH})



if(USE_CUDA)

	project(${PROJECT_NAME} CUDA)
	list(APPEND CUDA_SRC iip_matrix.cu iip_wav.cu iip_blas_lv2.cu iip_blas_lv3.cu)

	

	add_definitions(-DUSE_CUDA=1)

		

	# compile .c files by nvcc

	#SET_SOURCE_FILES_PROPERTIES(${COMMON_SRC} PROPERTIES LANGUAGE CUDA)

	if(WIN32)

	enable_language(CXX)

	SET_SOURCE_FILES_PROPERTIES(${MAIN_SRC} PROPERTIES LANGUAGE CXX)	

	else(WIN32)
	SET_SOURCE_FILES_PROPERTIES(${MAIN_SRC} PROPERTIES LANGUAGE CUDA)
	endif(WIN32)

	list(APPEND LINKLIBS ${CUDA_CUBLAS_LIBRARIES})  



	add_executable(out_main ${MAIN_SRC} ${CUDA_SRC})		



	target_link_libraries(out_main ${LINKLIBS} cublas)

	#	add_executable(out ${COMMON_SRC} ${CUDA_SRC})	

	#	list(APPEND LINKLIBS ${CUDA_CUBLAS_LIBRARIES})

	#	target_link_libraries(out ${LINKLIBS} cublas)

  



else()

	project(${PROJECT_NAME} C)

	

	#eliminate '/W3' : D8021  invalid numeric argument error

	string (REPLACE "/W3" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}") 

	

	message(STATUS "CMAKE_C_FLAGS : " ${CMAKE_C_FLAGS})	

	message(STATUS "CMAKE_EXE_LINKER_FLAGS : " ${CMAKE_EXE_LINKER_FLAGS})

	

	list(APPEND C_SRC iip_matrix.c iip_blas_lv1.c iip_wav.c iip_blas_lv2.c iip_blas_lv3.c)

		

	if(USE_OPEN)

		find_package(BLAS REQUIRED)

		if(NOT BLAS_FOUND)

			message(FATAL_ERROR "Can't find BLAS")

		endif()

		add_definitions(-DUSE_OPEN=1)

		add_definitions(-DUSE_CBLAS=1)

	

		message(STATUS ${BLAS_LIBRARIES})

		list(APPEND LINKLIBS ${BLAS_LIBRARIES})

	

	elseif(USE_MKL)

#linux 라면 mkl 사용시 컴파일 옵션을 구질구질하게 걸어줘야하는데

#Windows VS를 쓰면 mkl 플러그인이 있어서 쓰기 편하다

#그렇다면 Window일땐 그냥 냅두는것이?

	if(UNIX)

		#MKL 은 원래도 수동으로 경로를 잡아서 cmake가 해줄순 없다

		set(MKL_PATH)

		if("${MKL_PATH}" STREQUAL "")

				message(FATAL_ERROR "ADD AN ARGUMENT '-DMKL_PATH=<ALL REQUIRED MKL libs>'")

		endif()

		message(STATUS "MKL_PATH : " ${MKL_PATH})

		list(APPEND LINKLIBS ${MKL_PATH})

	endif()	

	

	add_definitions(-DUSE_MKL=1)

	add_definitions(-DUSE_CBLAS=1)



	endif()

	



	if(USE_OPENMP)

		find_package(OpenMP REQUIRED)

		if(NOT OPENMP_FOUND)

			message(FATAL_ERROR "Can't find OpenMP")

		endif()

		list(APPEND CMAKE_C_FLAGS ${OpenMP_C_FLAGS})

		list(APPEND CMAKE_EXE_LINKER_FLAGS ${OpenMP_EXE_LINKER_FLAGS})

	endif()



	add_executable(out_main ${MAIN_SRC} ${C_SRC})	



	target_link_libraries(out_main ${LINKLIBS})





endif()



#generic function 
