## ======================================================================== ##
## Copyright 2018-2019 Ingo Wald                                            ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##     http://www.apache.org/licenses/LICENSE-2.0                           ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
## ======================================================================== ##

find_package(OptiX)

find_program(BIN2C bin2c
  DOC "Path to the cuda-sdk bin2c executable.")

macro(cuda_compile_and_embed output_var cuda_file)
  set(var_name ${output_var})
  cuda_compile_ptx(ptx_files ${cuda_file})
  list(GET ptx_files 0 ptx_file)
  set(embedded_file ${ptx_file}_embedded.c)
  add_custom_command(
    OUTPUT ${embedded_file}
    COMMAND ${BIN2C} -c --padd 0 --type char --name ${var_name} ${ptx_file} > ${embedded_file}
    DEPENDS ${ptx_file}
    COMMENT "compiling (and embedding ptx from) ${cuda_file}"
    )
  set(${output_var} ${embedded_file})
endmacro()

#set(CUDA_HOST_COMPILER gcc)
set(CUDA_GENERATED_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cuda_output_dir)

include_directories(${OptiX_INCLUDE})
