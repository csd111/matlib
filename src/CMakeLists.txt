
add_library(matlib INTERFACE)

# - headers
install(FILES
  ${src}/vector.h
  ${src}/matrix.h
  ${src}/matrix_op.h
  DESTINATION include/matlib
)

# - external libraries
set(external_libraries "")
set(external_includes "")

# - if we want to build the tests
if (${matlib_enable_tests})

  add_executable(matlib_test
    ${src}/test.cc
  )
  add_dependencies(matlib_test
    external_googletest
  )
  if (UNIX AND NOT APPLE)
    target_link_libraries(matlib_test
      gtest
      gtest_main
      matlib
      ${external_libraries}
    )
  else ()
    target_link_libraries(matlib_test
      optimized gtest
      optimized gtest_main
      debug gtestd
      debug gtest_maind
      matlib
      ${external_libraries}
    )
  endif()
  target_include_directories(matlib_test
    PUBLIC
      ${gtest_install_dir}/include
      ${external_includes}
  )

endif()