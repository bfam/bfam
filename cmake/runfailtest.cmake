string(REPLACE " " ";" TEST_PROG_LIST ${TEST_PROG})
execute_process(COMMAND ${TEST_PROG_LIST}
                RESULT_VARIABLE HAD_ERROR
                OUTPUT_VARIABLE OUTPUT)
message("${OUTPUT}")

message("Expected Error: ${HAD_ERROR}")
if(HAD_ERROR)
  message("Test passed")
else()
  message(FATAL_ERROR "Test failed")
endif()
