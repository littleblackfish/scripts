cmake_minimum_required(VERSION 2.6)

project(gnuplotutils NONE)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if (NOT DEFINED DATA)
  set(DATA "share/gnuplotutils")
endif(NOT DEFINED DATA)
if (NOT DEFINED MAN)
  set(MAN "share/man")
endif(NOT DEFINED MAN)

find_package(TXT2TAGS)

configure_file(${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

foreach(PROG gp2eps gplot tricktex)
  if (TXT2TAGS_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/${PROG}.t2t)
    add_custom_command(OUTPUT ${PROG}.man 
      COMMAND ${TXT2TAGS_EXECUTABLE} -q -t man -i ${CMAKE_SOURCE_DIR}/${PROG}.t2t -o ${PROG}.man)
    list(APPEND MANPAGES ${PROG}.man)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROG}.man DESTINATION ${MAN}/man1 RENAME ${PROG}.1)
    add_custom_command(OUTPUT ${PROG}.wiki 
      COMMAND ${TXT2TAGS_EXECUTABLE} -q -t gwiki -i ${CMAKE_SOURCE_DIR}/${PROG}.t2t -o ${PROG}.wiki)
    list(APPEND WIKIPAGES ${PROG}.wiki)
  endif (TXT2TAGS_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/${PROG}.t2t)
  install(PROGRAMS ${PROG} DESTINATION bin)
endforeach(PROG)

if (TXT2TAGS_FOUND)
  add_custom_target(manpages ALL DEPENDS ${MANPAGES})
  add_custom_target(wiki DEPENDS ${WIKIPAGES})
endif (TXT2TAGS_FOUND)

install(PROGRAMS extras/demo.sh DESTINATION ${DATA})
