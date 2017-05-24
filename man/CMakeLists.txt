# CMake script for the CoMap package.
# Authors:
#   Julien Dutheil
#   Francois gindraud (2017)
# Created: 22/08/2009

# Build manpages.
# In practice, they are just compressed from the text files using gzip.
# Manpages are built and installed as part of "all" if gzip is found.

find_program (GZIP NAMES gzip DOC "gzip compression program")
if (NOT GZIP)
  message (STATUS "gzip program not found: 'man' target disabled (builds manpages)")
else (NOT GZIP)
  message (STATUS "Found gzip as '${GZIP}': 'man' target enabled (builds manpages)")

  # Create a list of manpage targets
  set (manpage-targets)

  # Take all manpages files in the directory
  file (GLOB manpage_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.1)
  foreach (manpage_file ${manpage_files})
    # Compress manpage, install, add to manpage list
    set (input ${CMAKE_CURRENT_SOURCE_DIR}/${manpage_file})
    set (output ${CMAKE_CURRENT_BINARY_DIR}/${manpage_file}.gz)
    add_custom_command (
      OUTPUT ${output}
      COMMAND ${GZIP} -c ${input} > ${output}
      DEPENDS ${input}
      COMMENT "Compressing manpage ${manpage_file}"
      VERBATIM
      )
    install (FILES ${output} DESTINATION share/man/man1)
    list (APPEND manpage-targets ${output})
    unset (input)
    unset (output)
  endforeach (manpage_file)

  # Add target "man", built with "all" (needed because install will fail if not built).
  add_custom_target (man ALL DEPENDS ${manpage-targets})
endif (NOT GZIP)
