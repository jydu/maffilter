# CMake script for MafFilter
# Authors:
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 22/08/2009

# Builds info, html, pdf doc.
# Info doc is built and install as part of "all" if makeinfo is found.
# Html doc is proposed as a "html" optional target if makeinfo is found.
# Pdf doc is proposed as a "pdf" optional target if makeinfo AND texi2dvi are found.

find_program (MAKEINFO NAMES makeinfo texi2any DOC "makeinfo doc generator program")
if (NOT MAKEINFO)
  message (STATUS "makeinfo program not found: 'info' and 'html' target disabled (builds info/html doc)")
else (NOT MAKEINFO)
  message (STATUS "Found makeinfo as '${MAKEINFO}': 'info' and 'html' target enabled (builds info/html doc)")

  set (input ${CMAKE_CURRENT_SOURCE_DIR}/maffilter.texi)

  # Build and install info page
  set (output ${CMAKE_CURRENT_BINARY_DIR}/maffilter.info)
  add_custom_command (
    OUTPUT ${output}
    COMMAND ${MAKEINFO} --no-split -o ${output} ${input}
    DEPENDS ${input}
    COMMENT "Generating info page"
    VERBATIM
    )
  install (FILES ${output} DESTINATION share/info)

  # Add "info" target, built with "all" (needed because install will fail if not built).
  add_custom_target (info ALL DEPENDS ${output})

  # Also provide a "html" target that builds html doc (not installed, and not part of "all").
  set (output ${CMAKE_CURRENT_BINARY_DIR}/maffilter.html)
  set (makeinfo-css "http://www.w3.org/StyleSheets/Core/Steely")
  add_custom_command (
    OUTPUT ${output}
    COMMAND ${MAKEINFO} --html --css-ref=${makeinfo-css} --no-split -o ${output} ${input}
    DEPENDS ${input}
    COMMENT "Generating html doc"
    VERBATIM
    )
  add_custom_target (html DEPENDS ${output})

  # Provide a "pdf" target that builds pdf doc (not installed, not part of "all").
  find_program (TEXIDVI NAMES texi2dvi)
  if (TEXIDVI)
    message (STATUS "Found texi2dvi: 'pdf' target enabled (builds pdf doc)")
    set (output ${CMAKE_CURRENT_BINARY_DIR}/maffilter.pdf)
    add_custom_command (
      OUTPUT ${output}
      COMMAND ${MAKEINFO} --pdf --Xopt=--clean -o ${output} ${input}
      DEPENDS ${input}
      COMMENT "Generating pdf doc"
      VERBATIM
      )
    add_custom_target (pdf DEPENDS ${output})
  endif (TEXIDVI)
endif (NOT MAKEINFO)

