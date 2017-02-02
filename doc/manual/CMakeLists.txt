
SET(TEX_FILES
  lifev-manual
  lifev-dev
)

SET(FILES_lifev-manual
  lifev-manual.tex lifev-install.tex lifev-examples.tex laplacian.tex cavity.tex navierstokes.tex fsi.tex BCCond.fig  BCHandler.fig  BasElSha.fig  Fem_nodes.fig  mixed_pattern.fig  theFE.fig BCCond.eps BasElSha.eps NSCode_mesh_types.eps theFE.eps BCHandler.eps Fem_nodes.eps mixed_pattern.eps BCCond.pdf BasElSha.pdf NSCode_mesh_types.pdf theFE.pdf BCHandler.pdf Fem_nodes.pdf mixed_pattern.pdf cavityFigure.pdf_t cavityFigure.pdf ALEmapping.pdf
)

SET(FILES_lifev-dev
 lifev-dev.tex lifev-install.tex lifev-examples.tex lifev-dev_conventions.tex lifev-dev_overview.tex lifev-dev_howto.tex BCCond.fig  BCHandler.fig  BasElSha.fig  Fem_nodes.fig  mixed_pattern.fig  theFE.fig BCCond.eps BasElSha.eps NSCode_mesh_types.eps theFE.eps BCHandler.eps Fem_nodes.eps mixed_pattern.eps BCCond.pdf BasElSha.pdf NSCode_mesh_types.pdf theFE.pdf BCHandler.pdf Fem_nodes.pdf mixed_pattern.pdf cavityFigure.pdf_t cavityFigure.pdf ALEmapping.pdf
)
SET(BIB_FILES 
  ../tex/lifev.bib
)



CONFIGURE_FILE(ltxcompile.cmake
               ltxcompile @ONLY)

FOREACH(TEX_FILE ${TEX_FILES})
  ADD_CUSTOM_TARGET(pdf_${TEX_FILE} echo
    DEPENDS ${TEX_FILE}.pdf
  )
  ADD_DEPENDENCIES(pdf pdf_${TEX_FILE})

  ADD_CUSTOM_COMMAND(
    OUTPUT    ${TEX_FILE}.pdf
    COMMAND   /bin/sh ltxcompile
    ARGS      ${TEX_FILE}.tex
    DEPENDS   ${FILES_${TEX_FILE}} ${BIB_FILES}
  )

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${TEX_FILE}.pdf
    DESTINATION share/doc/lifev
    COMPONENT doc
    OPTIONAL
  )
ENDFOREACH()
