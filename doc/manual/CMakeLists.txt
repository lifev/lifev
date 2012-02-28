add_custom_target (manual_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../tex/lifestyle.cls;lifestyle.cls
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../tex/rcs.sty;rcs.sty
                   VERBATIM)

add_latex_document (lifev-manual.tex
                    INPUTS lifev-manual.tex lifev-install.tex lifev-examples.tex cavity.tex navierstokes.tex fsi.tex BCCond.fig  BCHandler.fig  BasElSha.fig  Fem_nodes.fig  mixed_pattern.fig  theFE.fig BCCond.eps BasElSha.eps NSCode_mesh_types.eps theFE.eps BCHandler.eps Fem_nodes.eps mixed_pattern.eps BCCond.pdf BasElSha.pdf NSCode_mesh_types.pdf theFE.pdf BCHandler.pdf Fem_nodes.pdf mixed_pattern.pdf cavityFigure.pdf_t cavityFigure.pdf ALEmapping.pdf
                    BIBFILES ../tex/lifev.bib
                    DEFAULT_PDF
                    MANGLE_TARGET_NAMES
                    USE_INDEX)

add_latex_document (lifev-dev.tex
                    INPUTS lifev-dev.tex lifev-install.tex lifev-examples.tex lifev-dev_conventions.tex lifev-dev_overview.tex lifev-dev_howto.tex BCCond.fig  BCHandler.fig  BasElSha.fig  Fem_nodes.fig  mixed_pattern.fig  theFE.fig BCCond.eps BasElSha.eps NSCode_mesh_types.eps theFE.eps BCHandler.eps Fem_nodes.eps mixed_pattern.eps BCCond.pdf BasElSha.pdf NSCode_mesh_types.pdf theFE.pdf BCHandler.pdf Fem_nodes.pdf mixed_pattern.pdf cavityFigure.pdf_t cavityFigure.pdf ALEmapping.pdf
                    BIBFILES ../tex/lifev.bib
                    DEFAULT_PDF
                    MANGLE_TARGET_NAMES
                    USE_INDEX)

add_custom_target (manual)
add_custom_target (dev-guide)
add_dependencies (manual DEPENDS lifev-manual_pdf)
add_dependencies (dev-guide DEPENDS lifev-dev_pdf)
add_dependencies (lifev-manual_pdf DEPENDS manual_links)
add_dependencies (lifev-dev_pdf DEPENDS manual_links)

