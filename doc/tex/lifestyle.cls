\def\filedate{2003/10/29}
\def\fileversion{1}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{lifestyle}[2003/10/29]
\message{`lifestyle v\fileversion, \filedate\space (Christophe Prud'homme <christophe.prudhomme@epfl.ch>)}

\RequirePackage{ifthen}

\newboolean{life@article} \setboolean{life@article}{true}
\newboolean{life@report}  \setboolean{life@report}{false}
\newboolean{life@toc}     \setboolean{life@toc}{false}
\newboolean{life@stretch} \setboolean{life@stretch}{false}

\DeclareOption{article}  {\setboolean{life@article}{true}\setboolean{life@report}{false}}
\DeclareOption{report}   {\setboolean{life@report}{true}\setboolean{life@article}{false}}
\DeclareOption{toc}      {\setboolean{life@toc}{true}}
\DeclareOption{stretch}  {\setboolean{life@stretch}{true}}

\ifthenelse{\boolean{life@article}}
{
  \DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
}{
  \ifthenelse{\boolean{life@report}}
  {
    \DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
  }{ 
    %% default is article
    \DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
  }
}

%%
%% Execution of Options
%%

\ProcessOptions

\ifthenelse{\boolean{life@article}}
{
  \LoadClass[titlepage]{article}
}{
  \ifthenelse{\boolean{life@report}}
  {
    \LoadClass[titlepage]{report}
  }{ 
    %% default is article
    \LoadClass[titlepage]{article}
  }
}



\RequirePackage{a4wide}

\RequirePackage{fancyhdr}
\RequirePackage{xspace}
\RequirePackage{multicol}
\RequirePackage{makeidx}
\RequirePackage{subfigure}
\RequirePackage{verbatim}
\RequirePackage{rcs}

\RequirePackage[english]{babel}
\RequirePackage{tabularx}
 
\RequirePackage[pdftex,colorlinks,pdfpagemode=None,pdfview=Fit,pdfcreator={LIFE}]{hyperref}

% Graphics rules
%
\RequirePackage[pdftex]{color}
\RequirePackage[pdftex]{graphicx}
\DeclareGraphicsExtensions{.jpg,.png,.pdf}
\DeclareGraphicsRule{*}{mps}{*}{}

% Some maths
%
\RequirePackage{amsmath,amsfonts,amssymb,amscd,amsthm,amsbsy}

\newcommand{\lifedraftdate}{\small Revision: \RCSRevision, \RCSDate, \RCSTime\ UTC Printed: \today}


%%
%% Index
%%
\newcommand{\ix}[1]{#1\index{#1}}
\newcommand{\nindex}[1]{\textit{#1}\index[name]{#1}}
\newcommand{\cppobj}[1]{\index{C++ Objects!#1}}

\makeglossary
\makeindex

\AtBeginDocument{

  \maketitle

  \ifthenelse{\boolean{life@stretch}}{\baselineskip 0.33in}{} 
  
  \ifthenelse{\boolean{life@toc}}{
    \setlength{\columnseprule}{0.6pt}
    \begin{multicols}{2}
      \tableofcontents
    \end{multicols}
  }     
  
}