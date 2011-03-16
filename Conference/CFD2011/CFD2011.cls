%-------------------------------------------------------------------------%This LaTeX class file was prepared to honor the author instructions
%given for
%the CFD2008 conference in Trondheim, june 2008.
%
%The class was created by
%Sverre Gullikstad Johnsen, sverre.g.johnsen@sintef.no
%SINTEF Materials and Chemistry
%Process technology, flow technology
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
\ProvidesClass{CFD2011}[2011/10/01 Prepared for the CFD2011 conference
in Trondheim]
\LoadClass[twocolumn, a4paper, 10pt, twoside]{article}


\newcommand{\headlinecolor}{\normalcolor}

\RequirePackage{color}
\RequirePackage{mathptmx} 
\RequirePackage{times}
\RequirePackage[T1]{fontenc}
\RequirePackage{helvet} 
\RequirePackage{fancyhdr}
\RequirePackage{natbib}
\renewcommand{\cite}{\citep}
%\RequirePackage{layout}


%% to get the awfull linespacing of word
\large
\renewcommand{\baselinestretch}{0.95}
\normalsize

\setlength{\voffset}{-23mm}
\setlength{\topmargin}{0mm}
\setlength{\headsep}{3.5mm}
\setlength{\headheight}{1cm}
\setlength{\topskip}{0mm}
\setlength{\footskip}{15mm}
\setlength{\textheight}{260.5mm}
\setlength{\parindent}{0pt}
\setlength{\textwidth}{180mm}
\setlength{\linewidth}{80mm}
\setlength{\columnsep}{10mm}
%\setlength{\leftmargin}{20mm}
%\setlength{\rightmargin}{20mm}
\setlength{\hoffset}{-0.5cm}
%\setlength{\marginparwidth}{20mm}
%\setlength{\oddsidemargin}{-5mm}
\setlength{\evensidemargin}{-5mm}

\newcommand{\arialEight}{\fontsize{8}{10}\fontfamily{phv}\fontseries{m}\fontshape{n}\selectfont}
\newcommand{\arialNine}{\fontsize{9}{11}\fontfamily{phv}\fontseries{m}\fontshape{n}\selectfont}
\newcommand{\arialNineBold}{\fontsize{9}{11}\fontfamily{phv}\fontseries{b}\fontshape{n}\selectfont}
\newcommand{\arialNineTitle}{\fontsize{9}{15}\fontfamily{phv}\fontseries{m}\fontshape{n}\selectfont}
\newcommand{\arialNineBoldTitle}{\fontsize{9}{15}\fontfamily{phv}\fontseries{b}\fontshape{n}\selectfont}
\newcommand{\arialTen}{\fontsize{10}{12}\fontfamily{phv}\fontseries{m}\fontshape{n}\selectfont}
\newcommand{\arialTenBold}{\fontsize{10}{12}\fontfamily{phv}\fontseries{b}\fontshape{n}\selectfont}
\newcommand{\arialTwelveBold}{\fontsize{12}{13}\fontfamily{phv}\fontseries{b}\fontshape{n}\selectfont}

%define fnsymbols
\def\fnsymbol#1{% 
  % #1 is between 1 and 9 inclusive 
  \ifcase#1\or 
  *\or\dag\or\ddag\or\S\or\P\or 
  $\Vert$\or**\or\dag\dag\or\ddag\ddag 
  \fi}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%macro to extract the first letter of a word
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcount\firstindex
\firstindex=0
\newcount\firstoutindex
\firstoutindex=0
\newcommand{\GetFirst}[1]{%
  \global\advance\firstindex by 1
  \expandafter\gdef\csname firstletter\romannumeral\firstindex\endcsname{#1}%
  \csname firstletter\romannumeral\firstindex\endcsname
}
%Output the letter stored in firstletter
\newcommand{\PutFirst}{
  \global\advance\firstoutindex by 1
  \csname firstletter\romannumeral\firstoutindex\endcsname
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define the title page information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\maketitle{\par
  \begingroup
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 1em\noindent
           \hb@xt@1.8em{%
                  \hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
    \twocolumn[\@maketitle]%
    \newpage
    \global\@topnum\z@  %Prevents figures from going at the top of
                        %page
    \thispagestyle{plain}
    %\thispagestyle{empty}
    \@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  %\global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{%
  \newpage
  \begin{center}%
  \let \footnote \thanks
    {\arialTwelveBold\headlinecolor
      \vskip 2.em
      \MakeUppercase{\@title}
    }%
    \vskip 1.em%
    {
      \begin{tabular}[t]{c}%
        \arialNineBold\headlinecolor
          % \lineskip .5em%
          % \lineskip .5em%
        \newlength{\templength}
        \begin{minipage}[n][-1\templength][t]{144mm}%
          \begin{centering}%
            \arialNineBold\headlinecolor
            {\@authoraddress}%
          \end{centering}%
        \end{minipage}
      \end{tabular}
    }%
    \arialNine\headlinecolor
    \vskip 5em%
  \end{center}%
  \vskip 9em
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   title, author and address block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\paperID}[1]{\gdef\@paperID{#1}}
\newtoks\@temptokenb
\long\def\append@item#1\to#2{%
  \@temptokena={{#1}}%
  \@temptokenb=\expandafter{#2}%
  \xdef#2{\the\@temptokenb\the\@temptokena}}
\def\@authors{}%
\def\author@#1{\append@item#1\to\@authors}%
\def\theauthors{
                \firstoutindex=0
                \@authors
              }
\def\thetitle{\@title}
\def\theID{\@paperID}
\newcount\temptotal
\newcount\authortotal
\newcount\authorindex
\authortotal=0
\newcount\addressnumberindex
\newcount\emailnumberindex
\newcount\emailtotal
\newcount\emailindex
\emailtotal=0
\renewcommand{\author}[2]{%
  \global\advance\authortotal by 1
  \expandafter\gdef\csname author\romannumeral\authortotal\endcsname
  {\ignorespaces\GetFirst#1\MakeUppercase{ #2}}%
  \ifnum\authortotal=1
    \author@{\PutFirst.  #2}
  \else
    \author@{, \PutFirst.  #2}
  \fi
  \expandafter\gdef\csname addressbool\romannumeral\authortotal\endcsname{0}%
  \expandafter\gdef\csname emailbool\romannumeral\authortotal\endcsname{0}%
  \expandafter\gdef\csname speakbool\romannumeral\authortotal\endcsname{0}%
}
\newcommand{\presenting}{
  \expandafter\gdef\csname speakbool\romannumeral\authortotal\endcsname{1}%
}
\newcount\addresstotal
\newcount\addressindex
\addresstotal=0
\newcommand{\address}[1]{%
  \expandafter\gdef\csname address\romannumeral\authortotal\endcsname
  {\ignorespaces#1}%
   \expandafter\gdef\csname addressbool\romannumeral\authortotal\endcsname{1}%
}
\newcommand{\email}[1]{%
  \expandafter\gdef\csname email\romannumeral\authortotal\endcsname
  {\ignorespaces\MakeLowercase{#1}}%
   \expandafter\gdef\csname emailbool\romannumeral\authortotal\endcsname{1}%
}
\newcommand{\@authoraddress}{%
  \newcommand{\@putauthor}{\csname author\romannumeral\authorindex\endcsname}%
  \newcommand{\@putspeakbool}{\csname speakbool\romannumeral\authorindex\endcsname}
  \newcommand{\@putaddress}{\csname address\romannumeral\authorindex\endcsname}%
  \newcommand{\@putaddressbool}{\csname addressbool\romannumeral\authorindex\endcsname}%
  \newcommand{\@putemail}{\csname email\romannumeral\authorindex\endcsname}%
  \newcommand{\@putemailbool}{\csname emailbool\romannumeral\authorindex\endcsname}%
  \global\addressindex=0
  \global\authorindex=0
  \global\emailindex=0
  \global\temptotal=\authortotal
  \ifnum\temptotal>1 \advance\temptotal by -1 \fi
  \loop \ifnum\authorindex<\temptotal \global\advance\authorindex by
  1%
    %% add a comma behind the author if there are any left
    \ifnum\authorindex=0\else{, }\fi%
    \arialNineBoldTitle
    \ifodd\@putspeakbool
      \underline\@putauthor
    \else
      \@putauthor%
    \fi
    \ifnum\authortotal>1
      \ifodd\@putaddressbool
        \global\advance\addressindex by 1%
        ${}^{\number\addressindex}$%
      \fi
      \ifodd\@putemailbool
        \global\advance\emailindex by 1%
        ${}^{\fnsymbol\emailindex}$%
      \fi
    \fi
  \repeat%
  %% instead of the comma add the word "and" before the last author 
  \ifnum\authorindex<\authortotal
    \global\advance\authorindex by 1%
    %%{ and }%
    {, }%
    \ifodd\@putspeakbool
      \underline\@putauthor
    \else
      \@putauthor%
    \fi
    %% number is necessary only for more than one author
    \ifnum\authortotal>1
      \ifodd\@putaddressbool
        \global\advance\addressindex by 1%
        ${}^{\number\addressindex}$%
      \fi
      \ifodd\@putemailbool
        \global\advance\emailindex by 1%
        ${}^{\fnsymbol\emailindex}$%
      \fi
    \fi
  \fi%
  \par%
  \global\authorindex=0
  \global\addressindex=0
  \arialNineTitle
  \loop \ifnum\authorindex<\authortotal\global\advance\authorindex by
  1%
    \par
    % number is necessary only for more than one address
    \ifnum\authortotal>1
      \par
      \ifodd\@putaddressbool
        \global\advance\addressindex by 1%
        ${}^{\number\addressindex}$%
      \fi
    \fi
    \@putaddress
    \par%
  \repeat%
  \vskip 1em%
  \global\emailindex=0%
  \global\authorindex=0%
  \arialEight
   \loop \ifnum\authorindex<\authortotal\global\advance\authorindex by 1%
      \par 
      \ifodd\@putemailbool
        \par
         % number is necessary only for more than one author
         \ifnum\authortotal>1
           \global\advance\emailindex by 1%
           ${}^{\fnsymbol\emailindex}$
         \fi
         E-mail: \@putemail
       \fi
     \par%
   \repeat%
  \global\emailtotal=0%
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define headers/footers appearence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\headers}{
  \pagestyle{fancyplain}
  \renewcommand{\headrulewidth}{0pt}
  \renewcommand{\footrulewidth}{0pt}
  \cfoot{\fontsize{9}{9}\fontfamily{times}\fontseries{m}\selectfont
    \headlinecolor\thepage}
  \rhead[\fancyplain{}{}]%
  {\fancyplain{\fontsize{9}{9}\selectfont\headlinecolor\it\flushright{\theID\vskip .2em}}
    {}
  }%
  \lhead[\fancyplain{}{}]%
  {\fancyplain{\fontsize{9}{9}\selectfont\headlinecolor\it
      \quad 8th International Conference on CFD in Oil \& Gas, Metallurgical and Process Industries\\
      \quad SINTEF/NTNU, Trondheim NORWAY\\
      \quad 21-23 June 2011}
    {}
  }%
  \chead[\fancyplain{}{\fontsize{9}{9}\fontfamily{times}\fontseries{m}
    \selectfont\headlinecolor\it{\theauthors}\vskip 2em}
  ]
  {\fancyplain{}
    {\fontsize{9}{9}\fontfamily{times}\fontseries{m}
      \selectfont\headlinecolor\it{\thetitle / \theID\vskip 2em}
    }
  }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define heading appearence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\setcounter{secnumdepth}{0}
\renewcommand{\section}{%
  \@startsection
     {section}{1}{0pt}{-1.5ex plus -1ex minus -.2ex}{1ex plus .2ex}
     {\arialTenBold\headlinecolor\selectfont\MakeUppercase}}%
\renewcommand{\subsection}{%
  \@startsection
     {subsection}{2}{0pt}{-1.5ex plus -1ex minus -.2ex}{1ex plus .2ex}
     {\arialTenBold\headlinecolor\selectfont}}%
\renewcommand{\subsubsection}{%
  \@startsection
     {subsection}{3}{0pt}{-1.5ex plus -1ex minus -.2ex}{1ex plus .2ex}
     {\fontsize{10}{12}\fontfamily{times}\fontseries{m}\it
       \selectfont\headlinecolor\selectfont}}%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define abstract and keyword macros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand{\abstract}[1]{
  \arialTenBold\headlinecolor\selectfont\MakeUppercase{Abstract}\\
   \normalfont\fontsize{9}{11}\selectfont
   #1
   \normalsize
}
\newcommand{\keywords}[1]{
  \vskip 1.em
  \arialNineBold
  Keywords:   
  \normalfont\fontsize{9}{11}\selectfont
  #1.
  \normalsize
  \vskip 3.em
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%modify the bibliography appearence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\NAT@bibsetup%
   [1]{\setlength{\leftmargin}{0pt}\setlength{\itemindent}{1em}%
       \setlength{\itemsep}{0pt}\setlength{\parsep}{\z@}}

%The following lines were copied from phiaea.bst and pasted into the
%CFD2008.bst
%at line475 , to get capitalized author names.
%       t "others" =
% 	'skip$
%       { t "u" change.case$ 't := }	
% 	if$

