%%%
%%% Students: if you're making modifications to this file, you'd
%%% better know what you're doing. If everything is done right, you
%%% should *only* be making changes to your own thesis.tex and other
%%% contents of your thesis or dissertation. This file is most
%%% definitely *not* where you make changes to the names of your
%%% committee members or other things that are specific to your
%%% particular thesis.
%%%
%%% See Mike Renfro (email: renfro, office: Clement 224, phone: x3601)
%%% or one of the LaTeX-conversant faculty if you run into problems with
%%% this style as presented, or if the graduate school flags problems
%%% with your formatting. If you're doing something wrong, we need to
%%% figure out how to do it right. If this document class is wrong, it
%%% needs to be fixed for *everyone*, not just you.
%%%

\NeedsTeXFormat{LaTeX2e}[1994/06/01]
\ProvidesClass{ttuthesis}[2017/12/30 v1.0 Tennessee Tech Thesis/Dissertation Class]
\newcount\copyrightpagepresent \copyrightpagepresent=0
\DeclareOption{copyrighted}{\copyrightpagepresent=1}
\DeclareOption{oneside}{%
  \ClassWarning{ttuthesis}{oneside option is not suitable for final prints through ProQuest}%
  \PassOptionsToClass{\CurrentOption}{memoir}%
  }
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{memoir}}
\ProcessOptions\relax
\LoadClass{memoir}

% The ted package for patching the chapter command and PDF title
\RequirePackage{ted}

% The ifthen package provides some clean if/then/else structures, and works
% with PCTeX.
\RequirePackage{ifthen}

% The indent package indents all first lines of paragraphs, including
% those that occur right after a section or chapter header.
\RequirePackage{indentfirst}
% rotating: for rotated tables or figures
\RequirePackage{rotating}
% For including code listings
\RequirePackage{listings}

% For checking if ToC, LoT, LoF, etc. span multiple pages
\RequirePackage{refcount}
\RequirePackage{labelcas} % obsolete on CTAN, included in thesis achives as of 2017/12/30.

% For automatic lists of symbols
\RequirePackage[intoc]{nomencl} \makenomenclature

% For proper uppercasing of headings containing math (e.g.,
% \chapter{The $\phi_{ij}$ Matrix})
\usepackage[overload]{textcase}

%%%
%%% Define paper size and margins for all pages
%%%
\setstocksize{11in}{8.5in}
\settrims{0pt}{0pt}
\settypeblocksize{9in}{6in}{*}
\setlrmarginsandblock{1.5in}{1in}{*}
\setulmarginsandblock{1in}{1in}{*}
\addtolength{\uppermargin}{4\onelineskip}
\addtolength{\textheight}{-6\onelineskip}
\setheadfoot{3\baselineskip}{\footskip}
\setheaderspaces{1in}{*}{*}

%%% Define default text for elements that the student should fill out
%%% themselves.

%%% Note to students: if you're thinking of changing any of this in
%%% the ttuthesis.cls file, please change it in your thesis.tex
%%% instead. Something like:
%%%
%%% \abstract{%
%%%   This is my first paragraph.
%%%
%%%   And this is my second one.}
\newcommand{\doctype}[1]{\providecommand\ttu@doctype{}\renewcommand{\ttu@doctype}{#1}}
\doctype{(Is this a thesis or a dissertation? Fix your thesis.tex.)}

% Can't follow the same pattern with \abstract, since that's already defined
% in memoir.cls
\let\ttuabstract\abstract
\let\endttuabstract\endabstract
\renewcommand{\abstract}[1]{\providecommand\ttu@abstract{}\renewcommand{\ttu@abstract}{#1}}
\abstract{%
  If you can read this, you forgot to put an abstract in your thesis.tex file.
}

\newcommand{\degree}[1]{\providecommand\ttu@degree{}\renewcommand{\ttu@degree}{#1}}
\degree{(What's your degree? Fix your thesis.tex.)}
\newcommand{\department}[1]{\providecommand\ttu@department{}\renewcommand{\ttu@department}{#1}}
\department{(What's your department? Fix your thesis.tex.)}
\newcommand{\graduationmonth}[1]{\providecommand\ttu@month{}\renewcommand{\ttu@month}{#1}}
\graduationmonth{(What month? Fix your thesis.tex.)}
\newcommand{\graduationyear}[1]{\providecommand\ttu@year{}\renewcommand{\ttu@year}{#1}}
\graduationyear{(What month? Fix your thesis.tex.)}
\newcommand{\dedication}[1]{\providecommand\ttu@dedication{}\renewcommand{\ttu@dedication}{#1}}
%\dedication{%
%  If you can read this, you asked for a dedication page to be
%  printed, but failed to include any dedication text in your
%  thesis file. Fix your thesis.tex.
%}
\newcommand{\acknowledgments}[1]{\providecommand\ttu@acknowledgments{}\renewcommand{\ttu@acknowledgments}{#1}}
%\acknowledgments{%
%  If you can read this, you asked for an acknowledgments page to be
%  printed, but failed to include any acknowledgments text in your
%  thesis file. Fix your thesis.tex.
%}

%%%
%%% Define style for abstract page.
%%%
\renewcommand{\abstractnamefont}{\normalfont\bfseries}
\renewcommand{\abstracttextfont}{\normalfont}

\makepagestyle{ttuabstract}

\newcommand{\abstractpage}{%
  \ifdraftdoc
    \thispagestyle{empty}
    \begin{center}
      {\bfseries\MakeUppercase{\thetitle}}
    \end{center}

    \begin{center}
      \theauthor
    \end{center}

    \begin{center}
      Draft created \today{} at \printtime*
    \end{center}

    \begin{center}
      \begin{ttuabstract}
        \addtolength{\parskip}{\baselineskip}
        \thispagestyle{ttuabstract}
        \ttu@abstract
        \addtolength{\parskip}{-\baselineskip}
      \end{ttuabstract}
    \end{center}
  \else
    \pagenumbering{alph}
	  \setcounter{page}{1}
    \thispagestyle{empty}
    \renewcommand{\baselinestretch}{1}
    \begin{center}
      {\bfseries\MakeUppercase{An Abstract of a \ttu@doctype}}
      
      \vspace*{\onelineskip}
      
      \begin{Spacing}{1}
        {\bfseries\MakeUppercase{\thetitle}}
      \end{Spacing}
      
      \vspace*{3\onelineskip}
      
      {\theauthor}
      
      \vspace*{\onelineskip}
      
      {\ttu@degree{} in \ttu@department}
    \end{center}
    \renewcommand{\abstractname}{}
    \vspace*{-2\onelineskip}
    \begin{ttuabstract}
      \addtolength{\parskip}{\baselineskip}
      \thispagestyle{ttuabstract}
      \ttu@abstract
      \addtolength{\parskip}{-\baselineskip}
    \end{ttuabstract}
    \renewcommand{\baselinestretch}{2}
  	\newpage
  \fi
  \frontmatter
}

%%%
%%% Define style for title page.
%%%

\newcommand{\titlepage}{%
  \ifdraftdoc
  \begin{center}
  (Draft copy: normal title page suppressed)
  \end{center}
  \else
    \pagenumbering{roman}
    \setcounter{page}{1}
    \thispagestyle{empty}
    \beforepartskip
    \begin{center}
      {\bfseries\MakeUppercase{\thetitle}}
    \end{center}
    \hfil \hbox{\vrule width 1.75in height 0.02in} \hfil
    \begin{center}
      A \ttu@doctype \\
      Presented to \\
      the Faculty of the College of Graduate Studies \\
      Tennessee Technological University \\
      by \\
      \theauthor
    \end{center}
    \hfil \hbox{\vrule width 1.75in height 0.02in} \hfil
    \begin{center}
      In Partial Fulfillment \\
      of the Requirements for the Degree \\
      \MakeUppercase{\ttu@degree} \\
      \ttu@department
    \end{center}
    \hfil \hbox{\vrule width 1.75in height 0.02in} \hfil
    \begin{center}
      \ttu@month{} \ttu@year
    \end{center}
    \afterpartskip
  	\newpage
  \fi
}

%%%
%%% Define style for copyright page, if any.
%%%
%\newcount\copyrightpagepresent \copyrightpagepresent=0
\newcommand{\copyrightpage}{%
  %\copyrightpagepresent=1
  \ifdraftdoc
    \begin{center}
      (Draft copy: normal copyright page suppressed)
    \end{center}
  \else
    \thispagestyle{plain}
    \beforepartskip
    \begin{center}
      \begin{Spacing}{1}
        Copyright \copyright{} \theauthor, \ttu@year \\
        All rights reserved
      \end{Spacing}
    \end{center}
    \afterpartskip
    \newpage
  \fi
}

%%%
%%% Define style for committee approval page.
%%%

% Define default value for presence of cochair. This is set to the
% proper value automatically later.
\newcount\cochair \cochair=0
\newcommand{\committeechair}[1]{\providecommand\ttu@committeechair{}\renewcommand{\ttu@committeechair}{#1}}
\committeechair{Chair's Name\\(fix your thesis.tex)}
\newcommand{\committeecochair}[1]{\providecommand\ttu@committeecochair{}\renewcommand{\ttu@committeecochair}{#1}}
\committeecochair{Cochair's Name\\(fix your thesis.tex)}
\newcommand{\committeemembers}[1]{\providecommand\ttu@committeemembers{}\renewcommand{\ttu@committeemembers}{#1}}
\committeemembers{
Member 1 Name\\(fix your thesis.tex),
Member 2 Name\\(fix your thesis.tex),
Member 3 Name\\(fix your thesis.tex)}

\newcommand{\nameanddateline}[2]{%
\begin{Spacing}{1}
\par \noindent%
\parbox{2.25in}{\hrulefill}\hspace{.75in}\parbox{1.25in}{\hrulefill} \\
\vspace{-\onelineskip}
\par \noindent %
\parbox{2.25in}{\raggedright #1\ifthenelse{\equal{#2}{}}{\relax}{, #2}}%
\hspace{.75in}\makebox[1.25in][l]{Date} \\
\vspace{\onelineskip}
\end{Spacing}
}

\newcommand{\vpname}{Mark Stephens}
\newcommand{\vptitle}{Dean}
\newcommand{\vpdept}{College of Graduate Studies}
\newcommand{\approvalpage}{%
  \ifdraftdoc
    \begin{center}
      (Draft copy: normal approval page suppressed)
    \end{center}
  \else
    \vspace*{-7\onelineskip}
		\ifthenelse{\equal{\ttu@committeecochair}{Cochair's Name\\(fix your thesis.tex)}}{
		  \newcommand{\chairtitle}{Chairperson}
			\cochair=0
		}{
      \newcommand{\chairtitle}{Cochairperson}
			\cochair=1
		}
    \thispagestyle{plain}
    \begin{center}
      {\bfseries\MakeUppercase{Certificate of Approval of \ttu@doctype}}

      \vspace{\onelineskip}

      {\bfseries{\MakeUppercase{\thetitle}}} \\
      by \\
      \theauthor
    \end{center}

    \vfill

    \vspace*{-\onelineskip}
    \noindent Graduate Advisory Committee:
    \vspace*{\onelineskip}

    \nameanddateline{\ttu@committeechair}{\chairtitle}
    \ifcase\cochair
      \relax % \cochair=0
    \or
      \nameanddateline{\ttu@committeecochair}{\chairtitle} % \cochair=1
    \else
      \relax
    \fi
    \@for\themember:=\ttu@committeemembers\do{\nameanddateline{\themember}{}}
    \vfill
    \newlength{\vpindent}
    \setlength{\vpindent}{3.0in}
    \begin{Spacing}{1}
      \par\noindent\hspace{\vpindent}\parbox{2.25in}{Approved for the Faculty:}\\
      \mbox{} %\\
      \par\noindent\hspace{\vpindent}\parbox{2.25in}{\hrulefill}
      \par\noindent\hspace{\vpindent}\parbox{2.25in}{\vpname \\ \vptitle \\ \vpdept} \\ \par 
      \mbox{} \par
      \par\noindent\hspace{\vpindent}\parbox{2.25in}{\hrulefill}
      \par\noindent\hspace{\vpindent}\parbox{2.25in}{Date}
    \end{Spacing}
    \vfill

    \newpage
  \fi
}

%%%
%%% Define style for dedication page.
%%%
\newcommand{\dedicationpage}{%
  \ifdraftdoc
    \begin{center}
      {\bfseries\MakeUppercase{Dedication}}
    \end{center}
    \ttu@dedication
  \else
    \thispagestyle{plain}
    \beforepartskip
    \begin{center}
      {\bfseries\MakeUppercase{Dedication}}
    \end{center}
    \ttu@dedication
    \afterpartskip
  	\newpage
  \fi
}

%%%
%%% Define style for acknowledgments page.
%%%
\newcommand{\acknowledgmentspage}{%
  \ifdraftdoc
  \begin{center}
    {\bfseries\MakeUppercase{Acknowledgments}}
  \end{center}

    \ttu@acknowledgments
  \else
		\thispagestyle{plain}
    \begin{center}
      {\bfseries\MakeUppercase{Acknowledgments}}
    \end{center}
    \vspace*{\onelineskip}

    \ttu@acknowledgments
		\newpage
  \fi
}

%%%
%%% Define Table of Contents style
%%%
\maxsecnumdepth{subsubsection}
\settocdepth{subsubsection}
\copypagestyle{toc-initial}{plain}
\makepagestyle{toc}
\makeoddhead{toc}{\chaptername}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{%
  \enlargethispage{-3\onelineskip}\normalfont\thepage\par%
  \vspace{\onelineskip}Page}
\makeevenhead{toc}{%
  \enlargethispage{-3\onelineskip}\normalfont\thepage\par%
  \vspace{\onelineskip}\chaptername}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{Page}
\makeoddfoot{toc}{}{}{}
\makeevenfoot{toc}{}{}{}
\makepagestyle{toc-app}
\makeoddhead{toc-app}{\appendixname}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{%
  \normalfont\thepage\par%
  \vspace{\onelineskip}Page}
\makeevenhead{toc-app}{%
  \normalfont\thepage\par%
  \vspace{\onelineskip}\appendixname}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{Page}
\makeoddfoot{toc-app}{}{}{}
\makeevenfoot{toc-app}{}{}{}

\renewcommand{\contentsname}{Table of Contents}
\renewcommand{\printtoctitle}[1]{%
  \centering \normalfont\bfseries\MakeUppercase{#1}}
\renewcommand{\aftertoctitle}{%
  \par\nobreak\vspace{\onelineskip}\vspace{\onelineskip}
  \mbox{}\hfill{\normalfont Page}
  \par\nobreak
  % The listings package has a bad habit of running this command, causing the
  % misc:beforetoc label to be defined twice. Fix given in 
  % http://stackoverflow.com/questions/495841/is-there-a-way-to-check-if-a-label-is-already-defined-in-latex
  %
  % If counter for misc:beforetoc is not defined, create it and set it to zero
  \protect\@ifundefined{c@misc:beforetoc}{%
    \newcounter{misc:beforetoc}%
    \setcounter{misc:beforetoc}{0}}{}%
  % If the misc:beforetoc counter is zero, add a beforetoc label to the ToC,
  % and increment the counter. Otherwise, do nothing.
  \ifthenelse{\value{misc:beforetoc} > 0}{}{%
    \addtocontents{toc}{\protect\label{misc:beforetoc}}%
    \addtocounter{misc:beforetoc}{1}%
  }%
  \addtocontents{toc}{\protect\pagestyle{toc}%
    \protect\thispagestyle{toc-initial}}%
}
\maxtocdepth{subsubsection}
\setlength{\cftbeforechapterskip}{2\baselineskip}

\renewcommand{\addappheadtotoc}{%
}
\newcount\appendicesonfirstpageoftoc
\renewcommand{\appendixpagename}{\MakeUppercase{\appendixtocname}}
\addtodef{\appendix}{}{\appendixpage \chapterstyle{appendix}}
% Ensuring \tableofcontents* automatically runs \checktocheaders
% https://groups.google.com/group/comp.text.tex/msg/fa73a51a04884301?hl=en
\let\oldmem@tableofcontents\mem@tableofcontents
\renewcommand\mem@tableofcontents[1]{%
\oldmem@tableofcontents{#1}%
\checktocheaders
}

% Put ``Chapter'' ahead of Chapter 1 entry in TOC, and ``Appendix''
% ahead of Appendix A entry in TOC. This is a horrible abuse of a
% function that originally was just used to break up each chapter in
% the List of Tables, Figures, etc. But it's the only function I can
% find that is called for each chapter, before the chapter entry is
% added to the TOC.
\renewcommand{\insertchapterspace}{%
  \ifanappendix%
    \ifthenelse{\equal{\thechapter}{A}}{%
      \addtocontents{toc}{\protect\label{misc:beforeapp}}%
      \ifcase\appendicesonfirstpageoftoc%
        %\typeout{appendicesonfirstpageoftoc = 0}%
        % Appendix listings start on second-Nth page of ToC
        \addtocontents{toc}{\protect\vskip\baselineskip\noindent\protect\mbox{\MakeUppercase{\appendixtocname}}\protect\hfill\par}%
      \or%
        %\typeout{appendicesonfirstpageoftoc = 1}%
        % Appendix listings start on first page of ToC
        \addtocontents{toc}{\protect\vskip\baselineskip\noindent\protect\mbox{\MakeUppercase{\appendixtocname}}\protect\hfill\par}%
      \else%
        \typeout{*** Error: appendicesonfirstpageoftoc must equal 0 or 1}
      \fi%
    }
  \else % not an appendix, i.e., it's a regular chapter
    \ifthenelse{\equal{\thechapter}{1}}{%
        \addtocontents{toc}{\protect\vskip\baselineskip\noindent\protect\mbox{\chaptername}\protect\hfill\par}%
    }{%
    {}%
    }
  \fi%
\addtocontents{toc}{\protect\vskip-\baselineskip}
}

% Normal, non-bolded fonts for chapters and page numbers.
\renewcommand{\cftchapterfont}{\normalfont}
\renewcommand{\cftchapterpagefont}{\normalfont}

% Put a period after the chapter number.
\renewcommand{\cftchapteraftersnum}{. }

%
% This is a hack to help ensure that all chapter titles are
% put in all-caps in the ToC, and that all appendix titles are not put
% in all-caps in the ToC. But at least I don't have a complete copy of
% memoir's \chapter command anymore with one edit applied.
%
% http://tex.stackexchange.com/questions/13248/using-ted-or-another-package-instead-of-etoolbox-for-patching-commands-in-memo
% http://tex.stackexchange.com/questions/42181/stop-memoir-from-making-mathmode-section-names-uppercase-in-header
\Substitute*[{\gdef\@chapter[##1]##2}]{\@chapter[#1]{#2}}%
{\addcontentsline{toc}{chapter}}%
{\let\f@rtocold\f@rtoc%
 \def\f@rtoc{\texorpdfstring{\MakeTextUppercase{\f@rtocold}}{\f@rtocold}}%
 \addcontentsline{toc}{chapter}}

% Dot leaders on chapter entries, identical to the other headings.
\renewcommand{\cftchapterdotsep}{\cftsectiondotsep}
\renewcommand{\cftchapterleader}{\normalfont\cftdotfill{\cftchapterdotsep}}

% Disable numbering of divisions lower than chapters.
\let\oldcftsf\cftsectionfont
\let\oldcftspn\cftsectionafterpnum

% Indent each level of heading by 0.25in (for 12pt, works out to 1.5em)
\cftsetindents{chapter}{0em}{1.5em}
\cftsetindents{section}{1.50em}{2.25em}
\cftsetindents{subsection}{3.0em}{3.0em}
\cftsetindents{subsubsection}{4.5em}{3.75em}

%%%
%%% Define List of Tables style
%%%
\renewcommand{\listtablename}{LIST OF TABLES}
\copypagestyle{lot-initial}{plain}
\makepagestyle{lot}
\makeoddhead{lot}{\tablename}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{\enlargethispage{-3\onelineskip}\normalfont\thepage\par\vspace{\onelineskip}Page}
\makeevenhead{lot}{\enlargethispage{-3\onelineskip}\normalfont\thepage\par\vspace{\onelineskip}\tablename}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{Page}
\makeoddfoot{lot}{}{}{}
\makeevenfoot{lot}{}{}{}
\renewcommand{\lotheadstart}{\newpage \chapterheadstart}
\renewcommand{\afterlottitle}{%
  \let\numberline\oldnumberline
  \par\nobreak\vspace{\onelineskip}\vspace{\onelineskip}
  {\normalfont Table}\hfill{\normalfont Page}
  \par\nobreak\vspace{\onelineskip}
  \addtocontents{lot}{\protect\label{misc:beforelot}\protect\pagestyle{lot}\protect\thispagestyle{lot-initial}}
}
\cftsetindents{table}{0in}{2.4em}
\addtodef{\listoftables}{}{\checklotheaders}

%%%
%%% Define List of Figures style
%%%
\renewcommand{\listfigurename}{LIST OF FIGURES}
\copypagestyle{lof-initial}{plain}
\makepagestyle{lof}
\makeoddhead{lof}{\figurename}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{\enlargethispage{-3\onelineskip}\normalfont\thepage\par\vspace{\onelineskip}Page}
\makeevenhead{lof}{\enlargethispage{-3\onelineskip}\normalfont\thepage\par\vspace{\onelineskip}\figurename}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{Page}
\makeoddfoot{lof}{}{}{}
\makeevenfoot{lof}{}{}{}
\renewcommand{\lofheadstart}{\newpage \chapterheadstart}
\renewcommand{\afterloftitle}{%
  \par\nobreak\vspace{\onelineskip}\vspace{\onelineskip}
  {\normalfont Figure}\hfill{\normalfont Page}
  \par\nobreak\vspace{\onelineskip}
  \addtocontents{lof}{\protect\label{misc:beforelof}\protect\pagestyle{lof}\protect\thispagestyle{lof-initial}}
}
\cftsetindents{figure}{0in}{2.4em}
\addtodef{\listoffigures}{}{\checklofheaders}

%%%
%%% Define List of Symbols style
%%%
\copypagestyle{los-initial}{plain}
\makepagestyle{los}
\makeoddhead{los}{}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{\thepage}
\makeevenhead{los}{\thepage}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{}
\makeoddfoot{los}{}{}{}
\makeevenfoot{los}{}{}{}
\renewcommand{\nompreamble}{\protect\pagestyle{los}\protect\thispagestyle{los-initial}}

%%%
%%% Define List of Programs style
%\renewcommand{\lstlistlistingname}{LIST OF PROGRAMS}
%\renewcommand{\lstlistingname}{Program}
\newlistof{listoflistings}{lol}{\lstlistlistingname}
\newlistentry[chapter]{listing}{lol}{0}
\copypagestyle{lol-initial}{plain}
\makepagestyle{lol}
\makeoddhead{lol}{\lstlistingname}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{\enlargethispage{-3\onelineskip}\normalfont\thepage\par\vspace{\onelineskip}Page}
\makeevenhead{lol}{\lstlistingname}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{\enlargethispage{-3\onelineskip}\normalfont\thepage\par\vspace{\onelineskip}Page}
\makeoddfoot{lol}{}{}{}
\makeevenfoot{lol}{}{}{}
\renewcommand{\lolheadstart}{\newpage \chapterheadstart}
\renewcommand{\afterloltitle}{%
  \par\nobreak\vspace{\onelineskip}\vspace{\onelineskip}
  {\normalfont \lstlistingname}\hfill{\normalfont Page}
  \par\nobreak\vspace{\onelineskip}
  \addtocontents{lol}{\protect\label{misc:beforelol}\protect\pagestyle{lol}\protect\thispagestyle{lol-initial}}
}
\renewcommand*{\l@lstlisting}[2]{\@dottedtocline{0}{0em}{2.4em}{#1}{#2}}
\let\lstlistoflistings\listoflistings
\addtodef{\lstlistoflistings}{}{\checklolheaders}

%%%
%%% Final cleanup commands for the ToC, LoT, and LoF.
%%%
\newcommand{\toconmultiplepagesnoappendices}{%
  \addtocontents{toc}{\protect\thispagestyle{toc}}
}
\newcommand{\toconmultiplepagesappendicesononepage}{%
  \addtocontents{toc}{\protect\thispagestyle{toc}}
}
\newcommand{\appendixlistonmultiplepages}{%
  \addtocontents{toc}{\protect\thispagestyle{toc-app}}
}
\newcommand{\tablelistonmultiplepages}{%
  \addtocontents{lot}{\protect\thispagestyle{lot}}
}
\newcommand{\figurelistonmultiplepages}{%
  \addtocontents{lof}{\protect\thispagestyle{lof}}
}
\newcommand{\nomenclatureonmultiplepages}{%
  \renewcommand{\nompostamble}{\protect\thispagestyle{los}}
}

\newcounter{tocmultipage} \setcounter{tocmultipage}{0}
\newcounter{appmultipage} \setcounter{appmultipage}{0}
\newcommand{\checktocheaders}%
{%
  \ifonsamepage{misc:beforetoc}%
  {%
    \setcounter{tocmultipage}{0}%
    \typeout{Single-ToC (says checktocheaders)}%
  }%
  {%
    \setcounter{tocmultipage}{1}%
    \eachlabelcase{%
      {misc:beforeapp}%
      {%
        \ifonsamepage{misc:beforeapp}%
        {%
          \setcounter{appmultipage}{0}%
          \typeout{Multi-ToC, Single-App (says checktocheaders)}%
        }%
        {%
          \setcounter{appmultipage}{1}%
          \typeout{Multi-ToC, Multi-App (says checktocheaders)}%
        }%
      }%
      {%
        \typeout{Multi-ToC (No-App) (says checktocheaders)}%
      }%
    }%
  }%
}%

\newcommand{\fixtocheaders}{
  \ifthenelse{\equal{\value{tocmultipage}}{1}}%
  {%
    \ifthenelse{\equal{\value{appmultipage}}{1}}%
    {%
      \immediate
			\write\@auxout{\noexpand\@writefile{toc}{\noexpand\thispagestyle{toc-app}}}
      \typeout{Multi-ToC, Multi-App (says fixtocheaders)}%
    }%
    {%
      \immediate
			\write\@auxout{\noexpand\@writefile{toc}{\noexpand\thispagestyle{toc}}}
      \typeout{Multi-ToC, Single-App or No-App (says fixtocheaders)}%
    }%
  }{%
    \immediate
		\write\@auxout{\noexpand\@writefile{toc}{\noexpand\thispagestyle{toc-initial}}}
    \typeout{Single-ToC (says fixtocheaders)}%
  }%
}

\newcounter{lofmultipage} \setcounter{lofmultipage}{0}
\newcommand{\checklofheaders}{%
  \ifonsamepage{misc:beforelof}{
	}{%
    \setcounter{lofmultipage}{1}%
    \typeout{Multi-LoF (says checklofheaders)}%
  }%
}
\newcommand{\fixlofheaders}{
  \ifthenelse{\equal{\value{lofmultipage}}{1}}{
	\immediate
	\write\@auxout{\noexpand\@writefile{lof}{\noexpand\thispagestyle{lof}}}
	}{
	\immediate
	\write\@auxout{\noexpand\@writefile{lof}{\noexpand\thispagestyle{lof-initial}}}
	}
}

\newcounter{lotmultipage} \setcounter{lotmultipage}{0}
\newcommand{\checklotheaders}{%
  \ifonsamepage{misc:beforelot}{
	}{%
    \setcounter{lotmultipage}{1}%
    \typeout{Multi-LoT (says checklotheaders)}%
  }%
}
\newcommand{\fixlotheaders}{
  \ifthenelse{\equal{\value{lotmultipage}}{1}}%
  {\immediate
	\write\@auxout{\noexpand\@writefile{lot}{\noexpand\thispagestyle{lot}}}
	}{
	\write\@auxout{\noexpand\@writefile{lot}{\noexpand\thispagestyle{lot-initial}}}
	}
}

\newcounter{lolmultipage} \setcounter{lolmultipage}{0}
\newcommand{\checklolheaders}{%
  \ifonsamepage{misc:beforelol}{
	}{%
    \setcounter{lolmultipage}{1}%
    \typeout{Multi-LoL (says checklolheaders)}%
  }%
}
\newcommand{\fixlolheaders}{
  \ifthenelse{\equal{\value{lolmultipage}}{1}}{
	\immediate
	\write\@auxout{\noexpand\@writefile{lol}{\noexpand\thispagestyle{lol}}}
	}{
	\immediate
	\write\@auxout{\noexpand\@writefile{lol}{\noexpand\thispagestyle{lol-initial}}}
	}
}

\newcounter{cPage}
\newcommand\ifonsamepage[3]{%
 \refstepcounter{cPage}\label{current\thecPage}%
 \ifthenelse{\equal{\pageref{#1}}{\pageref{current\thecPage}}}%
  {\protect{#2}}%
  {\protect{#3}}%
  %\typeout{pageref 1 = \pageref{#1}, pageref current = \pageref{current\thecPage}}%
}

% Ensure ToC headers are corrected
\AtEndDocument{\fixtocheaders
\fixlofheaders
\fixlotheaders
\fixlolheaders
}

%%%
%%% Body text formatting
%%%

\addtodef{\mainmatter}{\pagestyle{thesis}}{}
% Double-space text, set default indentation to 3.0em (0.5in for 12pt).
\ifdraftdoc\else
\renewcommand{\baselinestretch}{2}
\fi
\parindent=3.0em

% Define headers and footers for each page
\copypagestyle{thesis-initial}{plain}
\copypagestyle{thesis}{myheadings}
\makeoddhead{thesis}{}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{\normalfont\thepage}
\makeevenhead{thesis}{\normalfont\thepage}{\ifdraftdoc Draft (\theauthor): \today{} at \printtime*\else\fi}{}

% Define styles for chapter headings
\makechapterstyle{default}{%
\renewcommand{\chapnamefont}{\normalfont\bfseries}
\renewcommand{\chapnumfont}{\normalfont\bfseries}
\renewcommand{\chaptitlefont}{\normalfont\bfseries}

% Remove extra space at top of chapter heading
\renewcommand{\chapterheadstart}{}
% Print ``CHAPTER'' in all caps, centered
\renewcommand{\printchaptername}{%
  \centering \chapnamefont \MakeUppercase{\@chapapp}
}
\renewcommand{\printchapternum}{%
  \centering \chapnumfont \thechapter\space
}
\setlength{\midchapskip}{0pt}
\renewcommand{\printchaptertitle}[1]{%
  \centering \chaptitlefont \MakeUppercase{##1}
}
\setlength{\afterchapskip}{2\baselineskip}
}
\renewcommand{\parttitlefont}{\normalfont\bfseries} % For Appendix Page

% Define styles for section headings (centered head)
\setsecheadstyle{\normalfont\bfseries\raggedright\centering}
\setbeforesecskip{\baselineskip}

% Define styles for subsection headings (freestanding sidehead)
\setsubsecheadstyle{\bfseries\raggedright}
\setbeforesubsecskip{\baselineskip}

% Define styles for subsubsection headings (paragraph sidehead)
\ifthenelse{\isundefined{\addperiod}}{
  \newcommand{\addperiod}[1]{#1.} % old memoir needs this, new memoir doesn't.
}{
}
\setsubsubsecheadstyle{\bfseries\raggedright\addperiod}
\setbeforesubsubsecskip{0pt}
\setsubsubsecindent{\parindent}
\setaftersubsubsecskip{-0.5em}

\makechapterstyle{appendix}{%
\renewcommand{\chapnamefont}{\normalfont\bfseries}
\renewcommand{\chapnumfont}{\normalfont\bfseries}
\renewcommand{\chaptitlefont}{\normalfont\bfseries}

% Remove extra space at top of chapter heading
\renewcommand{\chapterheadstart}{}
% Print ``CHAPTER'' in all caps, centered
\renewcommand{\printchaptername}{%
  \beforepartskip \centering \chapnamefont \MakeUppercase{\@chapapp}
}
\renewcommand{\printchapternum}{%
  \centering \chapnumfont \thechapter\space
}
\setlength{\midchapskip}{0pt}
\renewcommand{\printchaptertitle}[1]{%
  \centering \chaptitlefont \MakeUppercase{##1} \beforepartskip \clearpage
}
\setlength{\afterchapskip}{2\baselineskip}
}
% No single lines of a paragraph at the end of a page
\widowpenalty=9999
\clubpenalty=9999

%%%
%%% Bibliography styling
%%%
\RequirePackage{etoolbox}
\AfterEndPreamble{%
  \renewcommand\bibsection{%
  \typeout{bibsection redefinition goes here}%
  \thispagestyle{plain}%
  \ifnobibintoc\else%
  \phantomsection%
  \addtocontents{toc}{\protect\vskip-\baselineskip}%
  \addcontentsline{toc}{chapter}{\texorpdfstring{\uppercase\expandafter{\bibname}}{\bibname}}%
  \fi%
  \beforepartskip%
  \begin{center}%
    bibname goes here \bfseries{\MakeUppercase{\bibname}}%
  \end{center}%
  \afterpartskip%
  \bibmark%
  \clearpage%
  \prebibhook%
}}
%\renewcommand\prebibhook{\begin{Spacing}{1}}
%\renewcommand\postbibhook{\end{Spacing}}

%%%
%%% Vita fixes. A Vita ends up acting like an appendix, but should
%%% never be numbered. Easiest for now to just add a custom line to
%%% the ToC for it.
%%%
\newcommand{\vita}[1]{\InputIfFileExists{#1}{\setupvita}{}}
\newcommand\setupvita{%
  \chapterstyle{default}
  \chapter*{Vita}
  \addtocontents{toc}{\protect\vskip-\baselineskip}
  \addcontentsline{toc}{chapter}%
    {\texorpdfstring{\uppercase\expandafter{Vita}}{Vita}}
}


% For automatic subfigures and subtables
\AfterEndPreamble{
  \newsubfloat{figure}
  \newsubfloat{table}
}

\IfFileExists{\jobname-packages-loaded-before-hyperref.sty}{%
  \RequirePackage{\jobname-packages-loaded-before-hyperref}
}{}

% For all sorts of hyperlinks in the PDF output.
\RequirePackage[pdfcreator={pdfLaTeX with hyperref},
  plainpages=false,
  pdfpagelabels,
	hypertexnames=true,
  colorlinks=false]{hyperref}

\IfFileExists{\jobname-packages-loaded-after-hyperref.sty}{%
  \RequirePackage{\jobname-packages-loaded-after-hyperref}%
}{}

% Fixing the annoyance that references to floats anchor on
% the captions, not the floats themselves.
\RequirePackage[all]{hypcap}
% Make sure sidewaystables and sidewaysfigures are treated properly,
% too.
\hypcapredef{sidewaystable}
\hypcapredef{sidewaysfigure}

% The memoir class has some conflicts with hyperref. The memhfixc
% package fixes these conflicts. Do not remove this line.
\RequirePackage{memhfixc}

%\tracingpatches

\AfterEndPreamble{\let\thepdftitle\thetitle%
  \Substitute*[\renewcommand\thepdftitle]{\thepdftitle}%
  {\\}{ }%
  \hypersetup{pdfauthor=\theauthor,pdftitle=\thepdftitle}%
}

\AfterEndPreamble{
\abstractpage \newpage
\titlepage \newpage
\ifthenelse{\copyrightpagepresent=1}{\copyrightpage}{}
\approvalpage
\@ifundefined{ttu@dedication}{}{\dedicationpage}
\@ifundefined{ttu@acknowledgments}{}{\acknowledgmentspage}
\chapterstyle{default}
}
%%%
%%% All done. Let memoir finish its setup now.
%%%

\checkandfixthelayout

\endinput
