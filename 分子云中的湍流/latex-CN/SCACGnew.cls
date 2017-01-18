\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{SCAG}[2008/11/15 v1, by Liu Yun]

\let\CCTCJKfonts=1
\LoadClass[twoside,CJK]{cctart}
%% math packages
\RequirePackage{amsmath,amsthm,cct,amssymb,txfonts,times,bm,tabularx,geometry,fancyhdr,
color,multicol,array,pifont,amscd,latexsym,upgreek,cite}

\usepackage{mathptmx}
%% graphics packages
\RequirePackage{graphicx}
%% other packages
\RequirePackage{url,float}
%% check pdfTeX mode
\RequirePackage{ifpdf}


%\input psfig.sty

\columnsep 0.7cm  %两栏之间的距离


\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{SCGC}[2008/11/15 v1]



%------------------------ Page layout ----------------------------------

\textwidth=170truemm
\textheight=210truemm
\topmargin=-10pt
\headsep=7truemm
\oddsidemargin=-17pt
\evensidemargin=-17pt
\parindent=0.74pt
\footskip=9mm
\renewcommand{\baselinestretch}{1.31}
\renewcommand{\arraystretch}{1.21}
\newcommand{\bahao}{\fontsize{9pt}{\baselineskip}\selectfont}
\newcommand{\xiaobahao}{\fontsize{8pt}{\baselineskip}\selectfont}
\newcommand{\wuhao}{\fontsize{9.5pt}{0.76\baselineskip}\selectfont}
%\newcommand{\xiaosanhao}{\fontsize{15pt}{2\baselineskip}\selectfont}
\newcommand{\zaiyao}{\fontsize{9.2pt}{\baselineskip}\selectfont}
\newcommand{\xiaowuhao}{\fontsize{9pt}{0.76\baselineskip}\selectfont}
\newcommand{\xiaosihao}{\fontsize{12pt}{\baselineskip}\selectfont}
\newcommand{\liuhao}{\fontsize{8pt}{\baselineskip}\selectfont}
\newcommand{\xiaoliuhao}{\fontsize{7.5pt}{\baselineskip}\selectfont}

%----------------------- Title Page ----------------------------------

\let\@Year\@empty \def\Year#1{\def\@Year{#1}}
\let\@Vol\@empty \def\Vol#1{\def\@Vol{#1}}
\let\@No\@empty \def\No#1{\def\@No{#1}}
\let\@BeginPage\@empty \def\BeginPage#1{\def\@BeginPage{#1}}
\let\@DOI\@empty \def\DOI#1{\def\@DOI{#1}}
\let\@EndPage\@empty \def\EndPage#1{\def\@EndPage{#1}}
\let\@ReceivedDate\@empty \def\ReceivedDate#1{\def\@ReceivedDate{#1}}
\let\@AcceptedDate\@empty \def\AcceptedDate#1{\def\@AcceptedDate{#1}}
\let\@Email\@empty \def\Email#1{\def\@Email{#1}}
\let\@Foundation\@empty \def\Foundation#1{\def\@Foundation{#1}}
\def\REF#1{\par\hangindent\parindent\indent\llap{#1\enspace}\ignorespaces}
\let\@Authormark\@empty \def\Authormark#1{\def\@Authormark{#1}}
\let\@AuthormarkCite\@empty \def\AuthormarkCite#1{\def\@AuthormarkCite{#1}}
\let\@Correspondingauthor\@empty \def\Correspondingauthor#1{\def\@Correspondingauthor{#1}}
\let\@ArtNo\@empty \def\ArtNo#1{\def\@ArtNo{#1}}

%\let\@Thanks\@empty \def\Thanks#1{\def\@Thanks{#1}}
\let\@title\@empty
\renewcommand\title[2][]{
\vspace*{-2.04cm}\noindent
  \parbox[c]{12.8cm}{
    {\small{\bahao}\sf\heiti
中国科学\,: 物理学~~力学~~天文学\hfill\@Year\,年\hspace{3mm}第\,\@Vol\,卷\hspace{3mm}
第\,\@No\,期\,:\;\@ArtNo} \vspace{1.5truemm} \hrule depth0pt
height0.25truemm width 128 truemm \vspace{1.5truemm} {\sf\xiaowuhao
SCIENTIA SINICA~~~Physica, Mechanica \& Astronomica\hfill phys.scichina.com}}\parbox[c]{4.2cm}
{\hfill\includegraphics{SB-ZG.eps}}
  \vspace{12truemm}
  \def\@title{{\bfseries\huge\heiti\zihao{2}#2}}
  \def\@titlehead{#1}
  \ifx\@titlehead\@empty\def\@titlehead{#2}\fi
}
\let\@author\@empty
\renewcommand*{\author}[3][]{%
  \ifx \@author\@empty
    \def\@author{\kaishu\zihao{4}#2\@authornote{#1}{#3}}
    \def\@authorhead{#2}
  \else
    \expandafter\def\expandafter\@author\expandafter
      {\@author,\,\, #2\@authornote{#1}{#3}}
    \expandafter\def\expandafter\@authorhead\expandafter
      {\@authorhead}
  \fi
}
\def\@authornote#1#2{%
  \def\SCAC@author@address{#1}%
  \def\SCAC@author@footnote{#2}%
  \ifx \SCAC@author@address\@empty
    \ifx \SCAC@author@footnote\@empty
    \else
       \footnote{#2}%
    \fi
  \else
    \textsuperscript{\footnotesize #1}%
    \ifx \SCAC@author@footnote\@empty
    \else
      \textsuperscript{}%
        \footnote{#2}%
    \fi
  \fi}
 \def\@address{}
\newcommand*{\address}[2][]{%
  \ifx \@address\@empty
    \def\@address{\mbox{#1~~#2}}
  \else
    \expandafter\def\expandafter\@address\expandafter
      {\@address\\ \mbox{#1~~#2}}
  \fi}

\renewcommand\maketitle{\par
  \begingroup
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\hbox{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\noindent\@textsuperscript{\normalfont
    \@thefnmark}##1}%
    \@maketitle
\@thanks
    \long\def\@makefntext##1{\noindent##1}%
  \endgroup
  \setcounter{footnote}{0}%
  \setcounter{page}{\@BeginPage}
}
\def\@maketitle{%
  \begingroup\noindent%
    %\let\footnote\thanks
    \begin{flushleft}{\@title}%
    \begin{flushleft}\rmfamily\Large\fangsong\zihao{4}
    \vskip 0.35 true cm%
      \begin{tabular}[t]{@{}l@{}}%
        \@author
      \end{tabular}\vspace{-0.2true cm}%
    \rmfamily\footnotesize\songti\zihao{6}
    \begin{flushleft}
    \@address\\
    \end{flushleft}%
%    \vskip-0.3cm{\zihao{6} *E-mail:\ \@Email}
%    \vskip 0.4 true cm%
    {\zihao{6}收稿日期: \@ReceivedDate;\ \ 接受日期:
    \@AcceptedDate \ \ %*\,通信作者
    \\\vspace{-1.5mm}
\@Foundation}
  \endgroup
}

%----------------------- Running heads ---------------------------------
\def\ps@myheadings{
  \def\@evenfoot{\vbox{\hbox to\textwidth{\footnotesize\rm\hbox
{\hfill} \hfill\songti\rm\zihao{-5}\small\@ArtNo-\thepage \hfill}}}
\def\@oddfoot{\vbox{\hbox to \textwidth{\hfill\hbox{
\songti\rm\zihao{-5}\small\@ArtNo-\thepage}\hfill}}}

\def\@oddhead{\vbox{\hbox to\textwidth{{\hfill}
    \rmfamily\small\songti\zihao{-5} 
    \@Authormark.~~~中国科学\,: ~物理学~~力学~~天文学\hspace{3mm}\@Year\,年\hspace{3mm}第~\@Vol~卷\hspace{3mm}第~\@No~期\hfill }
    \vspace{1.5truemm}
    \hrule depth 0pt height 0.15truemm width\textwidth}
  }
 
 
  \def\@evenhead{\vbox{\hbox to \textwidth{{\hfill}
    \rmfamily\small\songti\zihao{-5}\@Authormark.~~~中国科学\,: ~物理学~~力学~~天文学\hspace{3mm}\@Year\,年\hspace{3mm}第~\@Vol~卷\hspace{3mm}第~\@No~期 \hfill}\vspace{1.5truemm}
    \hrule depth 0pt height 0.15truemm width\textwidth}
  }
  \thispagestyle{fancy}\fancyfoot{}
}
\pagestyle{myheadings}


%-------------------------- Abstract ----------------------------------
\renewenvironment{abstract}
  {\trivlist\itemindent\parindent\item\relax%
    \vskip-4mm
    \textcolor[rgb]{0.4,0.4,0.4}{\rule{17cm}{0.4pt}}
    \vskip2.5mm\noindent\heiti{摘要}\quad\ignorespaces\fangsong }
  {\endtrivlist}


\newcommand*\keywords[1]{
  {\vskip1.8mm\noindent {{\heiti 关键词}\quad\quad\fangsong #1}}
}
\newcommand*\PACS[1]{
  {\vskip3.5mm\noindent{\bf PACS:}~#1}
}

\newcommand*\AutoDOI{
  {\vskip3.5mm\noindent{\bf doi:}~\@DOI}
}

\newcommand*\Thanks[1]{
  {\vspace{6mm}\noindent\parbox[b]{1cm}{\heiti{致谢}}\vspace{-5.7mm}
\parbox[t]{15cm}{\fangsong{#1}\zihao{-5}}
 }
}

%----------------------------------图表相关-----------------------------
\usepackage{booktabs}%%(11月注释：表格宏包，可调整表格中的线宽)

\usepackage{ccaption} %%图表中英双标题
\captionnamefont{\bfseries\footnotesize}%%(11月注释：小号黑体标题编号)
\captiontitlefont{\footnotesize}%%(11月注释：小号标题)
\captiondelim{\quad}%%(11月注释：定义标题编号与标题间空一个汉字)

\renewcommand{\figurename}{\bfseries{图}}





%\usepackage{caption2}
%\renewcommand{\captionfont}{\footnotesize}%%(11月注释：小号标题)
%\renewcommand{\captionlabelfont}{\footnotesize\bfseries}%%(11月注释：小号黑体标题编号)
%\renewcommand{\captionlabeldelim}{\quad}%%(11月注释：定义标题编号与标题间空三格)




  %%(11月注释：适当放宽图形环境中的各种比例参数)
\renewcommand{\textfraction}{0.15}
\renewcommand{\topfraction}{0.85}
\renewcommand{\bottomfraction}{0.65}
\renewcommand{\floatpagefraction}{0.60}


%-------------------------- Sections ----------------------------------
%-------------------------- Sections ----------------------------------
\renewcommand\sectionname{\thesection}
\renewcommand\sectionformat{}
\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {-4.5ex \@plus -1ex \@minus -.2ex}%
                                   {2.3ex \@plus.2ex}%
                                   {\normalfont\large\bfseries\zihao{-4}}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\normalsize\bfseries\heiti}}

\def\sec#1{\vspace{7mm}\leftline{\bf\bfseries\large\zihao{-4}#1}\vspace{3mm}}%定义三级章节
\def\subsec#1{\vspace{4mm}\leftline{\bfseries\normalsize\bf{#1}\vspace{2mm}}} %定义二级章节
\def\subsubsec#1{\vspace{4mm}\leftline{\bfseries\normalsize\bf{#1}\vspace{1mm}}} %定义二级章节
\def\th#1{{\bf\heiti\zihao{5} #1}\quad} %定义,定理,引理,证明等

\definecolor{grey}{rgb}{0.5,0.5,0.5}%(11月注释：定义灰色)

\renewcommand\refname{\parbox[c]{2cm}{\heiti\zihao{5}参考文献}
\hspace{-6mm}\color{grey}\rule[0.01\baselineskip]{0.9\textwidth}{0.5pt}}
\def\@biblabel#1{#1~}
\renewenvironment{thebibliography}[1]
     {\section*{\refname}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \songti\leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \itemsep 0ex
      \parskip 0ex
      \footnotesize
      \sloppy}

%----------------------- Math ---------------------------------
\def\leq{\leqslant}
\def\geq{\geqslant}
\def\le{\leqslant}
\def\ge{\geqslant}
\abovedisplayskip=3pt plus 1pt minus 1pt
\belowdisplayskip=3pt plus 1pt minus 1pt
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}%
\newcommand{\Set}[2]{\ensuremath{\left\{ #1 \left| \right. #2 \right\}}}%
\newcommand{\abs}[1]{\ensuremath{\left| #1 \right|}}%
\newcommand{\map}[3]{\ensuremath{#1 : #2 \rightarrow #3}}%
\newcommand{\C}[1]{\ensuremath{\mathcal{#1}}}%
\newcommand{\restrict}[2]{\ensuremath{#1 \left|_{#2}\right.}}%
\newcommand{\NN}{\mathbb{N}^0}
\newcommand{\PN}{\mathbb{N}}
\newcommand{\syn}{\mathbf{S}}
\newcommand{\Syn}{\mathbf{M}}
\allowdisplaybreaks


%-----------------------cite ---------------------------------
%上标引用方法1. 定义"upcite"
\newcommand{\upcite}[1]{\textsuperscript{\cite{#1}}}

%标引用方法2. 重新定义"cite"，（重新定义cite后，不能使用upcite，否则两次上标
%因此，要使用upsite，必须将以下三行命令注释掉。）
%\makeatletter
%    \def\@cite#1#2{\textsuperscript{[{#1\if@tempswa , #2\fi}]}}
%\makeatother

%----------------------- Others ---------------------------------
\newenvironment{listnum}{%
\begin{list}{{\rm(\arabic{enumi})}}{%
\usecounter{enumi}%
\labelwidth=4em%
\labelsep=0.5em%
\leftmargin=0em%
\rightmargin=0em%
\parsep=\parskip%
\itemsep=0ex
\topsep=0em%
\itemindent=3.2em%
}%
}{\end{list}}%

\newenvironment{listroman}{%
\begin{list}{{\rm(\roman{enumi})}}{%
\usecounter{enumi}%
\labelwidth=4em%
\labelsep=0.5em%
\leftmargin=3em%
\rightmargin=2em%
\parsep=\parskip%
\itemsep=0ex
\topsep=0em%
}%
}{\end{list}}%

%----------------------------------------------------
\makeatletter
\newenvironment{tablehere}
  {\def\@captype{table}}
 {}

\newenvironment{figurehere}
 {\def\@captype{figure}}
 {}
\makeatother

\endinput
