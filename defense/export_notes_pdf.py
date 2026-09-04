#!/usr/bin/env python3
"""Render defense/slide_notes.ipynb to a readable PDF on the Desktop.

No pandoc or browser needed: the notebook's markdown is converted to LaTeX here
and compiled with tectonic, the same engine the deck uses.

    python3 defense/export_notes_pdf.py        (from the repo root)
"""
import json, os, re, shutil, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(ROOT, 'defense', 'slide_notes.ipynb')
BUILD = os.path.join(ROOT, 'defense', '.notes_build')
OUT  = os.path.join(BUILD, 'slide_notes.tex')
DEST = os.path.expanduser('~/Desktop/Defence talk notes.pdf')

UNI = {
 '\u2014':'---', '\u2013':'--', '\u2026':r'\ldots{}',
 '\u21d2':r'$\Rightarrow$', '\u2192':r'$\rightarrow$', '\u21b3':r'$\hookrightarrow$',
 '\u00b7':r'$\cdot$', '\u2713':r'\good{$\checkmark$}', '\u2605':r'$\bigstar$',
 '\u25b8':r'\bul', '\u26a0':r'\warn', '\u2753':r'\qmark',
 '\u2248':r'$\approx$', '\u00d7':r'$\times$', '\u2264':r'$\le$', '\u2265':r'$\ge$',
 '\u00b1':r'$\pm$', '\u221e':r'$\infty$',
}
SPECIAL = {'&':r'\&','%':r'\%','#':r'\#','_':r'\_','{':r'\{','}':r'\}',
           '~':r'\textasciitilde{}','^':r'\textasciicircum{}','$':r'\$',
           '\\':r'\textbackslash{}'}

def esc(t):
    """Escape LaTeX specials in plain text (math already split out)."""
    out=[]
    for ch in t:
        if ch in SPECIAL: out.append(SPECIAL[ch])
        elif ch in UNI:   out.append(UNI[ch])
        else:             out.append(ch)
    return ''.join(out)

def inline(t):
    """One line of markdown -> LaTeX. Styling is tokenised BEFORE the math split,
    so a bold or italic span that contains inline maths survives intact."""
    # 1. tokenise styling on the whole line
    t = re.sub(r'`([^`]+)`',            lambda m: '\x02'+m.group(1)+'\x03', t)
    t = re.sub(r'\*\*([^*]+?)\*\*',      lambda m: '\x04'+m.group(1)+'\x05', t)
    t = re.sub(r'(?<![*\w])\*([^*]+?)\*(?!\w)', lambda m: '\x06'+m.group(1)+'\x07', t)
    # 2. split on maths; escape only the prose
    parts = re.split(r'(\$[^$]*\$)', t)
    res=[]
    for k,seg in enumerate(parts):
        if k % 2:
            for u,v in (('\u2014','---'),('\u2192',r'\rightarrow ')):
                seg = seg.replace(u,v)
            res.append(seg); continue
        seg = re.sub(r'\\([#_*`{}\[\]()+\-.!|\\])', r'\1', seg)
        res.append(esc(seg))
    out = ''.join(res)
    # 3. put the styling back as LaTeX
    for a,b,cmd in (('\x02','\x03',r'\texttt{'),
                    ('\x04','\x05',r'\textbf{'),
                    ('\x06','\x07',r'\emph{')):
        out = out.replace(a, cmd).replace(b, '}')
    return out

def table(rows):
    body = [r for r in rows if not re.match(r'^\s*\|[\s:|-]+\|\s*$', r)]
    cells = [[c.strip() for c in r.strip().strip('|').split('|')] for r in body]
    n = max(len(c) for c in cells)
    cells = [c + ['']*(n-len(c)) for c in cells]
    out = [r'\begin{center}\small', r'\begin{tabular}{@{}' + 'l'*n + r'@{}}', r'\toprule']
    out.append(' & '.join(inline(c) for c in cells[0]) + r' \\')
    out.append(r'\midrule')
    for row in cells[1:]:
        out.append(' & '.join(inline(c) for c in row) + r' \\')
    out += [r'\bottomrule', r'\end{tabular}', r'\end{center}']
    return out

def convert(md):
    lines = md.split('\n'); out=[]; i=0; tbl=[]
    def flush():
        nonlocal tbl
        if tbl: out.extend(table(tbl)); tbl=[]
    while i < len(lines):
        ln = lines[i]; s = ln.strip()
        if s.startswith('|'):
            tbl.append(ln); i+=1; continue
        flush()
        if not s:
            out.append(''); i+=1; continue
        if s.startswith('$$'):                       # display math
            if s.count('$$') >= 2:                   # opens and closes on one line
                body = s[2:s.rindex('$$')]; i+=1
            else:                                    # spans several lines
                blk=[s]; i+=1
                while i < len(lines) and '$$' not in lines[i]: blk.append(lines[i]); i+=1
                if i < len(lines): blk.append(lines[i]); i+=1
                body = '\n'.join(blk).replace('$$','')
            out += [r'\begin{equation*}', body.strip(), r'\end{equation*}']; continue
        if re.match(r'^-{3,}$', s):
            out += [r'\vspace{0.4em}\hrule\vspace{0.6em}']; i+=1; continue
        m = re.match(r'^(#{1,4})\s+(.*)$', s)
        if m:
            lvl, txt = len(m.group(1)), inline(m.group(2))
            cmd = {1:r'\notesection', 2:r'\notesub', 3:r'\notesubsub', 4:r'\notesubsub'}[lvl]
            out += ['', cmd+'{'+txt+'}']; i+=1; continue
        if s.startswith('>'):                        # blockquote -> say-this box
            blk=[]
            while i < len(lines) and lines[i].strip().startswith('>'):
                blk.append(lines[i].strip()[1:].strip()); i+=1
            out += [r'\begin{sayblock}', inline(' '.join(blk)), r'\end{sayblock}']; continue
        indent = len(ln) - len(ln.lstrip())
        m = re.match(r'^(\d+)\.\s+(.*)$', s)
        if m:
            body=[m.group(2)]; i+=1
            # NB: $$ must stop a numbered item, or a display equation inside a list is
            # swallowed as continuation prose and emitted as raw LaTeX source.
            while i < len(lines) and lines[i].strip() and not re.match(r'^\s*(\d+\.|[▸↳|#>]|-{3,}|\$\$)', lines[i]):
                body.append(lines[i].strip()); i+=1
            out.append(r'\numitem{'+m.group(1)+'}{'+inline(' '.join(body))+'}'); continue
        # ordinary paragraph / bullet: gather continuation lines
        body=[s]; i+=1
        while i < len(lines) and lines[i].strip() and \
              not re.match(r'^\s*([▸↳|#>]|\d+\.|-{3,}|\$\$)', lines[i]):
            body.append(lines[i].strip()); i+=1
        txt = inline(' '.join(body))
        if indent >= 2 and s.startswith('\u21b3'):
            out.append(r'\subpt{'+txt+'}')
        else:
            out.append(r'\item[]\relax '+txt if False else r'\noindent '+txt+r'\par\smallskip')
        continue
    flush()
    return '\n'.join(out)

PRE = r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=2.1cm,top=2.0cm,bottom=2.0cm]{geometry}
\usepackage[T1]{fontenc}\usepackage{lmodern}
\usepackage{amsmath,amssymb,bm}
\usepackage{booktabs,array,longtable}
\usepackage{xcolor}
\usepackage{enumitem}
\usepackage{needspace}
\usepackage[colorlinks=true,linkcolor=accent,urlcolor=accent]{hyperref}
\definecolor{ink}{RGB}{31,42,54}
\definecolor{accent}{RGB}{5,141,199}
\definecolor{warnc}{RGB}{176,58,46}
\definecolor{goodc}{RGB}{17,122,101}
\color{ink}
% the deck's macros, so the maths in the notes compiles
\newcommand{\E}{\mathcal{E}}
\newcommand{\Lech}{\mathcal{L}}
\newcommand{\ket}[1]{\left|#1\right\rangle}
\newcommand{\bra}[1]{\left\langle#1\right|}
\newcommand{\braket}[2]{\left\langle#1\middle|#2\right\rangle}
% markers
\newcommand{\bul}{\textcolor{accent}{$\blacktriangleright$}\ }
\newcommand{\warn}{\textcolor{warnc}{\textbf{!}}\ }
\newcommand{\qmark}{\textcolor{accent}{\textbf{Q}}\ }
\newcommand{\good}[1]{\textcolor{goodc}{#1}}
\newcommand{\subpt}[1]{\par\smallskip\hangindent=1.6em\hspace{1.6em}\small #1\par\smallskip}
\newcommand{\numitem}[2]{\par\smallskip\hangindent=1.6em\noindent\textbf{#1.}\ #2\par\smallskip}
\newenvironment{sayblock}
  {\par\smallskip\noindent\begin{list}{}{\leftmargin=1.2em\rightmargin=0.6em}\item[]
   \itshape\color{ink!88}}
  {\end{list}\par\smallskip}
\newcommand{\notesection}[1]{\needspace{4\baselineskip}\par\bigskip
  {\Large\bfseries\color{accent!80!black}#1}\par\smallskip\hrule height 1pt\par\medskip}
\newcommand{\notesub}[1]{\needspace{4\baselineskip}\par\medskip
  {\large\bfseries\color{ink}#1}\par\smallskip}
\newcommand{\notesubsub}[1]{\needspace{3\baselineskip}\par\smallskip
  {\bfseries\color{ink!85}#1}\par\smallskip}
\setlength{\parindent}{0pt}
\begin{document}
"""

nb = json.load(open(SRC))
body = '\n\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type']=='markdown')
os.makedirs(BUILD, exist_ok=True)
open(OUT,'w').write(PRE + convert(body) + '\n\\end{document}\n')

env = dict(os.environ, PATH=os.path.expanduser('~/.venvs/thesis/bin') + ':' + os.environ['PATH'])
subprocess.run(['tectonic','-X','compile','slide_notes.tex','--outfmt','pdf'],
               cwd=BUILD, env=env, check=True, capture_output=True)
shutil.copy(os.path.join(BUILD,'slide_notes.pdf'), DEST)
print('wrote', DEST)
