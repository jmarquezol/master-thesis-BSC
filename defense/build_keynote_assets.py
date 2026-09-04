#!/usr/bin/env python3
"""Regenerate the Keynote scaffolding from main.tex / main.pdf.

Everything this writes is gitignored and disposable; run it after any change
that moves the page count, then run the AppleScripts.

    cd defense && python3 build_keynote_assets.py

Produces
  main_nofig.tex / .pdf        main.tex with every \\includegraphics blanked
                               (space reserved) so the real figures can be
                               dropped on top as movable Keynote objects
  keynote_pages/slide_NN.pdf   one vector page per PDF page
  keynote_pages/notes_NN.txt   the \\note{} text, on the FIRST page of a frame
                               only (build steps get an empty file)
  keynote_layers/…             the same pair, from the figure-blanked build
  figures_for_keynote/         every figure the deck includes, copied flat
"""
import os, re, shutil, subprocess, sys
import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
VENV = os.path.expanduser("~/.venvs/thesis/bin")
ENV = dict(os.environ, PATH=VENV + ":" + os.environ["PATH"])

NOFIG_PATCH = ("\\let\\KNincludegraphics\\includegraphics\n"
               "\\renewcommand{\\includegraphics}[2][]"
               "{\\phantom{\\KNincludegraphics[#1]{#2}}}\n")


def compile_tex(name):
    subprocess.run(["tectonic", "-X", "compile", name, "--outfmt", "pdf",
                    "-Z", "shell-escape", "--keep-intermediates"],
                   env=ENV, check=True, capture_output=True)


def notes_from_tex(path="main.tex"):
    """The \\note{} bodies, in frame order, as plain text."""
    src = open(path).read()
    # drop comment lines first: the preamble documents \note{...} in a comment,
    # and counting that as a real note shifts every page's script by one.
    src = "\n".join(re.sub(r"(?<!\\)%.*$", "", line) for line in src.splitlines())
    out, i = [], 0
    while True:
        j = src.find("\\note{", i)
        if j < 0:
            break
        k, depth = j + 6, 1
        while depth:
            if src[k] == "{":
                depth += 1
            elif src[k] == "}":
                depth -= 1
            k += 1
        body = src[j + 6:k - 1]
        body = body.replace("---", "\u2014").replace("--", "\u2013")
        body = re.sub(r"\\[a-zA-Z]+\s*", "", body)
        body = re.sub(r"[{}$]", "", body)
        out.append(re.sub(r"\s+", " ", body).strip())
        i = k
    return out


def first_pages(pdf):
    """Page indices (0-based) that start a new frame, read off the footer."""
    doc = fitz.open(pdf)
    starts, prev = [], None
    for n, page in enumerate(doc):
        text = page.get_text()
        # anchor on the footer, not on any "n / m" the slide itself may contain:
        # a figure carrying EIG^{1/2} used to match and merge two frames into one.
        m = re.search(r"September 2026\s+(\d+)\s*/\s*\d+", text)
        key = m.group(1) if m else ("backup:" + text[:40])
        if key != prev:
            starts.append(n)
        prev = key
    doc.close()
    return starts


def split(pdf, outdir, suffix, notes):
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    doc = fitz.open(pdf)
    starts = first_pages(pdf)
    frame_of = {p: i for i, p in enumerate(starts)}
    for n in range(doc.page_count):
        one = fitz.open()
        one.insert_pdf(doc, from_page=n, to_page=n)
        one.save(os.path.join(outdir, f"slide_{n+1:02d}{suffix}.pdf"))
        one.close()
        idx = frame_of.get(n)
        text = notes[idx] if idx is not None and idx < len(notes) else ""
        with open(os.path.join(outdir, f"notes_{n+1:02d}.txt"), "w") as fh:
            fh.write(text)
    doc.close()
    return doc.page_count if False else len(starts)


def copy_figures(outdir="figures_for_keynote"):
    src = open("main.tex").read()
    names = sorted(set(re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", src)
                       + re.findall(r"\\includegraphics\{([^}]+)\}", src)))
    # drop the macro-argument placeholders (\domecell defines \includegraphics{#2})
    names = [n for n in names if not n.startswith("#")]
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    missing = []
    for name in names:
        for root in ("imgs", "../thesis/imgs"):
            cand = os.path.join(root, name)
            if os.path.exists(cand):
                shutil.copy2(cand, os.path.join(outdir, os.path.basename(name)))
                break
        else:
            missing.append(name)
    return len(names) - len(missing), missing


def main():
    print("compiling main.tex …")
    compile_tex("main.tex")

    print("writing main_nofig.tex …")
    src = open("main.tex").read()
    marker = "\\begin{document}\n"
    assert src.count(marker) == 1
    open("main_nofig.tex", "w").write(src.replace(marker, marker + NOFIG_PATCH, 1))
    compile_tex("main_nofig.tex")

    notes = notes_from_tex()
    n1 = split("main.pdf", "keynote_pages", "", notes)
    n2 = split("main_nofig.pdf", "keynote_layers", "_nofig", notes)
    n3, missing = copy_figures()

    pages = fitz.open("main.pdf").page_count
    print(f"pages {pages}  frames {n1} (layers {n2})  figures {n3}")
    if missing:
        print("MISSING figures:", missing, file=sys.stderr)


if __name__ == "__main__":
    main()
