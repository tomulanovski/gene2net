"""Assemble the method-chapter drafts into one Word .docx.

Reads the section markdown files in order, strips each file's meta "DRAFT for the thesis..."
paragraph, and renders headings, paragraphs, bold, inline code, bullet lists, pipe tables, and
figures into a docx. The method-architecture section is salvaged from the archived method draft,
truncated before its old "From predictions to a network" section, which the decode section
supersedes.

Tables and figures are captioned, numbered in document order, and cross-referenced:

    Table: {#tab:oracle} Caption text.              directly above a pipe table
    Figure: {#fig:frac} figures/frac.png | Caption.  on its own line
    ... as @tab:oracle shows ...                    in running text, becomes "Table 4"

Table captions go above the table and figure captions below the figure, each with a live SEQ
field so Word can renumber them and build a list of tables or figures. The build stops if a
table has no caption, a caption has no label, a label is used twice or never referenced in the
running text, a reference names no label, or a figure file is missing.

No em-dashes: the drafts contain none, and the script warns if any en-dash or em-dash appears.

Run:  python scripts/build_chapter_docx.py
Out:  docs/gene2net_method_chapter.docx
"""
import os
import re

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.join(HERE, "..", "docs")

# (filename relative to docs/, truncate-before-marker or None). Order is the chapter order.
SECTIONS = [
    ("chapter_introduction_draft.md", None),
    ("archive/chapter_method_draft.md", "## From predictions to a network"),
    ("chapter_decode_draft.md", None),
    ("chapter_setup_draft.md", None),
    ("chapter_results_benchmark_draft.md", None),
    ("chapter_diagnostic_draft.md", None),
    ("chapter_retention_draft.md", None),
    ("chapter_feature_importance_draft.md", None),
    ("chapter_limitations_draft.md", None),
    ("chapter_future_work_draft.md", None),
    ("chapter_conclusion_draft.md", None),
    ("chapter_appendix_tables_draft.md", None),
]

CHAPTER_TITLE = "A learned detect-then-place method for polyploid network reconstruction"
OUT = os.path.join(DOCS, "gene2net_method_chapter.docx")
FIGURE_WIDTH = Inches(6.3)

CAPTION_RE = re.compile(r"^(Table|Figure):\s*\{#((?:tab|fig):[A-Za-z0-9_-]+)\}\s*(.*)$")
REF_RE = re.compile(r"@((?:tab|fig):[A-Za-z0-9_-]+)")


def strip_meta(text):
    """Remove the 'DRAFT for the thesis...' meta paragraph (until the next blank line)."""
    out, skip = [], False
    for line in text.split("\n"):
        if line.strip().startswith("DRAFT for the thesis"):
            skip = True
            continue
        if skip:
            if line.strip() == "":
                skip = False
            continue
        out.append(line)
    return "\n".join(out)


def truncate_before(text, marker):
    if marker and marker in text:
        return text[: text.index(marker)]
    return text


def add_runs(paragraph, text):
    """Add text to a paragraph, honoring **bold** and `code` inline markers."""
    pattern = re.compile(r"(\*\*.+?\*\*|`[^`]+`)")
    pos = 0
    for m in pattern.finditer(text):
        if m.start() > pos:
            paragraph.add_run(text[pos:m.start()])
        tok = m.group(0)
        if tok.startswith("**"):
            paragraph.add_run(tok[2:-2]).bold = True
        else:
            r = paragraph.add_run(tok[1:-1])
            r.font.name = "Consolas"
        pos = m.end()
    if pos < len(text):
        paragraph.add_run(text[pos:])


def clean_heading(s):
    return s.replace("**", "").replace("`", "").strip()


def parse_row(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def is_sep(line):
    return re.match(r"^\s*\|?[\s:|-]+\|?[\s:|-]*$", line) is not None and set(line.strip()) <= set("|-: ")


def is_caption(line):
    return line.strip().startswith(("Table:", "Figure:"))


def collect_labels(texts):
    """Number every captioned table and figure in document order: label -> (kind, number)."""
    labels, count = {}, {"Table": 0, "Figure": 0}
    for text in texts:
        for line in text.split("\n"):
            if not is_caption(line):
                continue
            m = CAPTION_RE.match(line.strip())
            if not m:
                raise SystemExit(f"Caption has no {{#tab:...}} or {{#fig:...}} label: {line.strip()[:80]}")
            kind, label, _ = m.groups()
            if (kind == "Table") != label.startswith("tab:"):
                raise SystemExit(f"Label prefix does not match the caption kind: {line.strip()[:80]}")
            if label in labels:
                raise SystemExit(f"Label used twice: {label}")
            count[kind] += 1
            labels[label] = (kind, count[kind])
    return labels


def check_references(texts, labels):
    """Every reference must name a label, and every label must be referenced in running text."""
    used = set()
    for text in texts:
        for line in text.split("\n"):
            refs = REF_RE.findall(line)
            for ref in refs:
                if ref not in labels:
                    raise SystemExit(f"Reference to an unknown label: @{ref}")
            if not is_caption(line):
                used.update(refs)
    unused = sorted(set(labels) - used)
    if unused:
        raise SystemExit("Never referenced in the running text: " + ", ".join("@" + u for u in unused))


def resolve(text, labels):
    return REF_RE.sub(lambda m: "{} {}".format(*labels[m.group(1)]), text)


def add_caption(doc, kind, number, text, labels, above):
    """Caption paragraph "Kind N. text", N as a live SEQ field so Word can renumber it."""
    p = doc.add_paragraph(style="Caption")
    p.paragraph_format.keep_with_next = above
    p.add_run(f"{kind} ")
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), f" SEQ {kind} \\* ARABIC ")
    run = OxmlElement("w:r")
    txt = OxmlElement("w:t")
    txt.text = str(number)
    run.append(txt)
    fld.append(run)
    p._p.append(fld)
    add_runs(p, ". " + resolve(text, labels))


def add_table(doc, rows):
    ncols = max(len(r) for r in rows)
    table = doc.add_table(rows=len(rows), cols=ncols)
    table.style = "Table Grid"
    for i, row in enumerate(rows):
        for j in range(ncols):
            cell_text = row[j] if j < len(row) else ""
            p = table.cell(i, j).paragraphs[0]
            add_runs(p, cell_text)
            if i == 0:
                for run in p.runs:
                    run.bold = True


def add_figure(doc, rel_path, caption, number, labels):
    path = os.path.normpath(os.path.join(DOCS, rel_path))
    if not os.path.exists(path):
        raise SystemExit(f"Figure file not found: {path}")
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.keep_with_next = True
    p.add_run().add_picture(path, width=FIGURE_WIDTH)
    add_caption(doc, "Figure", number, caption, labels, above=False)


def render_markdown(doc, text, labels):
    lines = text.split("\n")
    i, para, pending = 0, [], None   # pending: (number, caption) of the next table

    def flush():
        nonlocal para
        if para:
            add_runs(doc.add_paragraph(), resolve(" ".join(para).strip(), labels))
            para = []

    while i < len(lines):
        line = lines[i]
        s = line.strip()
        if not s:
            flush(); i += 1; continue
        if s.startswith("### "):
            flush(); doc.add_heading(clean_heading(s[4:]), level=3); i += 1; continue
        if s.startswith("## "):
            flush(); doc.add_heading(clean_heading(s[3:]), level=2); i += 1; continue
        if s.startswith("# "):
            flush(); doc.add_heading(clean_heading(s[2:]), level=1); i += 1; continue
        if is_caption(s):
            flush()
            kind, label, rest = CAPTION_RE.match(s).groups()
            number = labels[label][1]
            i += 1
            if kind == "Figure":
                if " | " not in rest:
                    raise SystemExit(f"Figure line needs 'path | caption': {s[:80]}")
                rel_path, caption = rest.split(" | ", 1)
                add_figure(doc, rel_path.strip(), caption.strip(), number, labels)
                continue
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i >= len(lines) or not lines[i].strip().startswith("|"):
                raise SystemExit(f"Table caption is not followed by a table: {s[:80]}")
            pending = (number, rest.strip())
            continue
        if s.startswith("|"):
            flush()
            if pending is None:
                raise SystemExit(f"Table has no caption line above it: {s[:80]}")
            block = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                block.append(lines[i]); i += 1
            rows = [parse_row(b) for b in block if not is_sep(b)]
            if rows:
                add_caption(doc, "Table", pending[0], pending[1], labels, above=True)
                add_table(doc, rows)
            pending = None
            continue
        if re.match(r"^\s*-\s+", line):
            flush()
            while i < len(lines) and re.match(r"^\s*-\s+", lines[i]):
                item = re.sub(r"^\s*-\s+", "", lines[i]).strip()
                add_runs(doc.add_paragraph(style="List Bullet"), resolve(item, labels))
                i += 1
            continue
        para.append(s); i += 1
    flush()


def main():
    texts = []
    for fname, marker in SECTIONS:
        path = os.path.join(DOCS, fname)
        if not os.path.exists(path):
            raise SystemExit(f"Missing section file: {path}")
        text = open(path, encoding="utf-8").read()
        texts.append((fname, strip_meta(truncate_before(text, marker)).strip()))

    labels = collect_labels([t for _, t in texts])
    check_references([t for _, t in texts], labels)

    doc = Document()
    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(12)
    doc.add_heading(CHAPTER_TITLE, level=0)

    warned = False
    for fname, text in texts:
        if ("—" in text or "–" in text) and not warned:
            print(f"WARNING: dash character found in {fname} (em/en dash)")
            warned = True
        render_markdown(doc, text, labels)

    doc.save(OUT)
    kinds = [k for k, _ in labels.values()]
    print(f"Wrote {OUT}")
    print(f"Sections: {len(SECTIONS)}, tables: {kinds.count('Table')}, figures: {kinds.count('Figure')}")


if __name__ == "__main__":
    main()
