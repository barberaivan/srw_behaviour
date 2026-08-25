# Insert figures into the manuscript .docx, in teal, following the house rules
# in CLAUDE.md: never touch an existing paragraph, no tracked changes, work on a
# copy, and verify programmatically before installing the result.
#
# This is the exact script that put the Results figures (Figures 3-6) into
# srw_behaviour_paper.docx on 2026-08-24. It is kept because the next batch of
# figures will need the same machinery — the fiddly parts are the image XML, the
# EMU sizing and the byte-identical verification, none of which is worth
# rediscovering.
#
# HOW TO REUSE IT
#   1. Point PAPER at the .docx to edit.
#   2. Rewrite the `after(...)` calls at the bottom: each one anchors on a unique
#      snippet of text of an EXISTING paragraph and appends new paragraphs right
#      after it.
#   3. Run it. It writes a new .docx next to the backup; nothing is installed
#      over the original until you copy it yourself.
#
# It refuses to run if the target already carries the paragraphs it would add,
# so re-running it by mistake cannot insert them twice.
#
# What it verifies before writing the output:
#   - word/document.xml and the rels file still parse as XML;
#   - every original paragraph is still present, byte-identical, in the same
#     order (that is what protects the co-authors' text, the [Claude, ...]
#     markers and the existing images).
#
# Render the result with
#   libreoffice --headless --convert-to pdf <file>.docx
# and look at the pages: it is the cheapest way to catch a broken image or table.

import os, re, shutil, subprocess, zipfile
from xml.etree import ElementTree as ET

PAPER = "/home/ivan/Insync/Whales/Behaviour/Paper/srw_behaviour_paper.docx"
REPO = "/home/ivan/dev/srw_behaviour"
WORK = "/tmp/claude-1000/-home-ivan-dev/d48a8568-d1b9-4a7e-971b-96cc551ecedd/scratchpad/docx"
SRC = os.path.join(WORK, "work")
BACKUP = os.path.join(WORK, "srw_behaviour_paper_BACKUP.docx")
OUT = os.path.join(WORK, "srw_behaviour_paper_new.docx")

CONTENT_W = 5731510  # EMU: page width minus margins, same as the existing images
MAX_H = 8400000

GUARD = "Claude, Figure 3, version A"

with zipfile.ZipFile(PAPER) as _z:
    if GUARD in _z.read("word/document.xml").decode("utf-8"):
        raise SystemExit(
            "%s already contains %r — it looks like this script has already run "
            "against it. Point PAPER at the pre-edit backup, or change GUARD."
            % (PAPER, GUARD))

shutil.copy2(PAPER, BACKUP)
if os.path.exists(SRC):
    shutil.rmtree(SRC)
os.makedirs(SRC)
with zipfile.ZipFile(PAPER) as z:
    names = z.namelist()
    z.extractall(SRC)

doc_path = os.path.join(SRC, "word/document.xml")
doc = open(doc_path, encoding="utf-8").read()
original_paras = re.findall(r"<w:p[ >].*?</w:p>", doc, re.S)

rels_path = os.path.join(SRC, "word/_rels/document.xml.rels")
rels = open(rels_path, encoding="utf-8").read()

# ---------------------------------------------------------------- helpers

def esc(t):
    return (t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))

def para(text, italic=False):
    rpr = "<w:rPr>" + ("<w:i/><w:iCs/>" if italic else "") + '<w:color w:val="008080"/></w:rPr>'
    return ('<w:p><w:pPr><w:pStyle w:val="normal1"/></w:pPr>'
            '<w:r>' + rpr + '<w:t xml:space="preserve">' + esc(text) + "</w:t></w:r></w:p>")

def png_size(path):
    with open(path, "rb") as f:
        head = f.read(33)
    w = int.from_bytes(head[16:20], "big")
    h = int.from_bytes(head[20:24], "big")
    return w, h

state = {"rid": 100, "docpr": 100, "img": 100}

def add_image(png_path):
    """copy the png into the package, register the relationship, return a <w:p>"""
    state["img"] += 1
    state["rid"] += 1
    state["docpr"] += 1
    name = "image%d.png" % state["img"]
    rid = "rId%d" % state["rid"]
    shutil.copy2(png_path, os.path.join(SRC, "word/media", name))

    global rels
    rels = rels.replace(
        "</Relationships>",
        '<Relationship Id="%s" Type="http://schemas.openxmlformats.org/officeDocument/'
        '2006/relationships/image" Target="media/%s"/></Relationships>' % (rid, name))

    pw, ph = png_size(png_path)
    cx = CONTENT_W
    cy = int(round(cx * ph / pw))
    if cy > MAX_H:
        cy = MAX_H
        cx = int(round(cy * pw / ph))
    did = state["docpr"]
    return (
        '<w:p><w:pPr><w:pStyle w:val="normal1"/></w:pPr><w:r><w:rPr></w:rPr><w:drawing>'
        '<wp:inline distT="0" distB="0" distL="0" distR="0">'
        '<wp:extent cx="%d" cy="%d"/><wp:effectExtent l="0" t="0" r="0" b="0"/>'
        '<wp:docPr id="%d" name="%s" descr=""></wp:docPr><wp:cNvGraphicFramePr>'
        '<a:graphicFrameLocks xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" noChangeAspect="1"/>'
        '</wp:cNvGraphicFramePr><a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
        '<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
        '<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
        '<pic:nvPicPr><pic:cNvPr id="%d" name="%s" descr=""></pic:cNvPr>'
        '<pic:cNvPicPr><a:picLocks noChangeAspect="1" noChangeArrowheads="1"/></pic:cNvPicPr></pic:nvPicPr>'
        '<pic:blipFill><a:blip r:embed="%s"></a:blip><a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
        '<pic:spPr bwMode="auto"><a:xfrm><a:off x="0" y="0"/><a:ext cx="%d" cy="%d"/></a:xfrm>'
        '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic></a:graphicData></a:graphic>'
        "</wp:inline></w:drawing></w:r></w:p>" % (cx, cy, did, name, did, name, rid, cx, cy))

def find_para(snippet):
    """the one existing paragraph whose text contains `snippet`"""
    hits = [p for p in original_paras
            if snippet in "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", p, re.S))]
    assert len(hits) == 1, (snippet, len(hits))
    return hits[0]

inserts = []   # (anchor paragraph, xml to put right after it)

def after(snippet, blocks):
    inserts.append((find_para(snippet), "".join(blocks)))

# ---------------------------------------------------------------- content

CAPTION_3 = (
    "Disturbed and undisturbed are not two groups of intervals but two reference "
    "points on the estimated disturbance scale z. The undisturbed distribution is "
    "what the model predicts for a whale at z = 0, one that has been free of "
    "attacks long enough to have recovered completely; the disturbed one, for a "
    "whale at z = 1, under sustained attack. For each of the two we show the "
    "steady-state distribution of that year's transition matrix, that is, the "
    "proportion of time a whale would spend in each behaviour category if it "
    "stayed at that level of disturbance indefinitely. Both are estimated, "
    "counterfactual quantities: no whale in our data is ever observed at exactly "
    "z = 0 or z = 1, and no interval is labelled as one or the other. Points are "
    "posterior means and vertical lines 95 % highest-density intervals. The "
    "observed proportions, computed over every interval of the year with the "
    "unobserved behaviours imputed by the model, are shown alongside them and "
    "carry no interval, being data rather than estimates."
)

after("Figure 3. Time series of observed and estimated behaviour", [
    para("Claude, Figure 3, version A — mothers and calves as separate blocks, "
         "merged with patchwork; the y scale is shared row by row between the two "
         "halves. File: plots/figure_3_behaviour_timeseries_faceted_v01.png", italic=True),
    add_image(REPO + "/plots/figure_3_behaviour_timeseries_faceted_v01.png"),
    para("Claude, Figure 3, version B — mothers and calves in the same panels, "
         "keyed by point shape and line type, with the three conditions dodged and "
         "joined by lines. File: plots/figure_3_behaviour_timeseries_dodged_v01.png",
         italic=True),
    add_image(REPO + "/plots/figure_3_behaviour_timeseries_dodged_v01.png"),
    para("Claude, proposed text for the caption, answering the marker above: " + CAPTION_3),
])

after("Produce both versions and we", [
    para("Claude, Figure 4, version A — mothers and calves as two plots merged with "
         "patchwork, tagged A and B, with heights in proportion to their number of "
         "periods (5 and 3). Single-year periods (1995 for mothers, 2011-2013 for "
         "calves, which only contributes 2013) are plain model predictions; "
         "multi-year periods average the years sampled in the period within each "
         "posterior sample, and the posterior is summarised only at the end. "
         "File: plots/figure_4_attack_scenarios_patchwork_v01.png", italic=True),
    add_image(REPO + "/plots/figure_4_attack_scenarios_patchwork_v01.png"),
    para("Claude, Figure 4, version B — a single plot with the rows nested as "
         "mother/calf x period. File: plots/figure_4_attack_scenarios_nested_v01.png",
         italic=True),
    add_image(REPO + "/plots/figure_4_attack_scenarios_nested_v01.png"),
])

after("rehac", [
    para("Claude, Figure 5 — remade with the current predictions and with the "
         "calves added to both short-term effects, potential and observed. Mothers "
         "and calves are keyed by colour. The long-term and total effects are "
         "mothers only, as they are measured against 1995 and the calves' series "
         "starts in 2013. File: plots/figure_5_effects_summary_v01.png", italic=True),
    add_image(REPO + "/plots/figure_5_effects_summary_v01.png"),
])

after("Do we have new mortality data up to 2025", [
    para("Claude, Figure 6, version A — panel B with mothers only, the four curves "
         "asked for (rest and travel, each split by exposure; in-place activity is "
         "left out, so the curves do not add up to one). Panel A shows the "
         "proportion of intervals with an attack on either member of the pair (left "
         "axis) and the mean number of attacks per interval on the mother and on the "
         "calf (right axis). Lines break at the years without sampling. "
         "File: plots/figure_6_attacks_behaviour_mortality_mothers_v01.png", italic=True),
    add_image(REPO + "/plots/figure_6_attacks_behaviour_mortality_mothers_v01.png"),
    para("Claude, Figure 6, version B — the same, with panel B showing mothers and "
         "calves together, keyed by point shape and line type. "
         "File: plots/figure_6_attacks_behaviour_mortality_both_v01.png", italic=True),
    add_image(REPO + "/plots/figure_6_attacks_behaviour_mortality_both_v01.png"),
    para("Meri: panel C is empty on purpose. Which calf mortality series should go "
         "there — the dead calves counted by the Southern Right Whale Health "
         "Monitoring Program, as in Piotto et al. (2024)? And do we have mortality "
         "data up to 2025, or does the series still stop earlier? The panel draws "
         "itself as soon as the numbers are dropped into the mortality data frame in "
         "plots/figure_6_attacks_behaviour_mortality.R."),
])

# ---------------------------------------------------------------- apply

for anchor, block in inserts:
    assert doc.count(anchor) == 1, "anchor not unique"
    doc = doc.replace(anchor, anchor + block)

open(doc_path, "w", encoding="utf-8").write(doc)
open(rels_path, "w", encoding="utf-8").write(rels)

# ---------------------------------------------------------------- verify

ET.fromstring(doc.encode("utf-8"))          # the XML must parse
ET.fromstring(rels.encode("utf-8"))

new_paras = re.findall(r"<w:p[ >].*?</w:p>", doc, re.S)
kept = [p for p in new_paras if p in original_paras]
assert len(original_paras) == len(kept), (len(original_paras), len(kept))
# every original paragraph, byte-identical and in the same order
j = 0
for p in new_paras:
    if j < len(original_paras) and p == original_paras[j]:
        j += 1
assert j == len(original_paras), (j, len(original_paras))
print("original paragraphs preserved:", j, "/", len(original_paras))
print("paragraphs added:", len(new_paras) - len(original_paras))

# repack, keeping the original entry order and adding the new media
with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
    for n in names:
        z.write(os.path.join(SRC, n), n)
    for n in sorted(os.listdir(os.path.join(SRC, "word/media"))):
        arc = "word/media/" + n
        if arc not in names:
            z.write(os.path.join(SRC, arc), arc)
print("written:", OUT)
