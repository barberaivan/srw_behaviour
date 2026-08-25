# Insert the Results text and three figure captions into the manuscript .docx,
# in teal, following the house rules in CLAUDE.md: never touch an existing
# paragraph, no tracked changes, work on a copy, and verify programmatically
# before installing the result.
#
# Written on 2026-08-25, answering these markers in the Results section:
#
#   [Claude, write the results sections. A paragraph or 2 per image ...]
#   [Claude, complete this image caption]            (Figure 4)
#   [Claude, explain again here the meaning of each effect ...]  (Figure 5)
#   [Claude, make this caption a bit more complete]  (Figure 6)
#
# The Figure 2 and Figure 3 captions say "do not edit", so they are left alone.
#
# Every number quoted in the text below was recomputed from the committed
# prediction tables, not read off the figures. See the ROADMAP entry for which
# file each family of numbers comes from.
#
# TWO CONVENTIONS WORTH KNOWING BEFORE EDITING THIS FILE
#   - The house rule forbids modifying a paragraph a co-author wrote, and the
#     captions are theirs. So a completed caption goes in as a NEW teal
#     paragraph right below the original, ready for Ivan to swap in by hand.
#   - Results text for a figure must sit BEFORE the figure, so each block is
#     anchored on the paragraph that precedes the figure's images.
#
# It refuses to run if the target already carries the text it would add, so
# re-running it by mistake cannot insert it twice.
#
# Render the result with
#   libreoffice --headless --convert-to pdf <file>.docx
# and look at the pages.

import os, re, shutil, zipfile
from xml.etree import ElementTree as ET

PAPER = "/home/ivan/Insync/Whales/Behaviour/Paper/srw_behaviour_paper.docx"
WORK = os.environ.get("CLAUDE_JOB_DIR", "/tmp") + "/tmp/docx_results"
SRC = os.path.join(WORK, "work")
BACKUP = os.path.join(
    os.path.dirname(PAPER), "srw_behaviour_paper_BACKUP_2026-08-25_pre-results.docx")
OUT = os.path.join(WORK, "srw_behaviour_paper_new.docx")

GUARD = "Behaviour composition and its change through the years"

os.makedirs(WORK, exist_ok=True)

with zipfile.ZipFile(PAPER) as _z:
    if GUARD in _z.read("word/document.xml").decode("utf-8"):
        raise SystemExit(
            "%s already contains %r: it looks like this script has already run "
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


# ---------------------------------------------------------------- helpers

def esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def rpr(italic=False, bold=False, underline=False):
    return ("<w:rPr>"
            + ("<w:b/><w:bCs/>" if bold else "")
            + ("<w:i/><w:iCs/>" if italic else "")
            + ('<w:u w:val="single"/>' if underline else "")
            + '<w:color w:val="008080"/></w:rPr>')


def para(text, italic=False, bold=False, underline=False):
    """one teal body paragraph"""
    return ('<w:p><w:pPr><w:pStyle w:val="normal1"/></w:pPr><w:r>'
            + rpr(italic, bold, underline)
            + '<w:t xml:space="preserve">' + esc(text) + "</w:t></w:r></w:p>")


def heading(text):
    """subsection heading: underlined, as in the Methods"""
    return para(text, underline=True)


def blank():
    return '<w:p><w:pPr><w:pStyle w:val="normal1"/></w:pPr></w:p>'


def find_para(snippet):
    """the one existing paragraph whose text contains `snippet`"""
    hits = [p for p in original_paras
            if snippet in "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", p, re.S))]
    assert len(hits) == 1, (snippet, len(hits))
    return hits[0]


inserts = []   # (anchor paragraph, xml to put right after it)


def after(snippet, blocks):
    inserts.append((find_para(snippet), "".join(blocks)))


# ---------------------------------------------------------------- the text

# ---- Figure 2: the etogram ------------------------------------------------

R2 = [
    heading("Behaviour composition and its change through the years"),
    blank(),
    para(
        "Pooled over the whole series, mother-calf pairs spent most of their "
        "time travelling slowly under water: this single category accounts for "
        "31.7 % of the mothers' intervals and 31.8 % of the calves'. Rest came "
        "next, taking 33.5 % of the mothers' time and 26.3 % of the calves', "
        "and the rest of the budget was split between fast travel (13.3 % and "
        "14.2 %) and in-place activity (8.3 % and 18.6 %). The two members of "
        "the pair therefore share the same broad routine, with two consistent "
        "differences: calves rested at the surface much less than their mothers "
        "(11.8 % against 20.0 % of the time) and were more than twice as active "
        "in place. Overall, mothers were visible at the surface in 42.8 % of the "
        "intervals and calves in 37.7 %."),
    blank(),
    para(
        "That composition is far from constant across years (Fig. 2). The "
        "clearest change is in exposure rather than in activity: mothers were at "
        "the surface in 70.0 % of the intervals in 1995 and in 74.5 % in 2005, "
        "but in only 24.3 % in 2019 and between 25.4 % and 43.4 % in the last "
        "three years sampled. Rest moved with it. Resting at the surface fell "
        "from 33.7 % of the time in 1995 to about 10 % in 2023 and 2025, while "
        "resting under water rose from 4.2 % in 1995 to a maximum of 29.4 % in "
        "2019, so that the total time spent resting changed far less than the "
        "way in which whales rested. Calves, followed since 2013, show the same "
        "displacement over a shorter span, with surface rest falling from 15.0 % "
        "in 2013 to 6.7 % in 2023. The series is noisy from year to year, and "
        "these are raw proportions, without any correction for the years, sites "
        "and observers behind them; the models in the next sections separate the "
        "part of this drift that is attributable to attacks."),
    blank(),
]

after("write the results sections", R2)


# ---- Figure 3: disturbed, undisturbed and observed behaviour --------------

R3 = [
    heading("Behaviour with and without gull attacks"),
    blank(),
    para(
        "In every year of the series the disturbed and undisturbed "
        "distributions differ in the same direction, and the difference is "
        "large (Fig. 3). A mother under sustained attack rests less, travels "
        "more and, above all, leaves the surface. In 1995 the model puts an "
        "undisturbed mother at the surface 74.2 % of the time and a disturbed "
        "one 24.5 %; between 2015 and 2019 the same contrast runs from about "
        "45 % to about 21 %. Fast travel under water is the category that "
        "absorbs most of the change: it rises from 8.0 % to 62.6 % of the time "
        "in 1995 and from roughly 4 % to between 20 % and 30 % in the recent "
        "years, while rest at the surface falls to about a third of its "
        "undisturbed value throughout. Calves respond in the same direction but "
        "more weakly, and their response is not the same one: their exposure "
        "hardly changes (49.3 % against 43.0 % at the surface in 2013), while "
        "their rest is roughly halved, their fast travel roughly tripled and "
        "their in-place activity increased, from 16.2 % to 26.9 % at the surface "
        "in the same year. A disturbed calf does not so much dive away as stop "
        "resting and start moving in place."),
    blank(),
    para(
        "The two references are not two parallel lines through time: the "
        "undisturbed one moves the most. The mothers' undisturbed time at the "
        "surface fell from 74.2 % in 1995 to between 35 % and 46 % in the last "
        "years, and their undisturbed rest at the surface from 39.0 % to 16.4 %, "
        "whereas the disturbed distribution changed comparatively little "
        "(surface exposure of 24.5 % in 1995 against 23.6 % to 34.1 % in 2021 to "
        "2025). The gap between the two therefore narrowed over three decades, "
        "and it narrowed because the undisturbed whale moved towards the "
        "disturbed one, not the other way round. The observed proportions lie "
        "between the two references and much closer to the undisturbed one, as "
        "expected from an attack occurring in about a fifth of the intervals."),
    blank(),
]

after("Etogram of southern right whales per year", R3)


# ---- Figure 4: the time course -------------------------------------------

R4 = [
    heading("How fast whales react, and how slowly they recover"),
    blank(),
    para(
        "Attacks act quickly. In every period, the estimated disturbance z "
        "covers half the distance to full disturbance within 5 to 10 min of the "
        "first attack, and exceeds 0.9 within 30 min, that is, well inside a "
        "normal follow. Behaviour follows almost as fast: under persistent "
        "attack, the mothers' probability of resting at the surface reaches "
        "90 % of its disturbed value in 20 to 25 min in every period, and the "
        "probability of travelling fast under water in 30 to 45 min (Fig. 4)."),
    blank(),
    para(
        "Recovery is the asymmetric half of the picture, and it has changed. In "
        "1995 a mother recovered as fast as she reacted: z fell below 0.5 within "
        "5 min of the last attack and surface rest was back to 90 % of its "
        "undisturbed value in 25 min. From 2011 onwards recovery is markedly "
        "slower, with a z half-time of about 20 min and 70 to 75 min needed for "
        "surface rest to recover 90 % of the difference, which is longer than a "
        "typical follow. Present-day mothers therefore react to gull attacks "
        "about as fast as they ever did, but return to their routine some three "
        "to four times more slowly, so that a given rate of attacks keeps them "
        "altered for much longer. Calves fall in between, recovering 90 % of "
        "their surface rest in 30 to 40 min. The 1995 comparison should be read "
        "with care: it rests on a single year, and its intervals are the widest "
        "in the figure."),
    blank(),
]

after("Whales behaviour across time", R4)


# ---- Figure 5: the four effects ------------------------------------------

R5 = [
    heading("How much behaviour is altered, and by which part of the effect"),
    blank(),
    para(
        "Whales are strongly reactive to attacks, and less so than they used to "
        "be (Fig. 5). The potential effect, the discrepancy between a fully "
        "disturbed and a fully undisturbed whale of the same year, was 0.585 "
        "(95 % HDI 0.467 to 0.697) for mothers in 1995, meaning that a "
        "continuously harassed mother would have spent nearly 60 % of her time "
        "doing something other than what she would have done undisturbed. From "
        "2004 onwards it settles between 0.28 and 0.43, and declines gently to "
        "0.281 (0.160 to 0.421) in 2025. The calves' potential effect declines "
        "over their shorter series in the same way, from 0.295 (0.206 to 0.390) "
        "in 2013 to 0.183 (0.108 to 0.267) in 2025, and stays below the "
        "mothers' throughout, consistent with the weaker response described "
        "above."),
    blank(),
    para(
        "Only a fraction of that potential is realised by the attack pressure "
        "whales actually experience. The short-term effect stays between 0.066 "
        "and 0.190 for mothers and between 0.058 and 0.112 for calves, roughly a "
        "third of the potential effect in both cases, and it is the only one of "
        "the four that tracks the attacks of the year (Spearman correlation of "
        "0.44 with attack occurrence). The long-term effect goes the other way: "
        "measured against 1995, the mothers' undisturbed behaviour is already "
        "0.113 to 0.271 away from its 1995 counterpart in the 2000s, and reaches "
        "0.410 (0.364 to 0.455) in 2019, 0.381 in 2023 and 0.365 in 2025. It "
        "rises with year (Spearman 0.67) and is uncorrelated with the attack "
        "pressure of the year in which it is measured (0.01). The total effect, "
        "which combines the two, peaks at 0.466 (0.419 to 0.512) in 2019: at "
        "that point a mother spent almost half of her time doing something "
        "different from what a mother in 1995 would have been doing, and the "
        "long-term component alone accounts for 0.410 of those 0.466, against a "
        "short-term effect of 0.133 in the same year. In 2025 the "
        "long-term effect (0.365) is still close to four times the short-term "
        "one (0.098). Whatever changed in these whales' behaviour, it is not "
        "something that ends when the gulls leave."),
    blank(),
]

after("Predictions of whales behaviour under two attack scenarios", R5)


# ---- Figure 6: attacks, behaviour and mortality together -----------------

R6 = [
    heading("Attack pressure, behaviour and calf mortality through the years"),
    blank(),
    para(
        "The proportion of intervals with an attack on either member of the "
        "pair ranges from 0.127 in 1995 to 0.363 in 2011, with no monotone "
        "trend after 2004 (Spearman correlation with year of -0.24). What did "
        "change is who receives the attacks (Fig. 6A). In 1995 the mother was "
        "pecked about twice as often as her calf (0.150 against 0.077 attacks "
        "per interval); from 2004 onwards the ratio is inverted, and the calf "
        "receives between two and five times more attacks than the mother in "
        "every year. Attacks on mothers decline steadily over the series "
        "(Spearman -0.62, from 0.10 to 0.27 attacks per interval up to 2013 down "
        "to 0.04 to 0.11 from 2015 onwards), whereas attacks on calves show no "
        "trend at all (-0.07) and remain between 0.16 and 0.75. Gull harassment "
        "at Península Valdés has not so much intensified as concentrated on "
        "calves."),
    blank(),
    para(
        "The behavioural series below it (Fig. 6B) does not follow the attacks "
        "of the year. Surface rest declines through the series (Spearman -0.50) "
        "and underwater rest increases (0.59), while underwater travel remains "
        "the most common category throughout, and none of these movements "
        "coincides with the peaks and troughs of panel A. The same dissociation "
        "appears in the effects: it is the long-term and total effects that grow "
        "with year (Spearman 0.67 and 0.62) and are uncorrelated with the attack "
        "pressure of the year (0.01 and 0.12), while attack pressure only "
        "explains the short-term effect. A behavioural change that persists "
        "while the attacks on mothers decline is not the immediate reaction to "
        "those attacks, and this is the pattern the discussion has to account "
        "for."),
    blank(),
]

after("Summary of effects of kelp gull attacks", R6)


# ---------------------------------------------------------------- captions

CAPTION_4 = (
    "Figure 4. Predictions of whales behaviour under two attack scenarios in "
    "different periods. Each curve is the behaviour probability predicted by the "
    "model along a follow in which the whale is attacked in every interval "
    "(persistent attacks, solid line), or in which the whale starts fully "
    "disturbed and is never attacked again (attacks cessation, dashed line). "
    "Both scenarios start from the steady state of the opposite condition, so "
    "the persistent-attack curves begin at the undisturbed distribution and the "
    "cessation curves at the disturbed one, and each pair of curves shows how "
    "fast the whale travels from one to the other. Time is measured from the "
    "first attack or from the last one, and the curves are cut at 70 min, about "
    "the length of the longest follows. Panel A shows mothers and panel B "
    "calves; rows are periods, chosen to pool years with similar attack levels "
    "and to keep enough follows in each. Periods contain only the years actually "
    "sampled, so 2014-2019 means 2015 to 2019, 2020-2025 means 2021, 2023 and "
    "2025, and the calves' 2011-2013 is 2013 alone. Multi-year periods are "
    "averaged within each posterior sample before the posterior is summarised. "
    "Lines are posterior means and bands 95 % highest-density intervals."
)

CAPTION_5_INSERT = (
    "The potential effect compares a fully disturbed whale with an undisturbed "
    "one of the same year: it is how much attacks could change behaviour, that "
    "is, the whale's reactivity. The short-term effect compares the behaviour "
    "predicted under the attacks actually received with the undisturbed "
    "behaviour of the same year: how much of that potential the real attack "
    "pressure realises. The long-term effect compares the undisturbed behaviour "
    "of the year with the undisturbed behaviour of 1995: how far the baseline "
    "itself has moved. The total effect compares the behaviour predicted under "
    "the observed attacks with the undisturbed behaviour of 1995, and combines "
    "the previous two."
)

CAPTION_6 = (
    "Figure 6. Time series of attacks, whale behaviour and calf mortality. "
    "(A) Gull attacks on the focal pair: the proportion of intervals of the year "
    "in which the mother, the calf or both were attacked (left axis), and the "
    "mean number of attacks received per interval by the mother and by the calf "
    "(right axis). Attacks are counted on the mothers' dataset, which covers the "
    "whole 1995-2025 series. (B) Proportion of time spent resting and "
    "travelling, each split by exposure, computed over every interval of the "
    "year with the unobserved behaviours imputed by the model. Travel pools slow "
    "and fast travel; in-place activity is not shown, so the four curves do not "
    "add up to one. (C) Calf mortality. Behaviour was not recorded between 1996 "
    "and 2003, nor in 2014, 2020, 2022 and 2024, and the lines break at those "
    "gaps instead of crossing them."
)

after("complete this image caption", [
    para("Claude, proposed caption, replacing the one above: " + CAPTION_4),
    blank(),
])

after("explain again here the meaning of each effect", [
    para("Claude, the sentences the marker asks for, to be inserted in the "
         "caption above right after the first sentence: " + CAPTION_5_INSERT),
    blank(),
])

after("make this caption a bit more complete", [
    para("Claude, proposed caption, replacing the one above: " + CAPTION_6),
    blank(),
    para("Claude, note on the figure file: the two versions of Figure 6 embedded "
         "above are _v01. There is now a _v02 of each, identical except that the "
         "placeholder text in the empty panel C no longer uses an em-dash. The "
         "panel is redrawn anyway as soon as the mortality series arrives, so the "
         "images were not replaced here.", italic=True),
    blank(),
])


# ---------------------------------------------------------------- apply

for anchor, block in inserts:
    assert doc.count(anchor) == 1, "anchor not unique"
    doc = doc.replace(anchor, anchor + block)

open(doc_path, "w", encoding="utf-8").write(doc)


# ---------------------------------------------------------------- verify

ET.fromstring(doc.encode("utf-8"))          # the XML must parse

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

# no em-dashes anywhere in what we wrote (CLAUDE.md rule)
added = [p for p in new_paras if p not in original_paras]
assert not any("—" in p for p in added), "em-dash in the inserted text"

# repack, keeping the original entry order
with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
    for n in names:
        z.write(os.path.join(SRC, n), n)
print("backup:", BACKUP)
print("written:", OUT)
