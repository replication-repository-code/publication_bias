"""Build the analysis samples used in the paper.

Inputs (this folder):
  FinEd_meta_update_effect_sizes_OSF.xlsx   our coding of the newer randomized trials (sheet 'Effect sizes'): the 873 estimates
                                            that enter the analyses (98 papers, 88 trials), one row per estimate, with the
                                            variables the analyses use
  Kaiser_et_al_2022_data_OSF.xlsx           the effect sizes of Kaiser, Lusardi, Menkhoff, and Urban (2022), from their
                                            replication files (sheet 'Sheet1'), with the variables the analyses use
Output: analysis_samples/*.csv, one row per estimate, with columns
  source   'Kaiser et al. (2022)' or 'Update (2026)'
  study    Kaiser study number (K#) or our paper identifier (P###)
  paper    paper identifier (Kaiser studies that are independent samples of one paper share it)
  cluster  trial identifier used for clustering (follow-ups of Kaiser trials share the Kaiser study's identifier)
  yi, sei  Hedges' g (positive = beneficial) and its standard error
  category outcome code (Kaiser A-G; ours A, sk, B-F, Other)
  role     outcome role in the newer trials (R1 knowledge, R2 subjective knowledge, R3 financial behavior,
           DQ decision quality, R4 financial well-being, R5 economic outcomes, R6 mechanisms)
  interv (Classroom/Online/Counseling/Nudge/Other), lowinc, topj, otype (subgroup codes)
Run: python3 build_samples.py
"""
import os, openpyxl, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'analysis_samples'); os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- newer trials
wb = openpyxl.load_workbook(os.path.join(HERE, 'FinEd_meta_update_effect_sizes_OSF.xlsx'), data_only=True)
ws = wb['Effect sizes']; h = [c.value for c in ws[1]]
R = [dict(zip(h, [c.value for c in r])) for r in ws.iter_rows(min_row=2)]
SE, AL = "SE(g)", "Hedges' g aligned incl. authors' view"
ROLE = "Outcome role (decision quality / positive outcomes)"
CATC = "Outcome category (A/sk/B/C/D/E/F/Other)"
# the file holds only the estimates that enter the analyses: randomized trials, estimates included in the analyses and not
# redundant with another estimate of the same trial, each with Hedges' g (positive = beneficial) and a positive standard error
assert len(R) == 873 and all(isinstance(r[AL], (int, float)) and isinstance(r[SE], (int, float)) and r[SE] > 0 for r in R)

# ---------------------------------------------------------------- Kaiser et al. (2022)
k = pd.read_excel(os.path.join(HERE, 'Kaiser_et_al_2022_data_OSF.xlsx'), 'Sheet1'); k = k[k.Study_ID.notna()]
kb = k[k.outcome_code.str[0].isin(list("BCDEFG"))]
kb_raw = kb.copy()   # before removing duplicates (used only by the single-level RoBMA check in folder RoBMA_single_level)
# rows entered twice (same study, estimate, and SE: Bhutoria & Vignoles C1/C6; Mills 2004 D; Modestino B1/B2 x2): keep one
print("Kaiser behavior duplicates dropped:", int(kb.duplicated(['Study_ID', 'hedges_g', 'se_g']).sum()))
kb = kb.drop_duplicates(['Study_ID', 'hedges_g', 'se_g'])
ka = k[k.outcome_code.str[0] == "A"]
# knowledge: one pair with the same study, estimate and SE (study 7, Barcellos et al., table 4 columns 1 and 2): keep one
print("Kaiser knowledge duplicates dropped:", int(ka.duplicated(['Study_ID', 'hedges_g', 'se_g']).sum()))
ka = ka.drop_duplicates(['Study_ID', 'hedges_g', 'se_g'])
# Kaiser 'studies' that are independent samples of one paper
KPAPER = {9: "KP-Batty2015", 10: "KP-Batty2015", 23: "KP-Brugiavini2015", 24: "KP-Brugiavini2015",
          32: "KP-Choi2010", 33: "KP-Choi2010", 34: "KP-Choi2010", 50: "KP-Gibson2012", 51: "KP-Gibson2012", 52: "KP-Gibson2012"}
# newer papers that follow up trials in Kaiser et al.'s data (Peruvian school program, Brazilian high-school program, Ugandan youth clubs)
TO_K = {"P039": "K47", "P040": "K47", "P041": "K47", "P017": "K26", "P047": "K60"}
def kai_otype(r):
    # outcome subgroup from the first letter of Kaiser et al.'s code; their two G1 estimates (study 1) belong to their
    # "Saving & Investing" domain, so they are grouped with saving and investing (D)
    c = str(r.outcome_code)[0]
    return "D" if c == "G" and r.outcome_domain == "Saving & Investing" else c
def kai_rows(kk, role):
    return [dict(source="Kaiser et al. (2022)", study="K%d" % int(r.Study_ID), paper=KPAPER.get(int(r.Study_ID), "K%d" % int(r.Study_ID)),
                 cluster="K%d" % int(r.Study_ID), yi=float(r.hedges_g), sei=float(r.se_g), category=r.outcome_code, role=role,
                 interv=r.intervention_channel,
                 lowinc=("Low-income" if int(r.low_income_dummy) == 1 else "General population"),
                 topj=("Top economics journal" if int(r.top_econ_journal) == 1 else "Other outlet"), otype=kai_otype(r)) for r in kk.itertuples()]
IMAP = {"Classroom course": "Classroom", "Online course": "Online", "Counseling": "Counseling", "Educational nudge": "Nudge"}
def interv(x):
    x = str(x or ""); return IMAP.get(x.split(";")[0].strip(), "Other")
def otype(r):
    ro = str(r.get(ROLE) or "")[:2]
    return str(r[CATC])[0] if ro == "R3" else {"DQ": "Decision quality", "R4": "Financial well-being", "R5": "Economic outcomes", "R6": "Mechanisms", "R1": "A"}.get(ro, ro)
def upd(f):
    return [dict(source="Update (2026)", study=r['Paper ID'], paper=r['Paper ID'], cluster=TO_K.get(r['Paper ID'], "U-" + str(r['Trial ID'])),
                 yi=r[AL], sei=r[SE], category=r[CATC], role=str(r.get(ROLE) or ""), interv=interv(r.get('Intervention type')),
                 lowinc=("General population" if str(r.get('Low-income sample') or '').startswith('No') else "Low-income"),
                 topj=("Top economics journal" if str(r.get('Top journal (your list)') or '').startswith('Yes') else "Other outlet"), otype=otype(r))
            for r in R if f(r)]
roles = lambda *rs: (lambda r: str(r.get(ROLE) or "")[:2] in rs)
kai_beh = kai_rows(kb, "Kaiser behaviour (B-G)"); kai_know = kai_rows(ka, "Kaiser knowledge (A)")
sets = {
    # behavior and downstream outcomes (body of the paper)
    "C_full": kai_beh + upd(roles("R3", "DQ", "R4", "R5", "R6")),
    # without mechanisms; without mechanisms and economic outcomes (web appendix A)
    "C_noR6": kai_beh + upd(roles("R3", "DQ", "R4", "R5")),
    "C_noR5R6": kai_beh + upd(roles("R3", "DQ", "R4")),
    # Kaiser et al.'s behavior categories only (web appendix F)
    "C_R3": kai_beh + upd(roles("R3")),
    # Kaiser et al.'s behavior data alone (web appendix G)
    "K_only": kai_beh,
    # objective and subjective financial knowledge (web appendix R); objective only; Kaiser et al.'s knowledge data alone
    "C_know": kai_know + upd(lambda r: str(r[CATC]) in ("A", "sk")),
    "C_knowobj": kai_know + upd(lambda r: str(r[CATC]) == "A"),
    "K_know": kai_know,
}
# low-income status is a property of the trial: follow-up papers must carry the code of the Kaiser trial they follow up
for n, rows in sets.items():
    df = pd.DataFrame(rows); bad = df.groupby('cluster').lowinc.nunique(); bad = bad[bad > 1]
    assert bad.empty, f"{n}: trials with more than one low-income code: {list(bad.index)}"
for n, rows in sets.items():
    df = pd.DataFrame(rows); df.to_csv(os.path.join(OUT, f"{n}.csv"), index=False)
    print(f"{n:10s} {len(df):5d} estimates  {df.cluster.nunique():4d} trials")
# Kaiser et al.'s behavior estimates before duplicates are removed (458 rows), for folder RoBMA_single_level
pd.DataFrame(kai_rows(kb_raw, "Kaiser behaviour (B-G)")).to_csv(os.path.join(OUT, "kaiser_raw.csv"), index=False)
print(f"kaiser_raw {len(kb_raw):5d} estimates (before removing duplicates)")

# ------------------------------------------------------------------ publication status of the newer papers (Table WA5)
# One row per paper (from the sheet 'Effect sizes'). Status as of our October 2026 check. A working paper that has since appeared
# in a journal counts as published: Bu et al. (P145, Journal of Financial and Quantitative Analysis 2022) and Custodio et al.
# (P031, CEPR DP15294 coded; Journal of Finance 2025).
LATER = {'P145', 'P031'}
pp = pd.DataFrame(R).drop_duplicates('Paper ID').sort_values('Paper ID')
ps = pd.DataFrame({'study': pp['Paper ID'], 'authors': pp['Authors'], 'year': pp['Year'], 'status_coded': pp['Publication status']})
ps['published'] = (ps.status_coded.astype(str).str.startswith('Published') | ps.study.isin(LATER)).astype(int)
ps.to_csv(os.path.join(OUT, 'paper_status.csv'), index=False); print('paper_status.csv:', len(ps), 'papers')
