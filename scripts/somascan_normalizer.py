"""Deterministic string-level pre-normalization for SomaScan protein targets.

Transforms raw SomaScan (Target, TargetFullName) pairs strictly at the string
level into normalized protein queries and multi-subunit component tables
(seqid_components.tsv) before plain harmonize_many(..., entity="proteins").

No external ID lookups are performed, and no database identifiers (UniProt,
Entrez, HGNC, or KG CURIEs) are ever emitted: every output component_query is
a human protein/gene name string resolved exclusively by BioMapper.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

COMPONENTS_HEADER: tuple[str, ...] = (
    "seq_id",
    "component_key",
    "component_query",
    "source_field",
)

MAPPED_INPUT_HEADER: tuple[str, ...] = (
    "component_key",
    "component_query",
    "seq_id",
    "source_field",
)

_ROMAN_TO_INT: dict[str, int] = {
    "I": 1,
    "II": 2,
    "III": 3,
    "IV": 4,
    "V": 5,
    "VI": 6,
    "VII": 7,
    "VIII": 8,
    "IX": 9,
    "X": 10,
    "XI": 11,
    "XII": 12,
    "XIII": 13,
}

# Single-protein entries that contain slashes, colons, or words like "complex" /
# "family" / "types I and II" and MUST NOT be split into multiple components.
_SINGLE_PROTEIN_GUARDS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^14-3-3 protein (?:beta/alpha|zeta/delta)$", re.I),
    re.compile(r"^14-3-3 protein family$", re.I),
    re.compile(r"^alpha-1-antichymotrypsin complex$", re.I),
    re.compile(r"^signal transducer and activator of transcription 1-alpha/beta$", re.I),
    re.compile(r"^interferon alpha/beta receptor 1$", re.I),
    re.compile(r"^heterogeneous nuclear ribonucleoproteins? a2?/b1?$", re.I),
    re.compile(r"^fractalkine/cx3cl-1$", re.I),
    re.compile(r"^bdnf/nt-3 growth factors receptor$", re.I),
    re.compile(r"^choline/ethanolamine kinase$", re.I),
    re.compile(r"^macrophage scavenger receptor types i and ii$", re.I),
    re.compile(r"^immunoglobulin [adegm]$", re.I),
    re.compile(r"^d-dimer$", re.I),
)

# Short or mnemonic Target strings that collide with unrelated gene symbols, are
# truncated, or are UniProt entry mnemonics rather than official gene symbols.
_UNSAFE_RAW_TARGETS: frozenset[str] = frozenset(
    {
        "14-3-3",
        "14-3-3E",
        "3HAO",
        "3HIDH",
        "41",
        "4EBP2",
        "ACTH",
        "AK1A1",
        "ALK-1",
        "ALT",
        "AMGO2",
        "AMNLS",
        "AMPM2",
        "AN32B",
        "ANGL3",
        "ANGL4",
        "ARGI1",
        "ARI3A",
        "ARMEL",
        "ARP19",
        "ART",
        "ARTS1",
        "ASAHL",
        "ASM3A",
        "ATPO",
        "ATS1",
        "ATS13",
        "ATS15",
        "B7",
        "B7-2",
        "B7-H1",
        "B7-H2",
        "BARK1",
        "BASI",
        "BGH3",
        "BLC",
        "BFL1",
        "BMP-14",
        "BNP-32",
        "BRF-1",
        "BSP",
        "BSSP4",
        "C1QR1",
        "C3a",
        "C3adesArg",
        "C3b",
        "C3d",
        "C4",
        "C4b",
        "C5a",
        "CAD15",
        "CATC",
        "CATE",
        "CATF",
        "CATZ",
        "CBG",
        "CBPE",
        "CD23",
        "CD30",
        "CD39",
        "CHIP",
        "CK-BB",
        "CK-MM",
        "CLC1B",
        "CLC4K",
        "CLC7A",
        "CLM6",
        "CO8A1",
        "CONA1",
        "COX-2",
        "CRDL1",
        "CREL1",
        "CRIS3",
        "CSH",
        "CSK21",
        "CTACK",
        "CTAP-III",
        "CTGF",
        "CYTD",
        "CYTF",
        "CYTN",
        "CYTT",
        "CgA",
        "DAF",
        "DAN",
        "DC-SIGN",
        "DC-SIGNR",
        "DERM",
        "DLC8",
        "DLRB1",
        "DPP2",
        "DR3",
        "DR6",
        "DRAK2",
        "DRG-1",
        "DRR1",
        "DUS3",
        "DcR3",
        "Dtk",
        "EG-VEGF",
        "EGFRvIII",
        "EMAP-2",
        "EMR2",
        "ENA-78",
        "ENTP3",
        "ENTP5",
        "EP15R",
        "EPHAA",
        "EPI",
        "ER",
        "ERAB",
        "ERK-1",
        "FABP",
        "FABPE",
        "FABPL",
        "FAK1",
        "FBLN3",
        "FCG2A",
        "FCG2B",
        "FCG3B",
        "FCGR1",
        "FN1.3",
        "FN1.4",
        "GCKR",
        "GCP-2",
        "GHC2",
        "GI24",
        "GIB",
        "GIIE",
        "GP114",
        "GPDA",
        "GPVI",
        "GV",
        "GX",
        "H2A3",
        "H2B2E",
        "H31",
        "H6ST1",
        "HAI-1",
        "HCC-1",
        "HCC-4",
        "HDGR2",
        "HEMK2",
        "HGH",
        "HGFA",
        "HMG-1",
        "HMGR",
        "HNRPQ",
        "HO-2",
        "HPG-",
        "HPLN1",
        "HSP 27",
        "HSP 40",
        "HSP 60",
        "HSP 70",
        "HSP 90b",
        "HVEM",
        "HXK1",
        "HXK2",
        "I-309",
        "I-TAC",
        "IF4A3",
        "IF4G2",
        "IFN-aA",
        "IFN10",
        "IL-17E",
        "IL-1F6",
        "IL-1F7",
        "IL-1F8",
        "IL-1Ra",
        "IL-22BP",
        "IL-8",
        "ILT-2",
        "ILT-4",
        "IMB1",
        "IMDH1",
        "IMDH2",
        "INGR2",
        "IP-10",
        "IR",
        "JAML1",
        "JNK2",
        "KI2L4",
        "KI3L2",
        "KI3S1",
        "KIRR3",
        "KPCI",
        "KPCT",
        "KREM2",
        "Ku70",
        "LAG-1",
        "LD78-beta",
        "LEAP-1",
        "LEG9",
        "LG3BP",
        "LIMP II",
        "LKHA4",
        "LPPL",
        "LRRT1",
        "LRRT3",
        "LYNB",
        "M2-PK",
        "MAPK2",
        "MAPK5",
        "MASP3",
        "MBL",
        "MCP-1",
        "MCP-2",
        "MCP-3",
        "MCP-4",
        "MDC",
        "MDHC",
        "MEK1",
        "MFGM",
        "MIC-1",
        "MIG",
        "MIP-1a",
        "MIP-3a",
        "MIP-3b",
        "MIP-5",
        "MIS",
        "MK01",
        "MK08",
        "MK11",
        "MK12",
        "MK13",
        "MMEL2",
        "MO2R1",
        "MOZ",
        "MP2K2",
        "MP2K3",
        "MP2K4",
        "MPIF-1",
        "MRCKB",
        "MSP",
        "MYPC1",
        "Met",
        "NAP-2",
        "NCAM-120",
        "NCC27",
        "NET1",
        "NET4",
        "NG36",
        "NKp30",
        "NKp44",
        "NKp46",
        "NLGNX",
        "NOTC2",
        "NPS-PLA2",
        "NRX1B",
        "NRX3B",
        "NSE",
        "NSF1C",
        "NUDC3",
        "NovH",
        "OBCAM",
        "OCAD1",
        "ON",
        "OPG",
        "OX2G",
        "PAFAH",
        "PAK7",
        "PARC",
        "PBEF",
        "PCI",
        "PDE11",
        "PEDF",
        "PERL",
        "PESC",
        "PF-4",
        "PFD5",
        "PGCB",
        "PGP9.5",
        "PGRP-S",
        "PH",
        "PHI",
        "PK3CG",
        "PKC-A",
        "PKC-B-II",
        "PKC-D",
        "PKC-G",
        "PKC-Z",
        "PLPP",
        "PLXB2",
        "PLXC1",
        "PPAC",
        "PPase",
        "PSA",
        "PSA1",
        "PSA2",
        "PSA6",
        "PSD7",
        "PSMA",
        "PSP",
        "PTHrP",
        "PTP-1B",
        "PTP-1C",
        "PUR8",
        "PlGF",
        "QORL1",
        "RANK",
        "RANTES",
        "RAP",
        "RBP",
        "RGM-C",
        "RS3",
        "RS3A",
        "RS7",
        "RUXF",
        "Rb",
        "SAA",
        "SAP",
        "SARP-2",
        "SCGF-alpha",
        "SCGF-beta",
        "SDF-1",
        "SE6L2",
        "SEM5A",
        "SEM6B",
        "SEPR",
        "SH21A",
        "SHP-2",
        "SHPS1",
        "SIG14",
        "SLAF5",
        "SLAF6",
        "SLAF7",
        "SLIK1",
        "SLIK5",
        "SMAC",
        "SNAA",
        "SNP25",
        "SOD",
        "SORC2",
        "SP-D",
        "SPTA2",
        "SRCN1",
        "ST4S6",
        "STRATIFIN",
        "TACI",
        "TAFI",
        "TAJ",
        "TARC",
        "TCCR",
        "TCPTP",
        "TCTP",
        "TECK",
        "TF",
        "TIG2",
        "TIMD3",
        "TMA",
        "TNF-b",
        "TNR4",
        "TRY3",
        "TS",
        "TSG-6",
        "TSP2",
        "TSP4",
        "TWEAK",
        "TWEAKR",
        "TXD12",
        "TrATPase",
        "TrkA",
        "TrkB",
        "TrkC",
        "Tpo",
        "UB2G2",
        "UB2L3",
        "UBC9",
        "UBP25",
        "UCRP",
        "URB",
        "VAV",
        "VEGF",
        "VEGF121",
        "WFKN1",
        "WFKN2",
        "XEDAR",
        "XTP3A",
        "YES",
        "b-ECGF",
        "b-Endorphin",
        "b-NGF",
        "bFGF",
        "bFGF-R",
        "c-Jun",
        "c-Myc",
        "cIAP-2",
        "iC3b",
        "p27Kip1",
        "p53",
        "pTEN",
        "sCD163",
        "sCD4",
        "sE-Selectin",
        "sFRP-3",
        "sICAM-1",
        "sICAM-2",
        "sICAM-3",
        "sICAM-5",
        "sL-Selectin",
        "sLRP1",
        "sLeptin R",
        "sRAGE",
        "sRANKL",
        "sTie-1",
        "sTie-2",
        "suPAR",
        "tPA",
        "tau",
        "uPA",
    }
)

# Named obligate heteromeric complexes whose SomaScan string does not contain a
# ":" or "/" delimiter between subunits. Each entry maps a lowercased Target or
# TargetFullName to its human protein subunit descriptions.
_NAMED_COMPLEX_SUBUNITS: dict[str, tuple[str, ...]] = {
    "interleukin-12": (
        "Interleukin-12 subunit alpha",
        "Interleukin-12 subunit beta",
    ),
    "interleukin-23": (
        "Interleukin-23 subunit alpha",
        "Interleukin-12 subunit beta",
    ),
    "interleukin-27": (
        "Interleukin-27 subunit alpha",
        "Interleukin-27 subunit beta",
    ),
    "follicle stimulating hormone": (
        "Glycoprotein hormones alpha chain",
        "Follitropin subunit beta",
    ),
    "luteinizing hormone": (
        "Glycoprotein hormones alpha chain",
        "Lutropin subunit beta",
    ),
    "thyroid stimulating hormone": (
        "Glycoprotein hormones alpha chain",
        "Thyrotropin subunit beta",
    ),
    "human chorionic gonadotropin": (
        "Glycoprotein hormones alpha chain",
        "Choriogonadotropin subunit beta 3",
    ),
    "complement c1q subcomponent": (
        "Complement C1q subcomponent subunit A",
        "Complement C1q subcomponent subunit B",
        "Complement C1q subcomponent subunit C",
    ),
    "complement component c8": (
        "Complement component C8 alpha chain",
        "Complement component C8 beta chain",
        "Complement component C8 gamma chain",
    ),
    "fibrinogen": (
        "Fibrinogen alpha chain",
        "Fibrinogen beta chain",
        "Fibrinogen gamma chain",
    ),
    "ferritin": (
        "Ferritin heavy chain",
        "Ferritin light chain",
    ),
    "calcineurin": (
        "Protein phosphatase 3 catalytic subunit alpha",
        "Calcineurin subunit B type 1",
    ),
    "calpain i": (
        "Calpain-1 catalytic subunit",
        "Calpain small subunit 1",
    ),
    "hemoglobin": (
        "Hemoglobin subunit alpha-1",
        "Hemoglobin subunit beta",
    ),
    "laminin": (
        "Laminin subunit alpha-1",
        "Laminin subunit beta-1",
        "Laminin subunit gamma-1",
    ),
    "complement c4": (
        "Complement C4-A",
        "Complement C4-B",
    ),
    "complement c4b": (
        "Complement C4-A",
        "Complement C4-B",
    ),
}

# Explicit multi-paralog cross-reactive SomaScan targets (slash-separated genes).
_PARALOG_SLASH_COMPONENTS: dict[str, tuple[str, ...]] = {
    "growth/differentiation factor 11/8": (
        "Growth/differentiation factor 11",
        "Growth/differentiation factor 8",
    ),
    "glycogen synthase kinase-3 alpha/beta": (
        "Glycogen synthase kinase-3 alpha",
        "Glycogen synthase kinase-3 beta",
    ),
    "gro-beta/gamma": (
        "C-X-C motif chemokine 2",
        "C-X-C motif chemokine 3",
    ),
    "heat shock protein hsp 90-alpha/beta": (
        "Heat shock protein HSP 90-alpha",
        "Heat shock protein HSP 90-beta",
    ),
    "rac-alpha/beta/gamma serine/threonine-protein kinase": (
        "RAC-alpha serine/threonine-protein kinase",
        "RAC-beta serine/threonine-protein kinase",
        "RAC-gamma serine/threonine-protein kinase",
    ),
}

# Cleavage products, active enzyme forms, and specific fragments mapped back to
# their parent human protein name (public biochemistry, zero IDs).
_CLEAVAGE_FRAGMENT_TO_PARENT: dict[str, str] = {
    "c3a anaphylatoxin": "Complement C3",
    "c3a anaphylatoxin des arginine": "Complement C3",
    "complement c3b": "Complement C3",
    "complement c3b, inactivated": "Complement C3",
    "complement c3d fragment": "Complement C3",
    "c5a anaphylatoxin": "Complement C5",
    "coagulation factor ixab": "Coagulation factor IX",
    "coagulation factor xa": "Coagulation factor X",
    "activated protein c": "Vitamin K-dependent protein C",
    "thrombin": "Prothrombin",
    "plasmin": "Plasminogen",
    "angiostatin": "Plasminogen",
    "endostatin": "Collagen alpha-1(XVIII) chain",
    "corticotropin": "Pro-opiomelanocortin",
    "beta-endorphin": "Pro-opiomelanocortin",
    "connective tissue-activating peptide iii": "Platelet basic protein",
    "neutrophil-activating peptide 2": "Platelet basic protein",
    "n-terminal pro-bnp": "Natriuretic peptides B",
    "brain natriuretic peptide 32": "Natriuretic peptides B",
    "pituitary adenylate cyclase-activating polypeptide 27": (
        "Pituitary adenylate cyclase-activating polypeptide"
    ),
    "pituitary adenylate cyclase-activating polypeptide 38": (
        "Pituitary adenylate cyclase-activating polypeptide"
    ),
    "somatostatin-28": "Somatostatin",
}

# Public human protein nomenclature canonicalizations (normalized lowercase
# protein name or SomaScan legacy alias -> official HGNC symbol or canonical
# unambiguous human protein name). Contains zero database/KG identifiers.
_NOMENCLATURE_ALIASES: dict[str, str] = {
    "1-phosphatidylinositol 4,5-bisphosphate phosphodiesterase gamma-1": "PLCG1",
    "14-3-3 protein beta/alpha": "YWHAB",
    "14-3-3 protein epsilon": "YWHAE",
    "14-3-3 protein sigma": "SFN",
    "14-3-3 protein theta": "YWHAQ",
    "14-3-3 protein zeta/delta": "YWHAZ",
    "15-hydroxyprostaglandin dehydrogenase [nad(+)]": "HPGD",
    "2'-5'-oligoadenylate synthase 1": "OAS1",
    "26s proteasome non-atpase regulatory subunit 7": "PSMD7",
    "3-hydroxy-3-methylglutaryl-coenzyme a reductase": "HMGCR",
    "3-hydroxyacyl-coa dehydrogenase type-2": "HSD17B10",
    "3-hydroxyanthranilate 3,4-dioxygenase": "HAAO",
    "3-hydroxyisobutyrate dehydrogenase, mitochondrial": "HIBADH",
    "3-phosphoinositide-dependent protein kinase 1": "PDPK1",
    "40s ribosomal protein s3": "RPS3",
    "40s ribosomal protein s3a": "RPS3A",
    "40s ribosomal protein s7": "RPS7",
    "40s ribosomal protein sa": "RPSA",
    "6-phosphogluconate dehydrogenase, decarboxylating": "PGD",
    "60 kda heat shock protein, mitochondrial": "HSPD1",
    "72 kda type iv collagenase": "MMP2",
    "abelson tyrosine-protein kinase 2": "ABL2",
    "acid sphingomyelinase-like phosphodiesterase 3a": "SMPDL3A",
    "acidic leucine-rich nuclear phosphoprotein 32 family member b": "ANP32B",
    "activin receptor type-1b": "ACVR1B",
    "adapter molecule crk": "CRK",
    "adenosine deaminase": "ADA",
    "adenylate kinase isoenzyme 1": "AK1",
    "adenylosuccinate lyase": "ADSL",
    "adhesion g protein-coupled receptor e2": "ADGRE2",
    "adhesion g-protein coupled receptor g5": "ADGRG5",
    "adiponectin": "ADIPOQ",
    "adp-ribosyl cyclase/cyclic adp-ribose hydrolase 1": "CD38",
    "adp-ribosyl cyclase/cyclic adp-ribose hydrolase 2": "BST1",
    "adrenomedullin": "ADM",
    "advanced glycosylation end product-specific receptor": "AGER",
    "afamin": "AFM",
    "aflatoxin b1 aldehyde reductase member 2": "AKR7A2",
    "aggrecan core protein": "ACAN",
    "agouti-related protein": "AGRP",
    "ah receptor-interacting protein": "AIP",
    "alanine aminotransferase 1": "GPT",
    "alcohol dehydrogenase [nadp(+)]": "AKR1A1",
    "alkaline phosphatase, tissue-nonspecific isozyme": "ALPL",
    "allograft inflammatory factor 1": "AIF1",
    "alpha-(1,3)-fucosyltransferase 5": "FUT5",
    "alpha-1-antichymotrypsin": "SERPINA3",
    "alpha-1-antitrypsin": "SERPINA1",
    "alpha-2-antiplasmin": "SERPINF2",
    "alpha-2-hs-glycoprotein": "AHSG",
    "alpha-2-macroglobulin": "A2M",
    "alpha-2-macroglobulin receptor-associated protein": "LRPAP1",
    "alpha-enolase": "ENO1",
    "alpha-fetoprotein": "AFP",
    "alpha-l-iduronidase": "IDUA",
    "alpha-soluble nsf attachment protein": "NAPA",
    "alpha-synuclein": "SNCA",
    "aminoacylase-1": "ACY1",
    "amp-activated protein kinase catalytic subunit alpha-1": "PRKAA1",
    "amp-activated protein kinase catalytic subunit alpha-2": "PRKAA2",
    "amp-activated protein kinase subunit beta-1": "PRKAB1",
    "amp-activated protein kinase subunit beta-2": "PRKAB2",
    "amp-activated protein kinase subunit gamma-1": "PRKAG1",
    "amphoterin-induced protein 2": "AMIGO2",
    "amphiregulin": "AREG",
    "amyloid beta a4 protein": "APP",
    "angiogenin": "ANG",
    "angiopoietin-1 receptor": "TEK",
    "angiotensin-converting enzyme 2": "ACE2",
    "angiotensinogen": "AGT",
    "ankyrin-2": "ANK2",
    "anterior gradient protein 2 homolog": "AGR2",
    "antileukoproteinase": "SLPI",
    "antithrombin-iii": "SERPINC1",
    "apolipoprotein a-i": "APOA1",
    "apolipoprotein b": "APOB",
    "apolipoprotein d": "APOD",
    "apolipoprotein e": "APOE",
    "apolipoprotein l1": "APOL1",
    "apolipoprotein m": "APOM",
    "apoptosis regulator bcl-2": "BCL2",
    "appetite-regulating hormone": "GHRL",
    "arginase-1": "ARG1",
    "aromatic-l-amino-acid decarboxylase": "DDC",
    "artemin": "ARTN",
    "arylsulfatase a": "ARSA",
    "arylsulfatase b": "ARSB",
    "asialoglycoprotein receptor 1": "ASGR1",
    "aspartate aminotransferase, cytoplasmic": "GOT1",
    "at-rich interactive domain-containing protein 3a": "ARID3A",
    "atp synthase subunit beta, mitochondrial": "ATP5F1B",
    "atp synthase subunit o, mitochondrial": "ATP5PO",
    "atp-dependent rna helicase ddx19b": "DDX19B",
    "aurora kinase a": "AURKA",
    "aurora kinase b": "AURKB",
    "azurocidin": "AZU1",
    "b-cell lymphoma 6 protein": "BCL6",
    "b-cell receptor cd22": "CD22",
    "bactericidal permeability-increasing protein": "BPI",
    "baculoviral iap repeat-containing protein 3": "BIRC3",
    "baculoviral iap repeat-containing protein 5": "BIRC5",
    "baculoviral iap repeat-containing protein 7": "BIRC7",
    "basal cell adhesion molecule": "BCAM",
    "basigin": "BSG",
    "bcl-2-like protein 1": "BCL2L1",
    "bcl-2-like protein 2": "BCL2L2",
    "bcl-2-related protein a1": "BCL2A1",
    "bcl2-associated agonist of cell death": "BAD",
    "bdnf/nt-3 growth factors receptor": "NTRK2",
    "beta-2-microglobulin": "B2M",
    "beta-adrenergic receptor kinase 1": "GRK2",
    "beta-ala-his dipeptidase": "CNDP1",
    "beta-nerve growth factor": "NGF",
    "betacellulin": "BTC",
    "bh3-interacting domain death agonist": "BID",
    "biglycan": "BGN",
    "bmp-binding endothelial regulator protein": "BMPER",
    "bone sialoprotein 2": "IBSP",
    "brain natriuretic peptide 32": "NPPB",
    "brain-derived neurotrophic factor": "BDNF",
    "brain-specific serine protease 4": "PRSS22",
    "breast cancer anti-estrogen resistance protein 3": "BCAR3",
    "brevican core protein": "BCAN",
    "brother of cdo": "BOC",
    "c-reactive protein": "CRP",
    "c-type mannose receptor 2": "MRC2",
    "c1q and tumor necrosis factor related protein": "C1QTNF1",
    "calcineurin subunit b type 1": "PPP3R1",
    "calcium/calmodulin-dependent 3',5'-cyclic nucleotide phosphodiesterase 1a": "PDE1A",
    "calcium/calmodulin-dependent protein kinase kinase 1": "CAMKK1",
    "calpain small subunit 1": "CAPNS1",
    "calpain-1 catalytic subunit": "CAPN1",
    "calpastatin": "CAST",
    "calreticulin": "CALR",
    "camp-dependent protein kinase catalytic subunit alpha": "PRKACA",
    "camp-regulated phosphoprotein 19": "ARPP19",
    "camp-specific 3',5'-cyclic phosphodiesterase 4d": "PDE4D",
    "carbohydrate sulfotransferase 15": "CHST15",
    "carbohydrate sulfotransferase 2": "CHST2",
    "carbohydrate sulfotransferase 6": "CHST6",
    "carbonic anhydrase-related protein 10": "CA10",
    "carboxypeptidase b2": "CPB2",
    "carboxypeptidase e": "CPE",
    "cardiotrophin-1": "CTF1",
    "cardiotrophin-like cytokine factor 1": "CLCF1",
    "casein kinase ii subunit alpha": "CSNK2A1",
    "casein kinase ii subunit alpha'": "CSNK2A2",
    "casein kinase ii subunit beta": "CSNK2B",
    "catalase": "CAT",
    "cathepsin a": "CTSA",
    "cathepsin l2": "CTSV",
    "cation-independent mannose-6-phosphate receptor": "IGF2R",
    "ccaat/enhancer-binding protein beta": "CEBPB",
    "cd40 ligand": "CD40LG",
    "cd5 antigen-like": "CD5L",
    "cell adhesion molecule 1": "CADM1",
    "cell adhesion molecule 3": "CADM3",
    "cell adhesion molecule-related/down-regulated by oncogenes": "CDON",
    "cell surface glycoprotein cd200 receptor 1": "CD200R1",
    "cellular tumor antigen p53": "TP53",
    "cerebral dopamine neurotrophic factor": "CDNF",
    "cgmp-dependent 3',5'-cyclic phosphodiesterase": "PDE2A",
    "cgmp-inhibited 3',5'-cyclic phosphodiesterase a": "PDE3A",
    "cgmp-specific 3',5'-cyclic phosphodiesterase": "PDE5A",
    "chitinase-3-like protein 1": "CHI3L1",
    "chitotriosidase-1": "CHIT1",
    "chloride intracellular channel protein 1": "CLIC1",
    "choline/ethanolamine kinase": "CHKB",
    "chordin-like protein 1": "CHRDL1",
    "choriogonadotropin subunit beta 3": "CGB3",
    "chorionic somatomammotropin hormone": "CSH1",
    "chromobox protein homolog 5": "CBX5",
    "chromogranin-a": "CHGA",
    "chymase": "CMA1",
    "ciliary neurotrophic factor": "CNTF",
    "ciliary neurotrophic factor receptor subunit alpha": "CNTFR",
    "clusterin": "CLU",
    "cmrf35-like molecule 6": "CD300C",
    "coactosin-like protein": "COTL1",
    "coiled-coil domain-containing protein 80": "CCDC80",
    "cofilin-1": "CFL1",
    "collagen alpha-1(viii) chain": "COL8A1",
    "collagen alpha-1(xviii) chain": "COL18A1",
    "collagen alpha-1(xxiii) chain": "COL23A1",
    "collagenase 3": "MMP13",
    "comm domain-containing protein 7": "COMMD7",
    "complement c1q subcomponent subunit a": "C1QA",
    "complement c1q subcomponent subunit b": "C1QB",
    "complement c1q subcomponent subunit c": "C1QC",
    "complement c1r subcomponent": "C1R",
    "complement c1s subcomponent": "C1S",
    "complement c4-a": "C4A",
    "complement c4-b": "C4B",
    "complement component 1 q subcomponent-binding protein, mitochondrial": "C1QBP",
    "complement component c1q receptor": "CD93",
    "complement component c8 alpha chain": "C8A",
    "complement component c8 beta chain": "C8B",
    "complement component c8 gamma chain": "C8G",
    "complement decay-accelerating factor": "CD55",
    "complement factor b": "CFB",
    "complement factor d": "CFD",
    "complement factor h": "CFH",
    "complement factor h-related protein 5": "CFHR5",
    "complement factor i": "CFI",
    "connective tissue growth factor": "CCN2",
    "copine-1": "CPNE1",
    "corticosteroid-binding globulin": "SERPINA6",
    "creatine kinase b-type": "CKB",
    "creatine kinase m-type": "CKM",
    "cryptic protein": "CFC1",
    "cyclin-a2": "CCNA2",
    "cyclin-c": "CCNC",
    "cyclin-dependent kinase 1": "CDK1",
    "cyclin-dependent kinase 2": "CDK2",
    "cyclin-dependent kinase 5": "CDK5",
    "cyclin-dependent kinase 5 activator 1": "CDK5R1",
    "cyclin-dependent kinase 8": "CDK8",
    "cyclin-dependent kinase inhibitor 1b": "CDKN1B",
    "cystatin-c": "CST3",
    "cystatin-d": "CST5",
    "cystatin-f": "CST7",
    "cystatin-m": "CST6",
    "cystatin-s": "CST4",
    "cystatin-sa": "CST2",
    "cystatin-sn": "CST1",
    "cysteine and glycine-rich protein 3": "CSRP3",
    "cysteine-rich secretory protein 3": "CRISP3",
    "cysteine-rich with egf-like domain protein 1": "CRELD1",
    "cytochrome c": "CYCS",
    "cytochrome p450 3a4": "CYP3A4",
    "cytokine receptor common subunit gamma": "IL2RG",
    "cytokine receptor-like factor 1": "CRLF1",
    "cytokine receptor-like factor 2": "CRLF2",
    "cytoplasmic protein nck1": "NCK1",
    "cytoplasmic tyrosine-protein kinase bmx": "BMX",
    "cytoskeleton-associated protein 2": "CKAP2",
    "cytosolic non-specific dipeptidase": "CNDP2",
    "cytotoxic and regulatory t-cell molecule": "CRTAM",
    "cytotoxic t-lymphocyte protein 4": "CTLA4",
    "d-dimer": "Fibrin D-dimer",
    "dctp pyrophosphatase 1": "DCTPP1",
    "death-associated protein kinase 2": "DAPK2",
    "decorin": "DCN",
    "dentin matrix acidic phosphoprotein 1": "DMP1",
    "dermatopontin": "DPT",
    "desert hedgehog protein": "DHH",
    "diablo homolog, mitochondrial": "DIABLO",
    "dickkopf-like protein 1": "DKKL1",
    "dipeptidyl peptidase 1": "CTSC",
    "dipeptidyl peptidase 2": "DPP7",
    "discoidin domain-containing receptor 2": "DDR2",
    "dna repair protein rad51 homolog 1": "RAD51",
    "dna topoisomerase 1": "TOP1",
    "dnaj homolog subfamily b member 1": "DNAJB1",
    "down syndrome cell adhesion molecule": "DSCAM",
    "drebrin-like protein": "DBNL",
    "dual 3',5'-cyclic-amp and -gmp phosphodiesterase 11a": "PDE11A",
    "dual specificity protein phosphatase 3": "DUSP3",
    "dual specificity tyrosine-phosphorylation-regulated kinase 3": "DYRK3",
    "dynactin subunit 2": "DCTN2",
    "dynein light chain 1, cytoplasmic": "DYNLL1",
    "dynein light chain roadblock-type 1": "DYNLRB1",
    "e3 sumo-protein ligase pias4": "PIAS4",
    "e3 ubiquitin-protein ligase chip": "STUB1",
    "e3 ubiquitin-protein ligase mdm2": "MDM2",
    "e3 ubiquitin-protein ligase rnf43": "RNF43",
    "e3 ubiquitin-protein ligase znrf3": "ZNRF3",
    "ectodysplasin-a": "EDA",
    "ectonucleoside triphosphate diphosphohydrolase 1": "ENTPD1",
    "ectonucleoside triphosphate diphosphohydrolase 3": "ENTPD3",
    "ectonucleoside triphosphate diphosphohydrolase 5": "ENTPD5",
    "ectonucleotide pyrophosphatase/phosphodiesterase family member 7": "ENPP7",
    "egf-containing fibulin-like extracellular matrix protein 1": "EFEMP1",
    "elafin": "PI3",
    "elongation factor 1-beta": "EEF1B2",
    "endoglin": "ENG",
    "endoplasmic reticulum aminopeptidase 1": "ERAP1",
    "endoplasmic reticulum resident protein 29": "ERP29",
    "endothelial cell-selective adhesion molecule": "ESAM",
    "endothelial cell-specific molecule 1": "ESM1",
    "endothelial monocyte-activating polypeptide 2": "AIMP1",
    "endothelin-converting enzyme 1": "ECE1",
    "enteropeptidase": "TMPRSS15",
    "eotaxin": "CCL11",
    "epidermal growth factor": "EGF",
    "epidermal growth factor receptor": "EGFR",
    "epidermal growth factor receptor substrate 15-like 1": "EPS15L1",
    "epiregulin": "EREG",
    "epithelial discoidin domain-containing receptor 1": "DDR1",
    "erythropoietin": "EPO",
    "erythropoietin receptor": "EPOR",
    "estradiol 17-beta-dehydrogenase 1": "HSD17B1",
    "estrogen receptor": "ESR1",
    "eukaryotic initiation factor 4a-iii": "EIF4A3",
    "eukaryotic translation initiation factor 4 gamma 2": "EIF4G2",
    "eukaryotic translation initiation factor 4e-binding protein 2": "EIF4EBP2",
    "eukaryotic translation initiation factor 4h": "EIF4H",
    "eukaryotic translation initiation factor 5": "EIF5",
    "eukaryotic translation initiation factor 5a-1": "EIF5A",
    "extracellular matrix protein 1": "ECM1",
    "extracellular superoxide dismutase [cu-zn]": "SOD3",
    "fact complex subunit ssrp1": "SSRP1",
    "fatty acid-binding protein, epidermal": "FABP5",
    "fatty acid-binding protein, heart": "FABP3",
    "fatty acid-binding protein, liver": "FABP1",
    "fc receptor-like protein 3": "FCRL3",
    "ferritin heavy chain": "FTH1",
    "ferritin light chain": "FTL",
    "fetuin-b": "FETUB",
    "fibrinogen alpha chain": "FGA",
    "fibrinogen beta chain": "FGB",
    "fibrinogen gamma chain": "FGG",
    "fibronectin": "FN1",
    "fms-related tyrosine kinase 3 ligand": "FLT3LG",
    "focal adhesion kinase 1": "PTK2",
    "follistatin": "FST",
    "follitropin subunit beta": "FSHB",
    "formimidoyltransferase-cyclodeaminase": "FTCD",
    "fractalkine": "CX3CL1",
    "fructose-bisphosphate aldolase a": "ALDOA",
    "g2/mitotic-specific cyclin-b1": "CCNB1",
    "galactoside 3(4)-l-fucosyltransferase": "FUT3",
    "galectin-10": "CLC",
    "galectin-3-binding protein": "LGALS3BP",
    "gamma-enolase": "ENO2",
    "gastrin-releasing peptide": "GRP",
    "gelsolin": "GSN",
    "glia-derived nexin": "SERPINE2",
    "glial cell line-derived neurotrophic factor": "GDNF",
    "glial fibrillary acidic protein": "GFAP",
    "glucagon": "GCG",
    "glucocorticoid receptor": "NR3C1",
    "glucokinase regulatory protein": "GCKR",
    "glucose-6-phosphate isomerase": "GPI",
    "glutamate carboxypeptidase 2": "FOLH1",
    "glutathione s-transferase a3": "GSTA3",
    "glutathione s-transferase p": "GSTP1",
    "glyceraldehyde-3-phosphate dehydrogenase": "GAPDH",
    "glycerol-3-phosphate dehydrogenase [nad(+)], cytoplasmic": "GPD1",
    "glycogen synthase kinase-3 alpha": "GSK3A",
    "glycogen synthase kinase-3 beta": "GSK3B",
    "glycoprotein hormones alpha chain": "CGA",
    "glycylpeptide n-tetradecanoyltransferase 1": "NMT1",
    "granulocyte colony-stimulating factor": "CSF3",
    "granulocyte colony-stimulating factor receptor": "CSF3R",
    "granulocyte-macrophage colony-stimulating factor": "CSF2",
    "granulins": "GRN",
    "granulysin": "GNLY",
    "granzyme a": "GZMA",
    "granzyme b": "GZMB",
    "granzyme h": "GZMH",
    "grb2-related adapter protein 2": "GRAP2",
    "gremlin-1": "GREM1",
    "group 10 secretory phospholipase a2": "PLA2G10",
    "group ib secretory phospholipase a2": "PLA2G1B",
    "group iie secretory phospholipase a2": "PLA2G2E",
    "group v phospholipase a2": "PLA2G5",
    "growth arrest-specific protein 1": "GAS1",
    "growth factor receptor-bound protein 2": "GRB2",
    "growth hormone receptor": "GHR",
    "growth-regulated alpha protein": "CXCL1",
    "gtp-binding nuclear protein ran": "RAN",
    "gtpase kras": "KRAS",
    "haptoglobin": "HP",
    "heat shock 70 kda protein 1a": "HSPA1A",
    "heat shock cognate 71 kda protein": "HSPA8",
    "heat shock protein beta-1": "HSPB1",
    "heat shock protein hsp 90-alpha": "HSP90AA1",
    "heat shock protein hsp 90-beta": "HSP90AB1",
    "heme oxygenase 2": "HMOX2",
    "hemk methyltransferase family member 2": "N6AMT1",
    "hemoglobin subunit beta": "HBB",
    "hemojuvelin": "HJV",
    "hemopexin": "HPX",
    "heparan-sulfate 6-o-sulfotransferase 1": "HS6ST1",
    "heparin cofactor 2": "SERPIND1",
    "heparin-binding egf-like growth factor": "HBEGF",
    "hepatitis a virus cellular receptor 2": "HAVCR2",
    "hepatocyte growth factor": "HGF",
    "hepatocyte growth factor activator": "HGFAC",
    "hepatocyte growth factor receptor": "MET",
    "hepatocyte growth factor-like protein": "MST1",
    "hepatoma-derived growth factor-related protein 2": "HDGFL2",
    "hepcidin": "HAMP",
    "herv-h ltr-associating protein 2": "HHLA2",
    "heterogeneous nuclear ribonucleoprotein a/b": "HNRNPAB",
    "heterogeneous nuclear ribonucleoprotein k": "HNRNPK",
    "heterogeneous nuclear ribonucleoprotein q": "SYNCRIP",
    "heterogeneous nuclear ribonucleoproteins a2/b1": "HNRNPA2B1",
    "hexokinase-1": "HK1",
    "hexokinase-2": "HK2",
    "high affinity cgmp-specific 3',5'-cyclic phosphodiesterase 9a": "PDE9A",
    "high affinity camp-specific 3',5'-cyclic phosphodiesterase 7a": "PDE7A",
    "high affinity immunoglobulin gamma fc receptor i": "FCGR1A",
    "high affinity nerve growth factor receptor": "NTRK1",
    "high mobility group protein b1": "HMGB1",
    "histidine triad nucleotide-binding protein 1": "HINT1",
    "histidine-rich glycoprotein": "HRG",
    "histone acetyltransferase kat6a": "KAT6A",
    "histone acetyltransferase type b catalytic subunit": "HAT1",
    "histone deacetylase 8": "HDAC8",
    "histone h1.2": "H1-2",
    "histone h2a type 3": "H2AW",
    "histone h2a.z": "H2AZ1",
    "histone h2b type 2-e": "H2BC21",
    "histone-lysine n-methyltransferase ehmt2": "EHMT2",
    "homeobox protein nanog": "NANOG",
    "homeodomain-interacting protein kinase 3": "HIPK3",
    "hsp90 co-chaperone cdc37": "CDC37",
    "hyaluronan and proteoglycan link protein 1": "HAPLN1",
    "hypoxia-inducible factor 1-alpha": "HIF1A",
    "icos ligand": "ICOSLG",
    "iduronate 2-sulfatase": "IDS",
    "immunoglobulin alpha fc receptor": "FCAR",
    "immunoglobulin superfamily containing leucine-rich repeat protein 2": "ISLR2",
    "importin subunit alpha-1": "KPNA2",
    "importin subunit beta-1": "KPNB1",
    "induced myeloid leukemia cell differentiation protein mcl-1": "MCL1",
    "inducible t-cell costimulator": "ICOS",
    "inhibin beta a chain": "INHBA",
    "inhibin beta b chain": "INHBB",
    "inhibitor of growth protein 1": "ING1",
    "inorganic pyrophosphatase": "PPA1",
    "inosine-5'-monophosphate dehydrogenase 1": "IMPDH1",
    "inosine-5'-monophosphate dehydrogenase 2": "IMPDH2",
    "insulin": "INS",
    "insulin receptor": "INSR",
    "insulin-degrading enzyme": "IDE",
    "insulin-like growth factor 1 receptor": "IGF1R",
    "insulin-like growth factor i": "IGF1",
    "integrin alpha-1": "ITGA1",
    "integrin alpha-iib": "ITGA2B",
    "integrin alpha-v": "ITGAV",
    "integrin beta-1": "ITGB1",
    "integrin beta-3": "ITGB3",
    "integrin beta-5": "ITGB5",
    "inter-alpha-trypsin inhibitor heavy chain h4": "ITIH4",
    "intercellular adhesion molecule 1": "ICAM1",
    "intercellular adhesion molecule 2": "ICAM2",
    "intercellular adhesion molecule 3": "ICAM3",
    "intercellular adhesion molecule 5": "ICAM5",
    "interferon alpha/beta receptor 1": "IFNAR1",
    "interferon beta": "IFNB1",
    "interferon gamma": "IFNG",
    "interferon gamma receptor 1": "IFNGR1",
    "interferon gamma receptor 2": "IFNGR2",
    "interferon regulatory factor 1": "IRF1",
    "interleukin-1 alpha": "IL1A",
    "interleukin-1 beta": "IL1B",
    "interleukin-1 receptor accessory protein": "IL1RAP",
    "interleukin-1 receptor antagonist protein": "IL1RN",
    "interleukin-1 receptor type 1": "IL1R1",
    "interleukin-1 receptor type 2": "IL1R2",
    "interleukin-1 receptor-like 1": "IL1RL1",
    "interleukin-1 receptor-like 2": "IL1RL2",
    "interleukin-10 receptor subunit alpha": "IL10RA",
    "interleukin-10 receptor subunit beta": "IL10RB",
    "interleukin-11 receptor subunit alpha": "IL11RA",
    "interleukin-12 receptor subunit beta-1": "IL12RB1",
    "interleukin-12 receptor subunit beta-2": "IL12RB2",
    "interleukin-12 subunit alpha": "IL12A",
    "interleukin-12 subunit beta": "IL12B",
    "interleukin-13 receptor subunit alpha-1": "IL13RA1",
    "interleukin-15 receptor subunit alpha": "IL15RA",
    "interleukin-17 receptor a": "IL17RA",
    "interleukin-17 receptor b": "IL17RB",
    "interleukin-17 receptor c": "IL17RC",
    "interleukin-17 receptor d": "IL17RD",
    "interleukin-18 receptor 1": "IL18R1",
    "interleukin-18 receptor accessory protein": "IL18RAP",
    "interleukin-18-binding protein": "IL18BP",
    "interleukin-2 receptor subunit alpha": "IL2RA",
    "interleukin-20 receptor subunit alpha": "IL20RA",
    "interleukin-22 receptor subunit alpha-1": "IL22RA1",
    "interleukin-22 receptor subunit alpha-2": "IL22RA2",
    "interleukin-23 receptor": "IL23R",
    "interleukin-23 subunit alpha": "IL23A",
    "interleukin-27 receptor subunit alpha": "IL27RA",
    "interleukin-27 subunit alpha": "IL27",
    "interleukin-27 subunit beta": "EBI3",
    "interleukin-3 receptor subunit alpha": "IL3RA",
    "interleukin-36 alpha": "IL36A",
    "interleukin-36 beta": "IL36B",
    "interleukin-4 receptor subunit alpha": "IL4R",
    "interleukin-5 receptor subunit alpha": "IL5RA",
    "interleukin-6 receptor subunit alpha": "IL6R",
    "interleukin-6 receptor subunit beta": "IL6ST",
    "interleukin-7 receptor subunit alpha": "IL7R",
    "interstitial collagenase": "MMP1",
    "junctional adhesion molecule-like": "JAML",
    "kallistatin": "SERPINA4",
    "kelch-like ech-associated protein 1": "KEAP1",
    "keratin, type i cytoskeletal 18": "KRT18",
    "kin of irre-like protein 3": "KIRREL3",
    "kinesin-like protein kif23": "KIF23",
    "kininogen-1": "KNG1",
    "kremen protein 2": "KREMEN2",
    "kunitz-type protease inhibitor 1": "SPINT1",
    "kunitz-type protease inhibitor 2": "SPINT2",
    "kynureninase": "KYNU",
    "l-lactate dehydrogenase b chain": "LDHB",
    "l-selectin": "SELL",
    "lactadherin": "MFGE8",
    "lactoperoxidase": "LPO",
    "lactotransferrin": "LTF",
    "lamin-b1": "LMNB1",
    "laminin subunit alpha-1": "LAMA1",
    "laminin subunit beta-1": "LAMB1",
    "laminin subunit gamma-1": "LAMC1",
    "latent-transforming growth factor beta-binding protein 4": "LTBP4",
    "layilin": "LAYN",
    "legumain": "LGMN",
    "leptin": "LEP",
    "leptin receptor": "LEPR",
    "leucine carboxyl methyltransferase 1": "LCMT1",
    "leucine-rich repeat serine/threonine-protein kinase 2": "LRRK2",
    "leucine-rich repeat transmembrane neuronal protein 1": "LRRTM1",
    "leucine-rich repeat transmembrane neuronal protein 3": "LRRTM3",
    "leucine-rich repeat transmembrane protein flrt1": "FLRT1",
    "leucine-rich repeat transmembrane protein flrt2": "FLRT2",
    "leucine-rich repeat transmembrane protein flrt3": "FLRT3",
    "leucine-rich repeats and immunoglobulin-like domains protein 3": "LRIG3",
    "leukemia inhibitory factor receptor": "LIFR",
    "leukocyte surface antigen cd47": "CD47",
    "leukotriene a-4 hydrolase": "LTA4H",
    "ligand-dependent nuclear receptor corepressor-like protein": "LCORL",
    "limbic system-associated membrane protein": "LSAMP",
    "lipopolysaccharide-binding protein": "LBP",
    "lithostathine-1-alpha": "REG1A",
    "low affinity immunoglobulin epsilon fc receptor": "FCER2",
    "low affinity immunoglobulin gamma fc region receptor ii-a": "FCGR2A",
    "low affinity immunoglobulin gamma fc region receptor ii-b": "FCGR2B",
    "low affinity immunoglobulin gamma fc region receptor iii-b": "FCGR3B",
    "low molecular weight phosphotyrosine protein phosphatase": "ACP1",
    "low-density lipoprotein receptor": "LDLR",
    "low-density lipoprotein receptor-related protein 1": "LRP1",
    "low-density lipoprotein receptor-related protein 1b": "LRP1B",
    "low-density lipoprotein receptor-related protein 8": "LRP8",
    "lumican": "LUM",
    "lutropin subunit beta": "LHB",
    "ly6/plaur domain-containing protein 3": "LYPD3",
    "lymphatic vessel endothelial hyaluronic acid receptor 1": "LYVE1",
    "lymphocyte activation gene 3 protein": "LAG3",
    "lymphocyte antigen 86": "LY86",
    "lymphocyte antigen 96": "LY96",
    "lymphotactin": "XCL1",
    "lymphotoxin-alpha": "LTA",
    "lymphotoxin-beta": "LTB",
    "lysosomal protective protein": "CTSA",
    "lysosome membrane protein 2": "SCARB2",
    "lysozyme c": "LYZ",
    "macrophage colony-stimulating factor 1": "CSF1",
    "macrophage colony-stimulating factor 1 receptor": "CSF1R",
    "macrophage mannose receptor 1": "MRC1",
    "macrophage metalloelastase": "MMP12",
    "macrophage migration inhibitory factor": "MIF",
    "macrophage scavenger receptor": "MSR1",
    "macrophage-capping protein": "CAPG",
    "macrophage-stimulating protein receptor": "MST1R",
    "malate dehydrogenase, cytoplasmic": "MDH1",
    "mammaglobin-b": "SCGB2A1",
    "mannan-binding lectin serine protease 1": "MASP1",
    "mannose-binding protein c": "MBL2",
    "mast/stem cell growth factor receptor kit": "KIT",
    "matrilysin": "MMP7",
    "matrix extracellular phosphoglycoprotein": "MEPE",
    "mediator of rna polymerase ii transcription subunit 1": "MED1",
    "megakaryocyte-associated tyrosine-protein kinase": "MATK",
    "melanoma-derived growth regulatory protein": "MIA",
    "membrane frizzled-related protein": "MFRP",
    "membrane metallo-endopeptidase-like 1": "MMEL1",
    "mesothelin": "MSLN",
    "methionine aminopeptidase 1": "METAP1",
    "methionine aminopeptidase 2": "METAP2",
    "methyl-cpg-binding domain protein 4": "MBD4",
    "mhc class i polypeptide-related sequence a": "MICA",
    "mhc class i polypeptide-related sequence b": "MICB",
    "microtubule-associated protein tau": "MAPT",
    "midkine": "MDK",
    "mitochondrial glutamate carrier 2": "SLC25A18",
    "mitochondrial import inner membrane translocase subunit tim14": "DNAJC19",
    "mitogen-activated protein kinase kinase kinase 7": "MAP3K7",
    "moesin": "MSN",
    "mucin-1": "MUC1",
    "mullerian-inhibiting factor": "AMH",
    "myc proto-oncogene protein": "MYC",
    "myeloblastin": "PRTN3",
    "myeloid cell surface antigen cd33": "CD33",
    "myeloperoxidase": "MPO",
    "myoglobin": "MB",
    "myosin-binding protein c, slow-type": "MYBPC1",
    "neuroblastoma suppressor of tumorigenicity 1": "NBL1",
    "neuroligin-4, x-linked": "NLGN4X",
    "n-acetyl-d-glucosamine kinase": "NAGK",
    "n-acetylglucosamine-6-sulfatase": "GNS",
    "n-acylethanolamine-hydrolyzing acid amidase": "NAAA",
    "nad-dependent protein deacetylase sirtuin-2": "SIRT2",
    "nadph--cytochrome p450 reductase": "POR",
    "nascent polypeptide-associated complex subunit alpha": "NACA",
    "natriuretic peptides b": "NPPB",
    "natural cytotoxicity triggering receptor 1": "NCR1",
    "natural cytotoxicity triggering receptor 2": "NCR2",
    "natural cytotoxicity triggering receptor 3": "NCR3",
    "natural killer cell receptor 2b4": "CD244",
    "netrin receptor unc5c": "UNC5C",
    "netrin receptor unc5d": "UNC5D",
    "neural cell adhesion molecule 1": "NCAM1",
    "neural cell adhesion molecule l1": "L1CAM",
    "neural cell adhesion molecule l1-like protein": "CHL1",
    "neurexophilin-1": "NXPH1",
    "neuronal cell adhesion molecule": "NRCAM",
    "neuronal growth regulator 1": "NEGR1",
    "neuropilin-1": "NRP1",
    "neurotrophin-3": "NTF3",
    "neurotrophin-4": "NTF4",
    "neutral ceramidase": "ASAH2",
    "neutrophil collagenase": "MMP8",
    "neutrophil elastase": "ELANE",
    "neutrophil gelatinase-associated lipocalin": "LCN2",
    "nicotinamide phosphoribosyltransferase": "NAMPT",
    "nidogen-1": "NID1",
    "nidogen-2": "NID2",
    "nkg2-d type ii integral membrane protein": "KLRK1",
    "noggin": "NOG",
    "non-histone chromosomal protein hmg-14": "HMGN1",
    "non-receptor tyrosine-protein kinase tyk2": "TYK2",
    "nsfl1 cofactor p47": "NSFL1C",
    "nt-3 growth factor receptor": "NTRK3",
    "nuclear receptor subfamily 1 group d member 1": "NR1D1",
    "nucleoside diphosphate kinase a": "NME1",
    "nucleoside diphosphate kinase b": "NME2",
    "nudc domain-containing protein 3": "NUDCD3",
    "ocia domain-containing protein 1": "OCIAD1",
    "olfactomedin-4": "OLFM4",
    "oncostatin-m": "OSM",
    "opioid-binding protein/cell adhesion molecule": "OPCML",
    "osteocalcin": "BGLAP",
    "osteomodulin": "OMD",
    "osteopontin": "SPP1",
    "ox-2 membrane glycoprotein": "CD200",
    "oxidized low-density lipoprotein receptor 1": "OLR1",
    "p-selectin": "SELP",
    "pancreatic hormone": "PPY",
    "pappalysin-1": "PAPPA",
    "parathyroid hormone": "PTH",
    "parathyroid hormone-related protein": "PTHLH",
    "peptide yy": "PYY",
    "peptidoglycan recognition protein 1": "PGLYRP1",
    "periostin": "POSTN",
    "peroxiredoxin-5, mitochondrial": "PRDX5",
    "peptidyl-prolyl cis-trans isomerase f, mitochondrial": "PPIF",
    "peroxisomal targeting signal 1 receptor": "PEX5",
    "persephin": "PSPN",
    "persulfide dioxygenase ethe1, mitochondrial": "ETHE1",
    "pescadillo homolog": "PES1",
    "phosphatidylethanolamine-binding protein 1": "PEBP1",
    "phosphatidylinositol 3,4,5-trisphosphate 3-phosphatase and dual-specificity protein phosphatase pten": "PTEN",
    "phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit gamma isoform": "PIK3CG",
    "phosphoglucomutase-1": "PGM1",
    "phosphoglycerate kinase 1": "PGK1",
    "phosphoglycerate mutase 1": "PGAM1",
    "phospholipase a2, membrane associated": "PLA2G2A",
    "pigment epithelium-derived factor": "SERPINF1",
    "pilr alpha-associated neural protein": "PIANP",
    "pituitary adenylate cyclase-activating polypeptide": "ADCYAP1",
    "placenta growth factor": "PGF",
    "plasma kallikrein": "KLKB1",
    "plasma protease c1 inhibitor": "SERPING1",
    "plasma serine protease inhibitor": "SERPINA5",
    "plasminogen": "PLG",
    "plasminogen activator inhibitor 1": "SERPINE1",
    "platelet basic protein": "PPBP",
    "platelet endothelial cell adhesion molecule": "PECAM1",
    "platelet factor 4": "PF4",
    "platelet glycoprotein 4": "CD36",
    "platelet glycoprotein ib alpha chain": "GP1BA",
    "platelet glycoprotein vi": "GP6",
    "platelet receptor gi24": "VSIR",
    "platelet-activating factor acetylhydrolase": "PLA2G7",
    "platelet-activating factor acetylhydrolase ib subunit beta": "PAFAH1B2",
    "platelet-derived growth factor c": "PDGFC",
    "platelet-derived growth factor receptor alpha": "PDGFRA",
    "platelet-derived growth factor receptor beta": "PDGFRB",
    "platelet-derived growth factor subunit a": "PDGFA",
    "platelet-derived growth factor subunit b": "PDGFB",
    "pleiotrophin": "PTN",
    "polymeric immunoglobulin receptor": "PIGR",
    "polyubiquitin k48-linked": "UBC",
    "polyubiquitin k63-linked": "UBC",
    "prefoldin subunit 5": "PFDN5",
    "pro-opiomelanocortin": "POMC",
    "programmed cell death 1 ligand 1": "CD274",
    "programmed cell death 1 ligand 2": "PDCD1LG2",
    "prokineticin-1": "PROK1",
    "prolactin": "PRL",
    "prolactin receptor": "PRLR",
    "proliferating cell nuclear antigen": "PCNA",
    "proliferation-associated protein 2g4": "PA2G4",
    "prolyl endopeptidase fap": "FAP",
    "properdin": "CFP",
    "proprotein convertase subtilisin/kexin type 7": "PCSK7",
    "proprotein convertase subtilisin/kexin type 9": "PCSK9",
    "protein amnionless": "AMN",
    "prostaglandin g/h synthase 2": "PTGS2",
    "prostate-specific antigen": "KLK3",
    "protein 4.1": "EPB41",
    "protein deglycase dj-1": "PARK7",
    "protein disulfide-isomerase": "P4HB",
    "protein disulfide-isomerase a3": "PDIA3",
    "protein kinase c alpha type": "PRKCA",
    "protein kinase c beta type": "PRKCB",
    "protein kinase c delta type": "PRKCD",
    "protein kinase c gamma type": "PRKCG",
    "protein kinase c iota type": "PRKCI",
    "protein kinase c theta type": "PRKCQ",
    "protein kinase c zeta type": "PRKCZ",
    "protein lin-7 homolog b": "LIN7B",
    "protein nov homolog": "CCN3",
    "protein phosphatase 3 catalytic subunit alpha": "PPP3CA",
    "protein set": "SET",
    "protein wnt-7a": "WNT7A",
    "protein z-dependent protease inhibitor": "SERPINA10",
    "protein-glutamine gamma-glutamyltransferase e": "TGM3",
    "protein-tyrosine kinase 6": "PTK6",
    "prothrombin": "F2",
    "proto-oncogene tyrosine-protein kinase receptor ret": "RET",
    "proto-oncogene tyrosine-protein kinase src": "SRC",
    "proto-oncogene vav": "VAV1",
    "pulmonary surfactant-associated protein d": "SFTPD",
    "pyridoxal kinase": "PDXK",
    "pyridoxal phosphate phosphatase": "PDXP",
    "pyruvate dehydrogenase (acetyl-transferring) kinase isozyme 1, mitochondrial": "PDK1",
    "pyruvate kinase pkm": "PKM",
    "quinone oxidoreductase-like protein 1": "CRYZL1",
    "rab gdp dissociation inhibitor beta": "GDI2",
    "rac-alpha serine/threonine-protein kinase": "AKT1",
    "rac-beta serine/threonine-protein kinase": "AKT2",
    "rac-gamma serine/threonine-protein kinase": "AKT3",
    "ras gtpase-activating protein 1": "RASA1",
    "ras-related c3 botulinum toxin substrate 1": "RAC1",
    "ras-related c3 botulinum toxin substrate 3": "RAC3",
    "receptor-type tyrosine-protein kinase flt3": "FLT3",
    "regenerating islet-derived protein 4": "REG4",
    "relaxin receptor 1": "RXFP1",
    "renin": "REN",
    "repulsive guidance molecule a": "RGMA",
    "resistin": "RETN",
    "reticulon-4": "RTN4",
    "reticulon-4 receptor": "RTN4R",
    "retinoblastoma-associated protein": "RB1",
    "retinoic acid receptor responder protein 2": "RARRES2",
    "retinol-binding protein 4": "RBP4",
    "rgm domain family member b": "RGMB",
    "ribonuclease h1": "RNASEH1",
    "ribosomal protein s6 kinase alpha-3": "RPS6KA3",
    "ribosomal protein s6 kinase alpha-5": "RPS6KA5",
    "ribosome maturation protein sbds": "SBDS",
    "rna-binding protein 39": "RBM39",
    "s-formylglutathione hydrolase": "ESD",
    "s-phase kinase-associated protein 1": "SKP1",
    "scavenger receptor class f member 1": "SCARF1",
    "scavenger receptor class f member 2": "SCARF2",
    "scavenger receptor cysteine-rich type 1 protein m130": "CD163",
    "sclerostin": "SOST",
    "secreted and transmembrane protein 1": "SECTM1",
    "secreted frizzled-related protein 1": "SFRP1",
    "secreted frizzled-related protein 3": "FRZB",
    "secretin": "SCT",
    "seizure 6-like protein 2": "SEZ6L2",
    "e-selectin": "SELE",
    "serine protease 27": "PRSS27",
    "serine protease htra2, mitochondrial": "HTRA2",
    "serine/threonine-protein kinase 16": "STK16",
    "serine/threonine-protein kinase 17b": "STK17B",
    "serine/threonine-protein kinase chk1": "CHEK1",
    "serine/threonine-protein kinase chk2": "CHEK2",
    "serine/threonine-protein kinase mrck beta": "CDC42BPB",
    "serine/threonine-protein kinase pak 7": "PAK5",
    "serine/threonine-protein kinase pim-1": "PIM1",
    "serine/threonine-protein kinase plk1": "PLK1",
    "serine/threonine-protein kinase receptor r3": "ACVRL1",
    "serine/threonine-protein kinase tbk1": "TBK1",
    "serine/threonine-protein kinase wnk3": "WNK3",
    "serotransferrin": "TF",
    "serum albumin": "ALB",
    "serum amyloid a-1 protein": "SAA1",
    "serum amyloid p-component": "APCS",
    "serum paraoxonase/arylesterase 1": "PON1",
    "sex hormone-binding globulin": "SHBG",
    "sh2 domain-containing protein 1a": "SH2D1A",
    "shc-transforming protein 1": "SHC1",
    "sialoadhesin": "SIGLEC1",
    "small glutamine-rich tetratricopeptide repeat-containing protein alpha": "SGTA",
    "small nuclear ribonucleoprotein f": "SNRPF",
    "small ubiquitin-related modifier 3": "SUMO3",
    "somatostatin": "SST",
    "somatotropin": "GH1",
    "sonic hedgehog protein": "SHH",
    "sorting nexin-4": "SNX4",
    "sparc": "SPARC",
    "sparc-like protein 1": "SPARCL1",
    "sparc-related modular calcium-binding protein 1": "SMOC1",
    "spectrin alpha chain, non-erythrocytic 1": "SPTAN1",
    "sphingosine kinase 1": "SPHK1",
    "sphingosine kinase 2": "SPHK2",
    "spondin-1": "SPON1",
    "stabilin-2": "STAB2",
    "stanniocalcin-1": "STC1",
    "stem cell growth factor-alpha": "CLEC11A",
    "stem cell growth factor-beta": "CLEC11A",
    "stress-induced-phosphoprotein 1": "STIP1",
    "stromal cell-derived factor 1": "CXCL12",
    "stromelysin-1": "MMP3",
    "stromelysin-2": "MMP10",
    "sumo-conjugating enzyme ubc9": "UBE2I",
    "superoxide dismutase [cu-zn]": "SOD1",
    "superoxide dismutase [mn], mitochondrial": "SOD2",
    "synaptosomal-associated protein 25": "SNAP25",
    "syntaxin-1a": "STX1A",
    "t-cell surface glycoprotein cd4": "CD4",
    "t-lymphocyte activation antigen cd80": "CD80",
    "t-lymphocyte activation antigen cd86": "CD86",
    "t-lymphocyte surface antigen ly-9": "LY9",
    "tartrate-resistant acid phosphatase type 5": "ACP5",
    "tata-box-binding protein": "TBP",
    "tenascin": "TNC",
    "teratocarcinoma-derived growth factor 1": "TDGF1",
    "testican-1": "SPOCK1",
    "testican-2": "SPOCK2",
    "tgf-beta receptor type-2": "TGFBR2",
    "tgf-beta-activated kinase 1 and map3k7-binding protein 1": "TAB1",
    "thioredoxin domain-containing protein 12": "TXNDC12",
    "thrombopoietin": "THPO",
    "thrombopoietin receptor": "MPL",
    "thymic stromal lymphopoietin": "TSLP",
    "thymidine kinase, cytosolic": "TK1",
    "thymidylate synthase": "TYMS",
    "thyroglobulin": "TG",
    "thyroid peroxidase": "TPO",
    "thyrotropin subunit beta": "TSHB",
    "thyroxine-binding globulin": "SERPINA7",
    "tissue factor": "F3",
    "tissue factor pathway inhibitor": "TFPI",
    "tissue-type plasminogen activator": "PLAT",
    "transcription factor ap-1": "JUN",
    "transcription factor iiib 90 kda subunit": "BRF1",
    "transforming growth factor beta receptor type 3": "TGFBR3",
    "transforming growth factor beta-1": "TGFB1",
    "transforming growth factor beta-2": "TGFB2",
    "transforming growth factor beta-3": "TGFB3",
    "transforming growth factor-beta-induced protein ig-h3": "TGFBI",
    "transgelin-2": "TAGLN2",
    "transketolase": "TKT",
    "translationally-controlled tumor protein": "TPT1",
    "triosephosphate isomerase": "TPI1",
    "tropomyosin alpha-1 chain": "TPM1",
    "tropomyosin alpha-4 chain": "TPM4",
    "tropomyosin beta chain": "TPM2",
    "troponin t, cardiac muscle": "TNNT2",
    "tumor necrosis factor": "TNF",
    "tumor necrosis factor receptor superfamily member 19l": "RELT",
    "tumor necrosis factor receptor superfamily member 27": "EDA2R",
    "tumor necrosis factor receptor superfamily member 3": "LTBR",
    "tumor necrosis factor receptor superfamily member 6": "FAS",
    "tumor necrosis factor receptor superfamily member edar": "EDAR",
    "tumor necrosis factor ligand superfamily member 6": "FASLG",
    "tumor necrosis factor-inducible gene 6 protein": "TNFAIP6",
    "tumor-associated calcium signal transducer 2": "TACSTD2",
    "tyrosine-protein kinase abl1": "ABL1",
    "tyrosine-protein kinase btk": "BTK",
    "tyrosine-protein kinase csk": "CSK",
    "tyrosine-protein kinase fer": "FER",
    "tyrosine-protein kinase fgr": "FGR",
    "tyrosine-protein kinase fyn": "FYN",
    "tyrosine-protein kinase hck": "HCK",
    "tyrosine-protein kinase jak2": "JAK2",
    "tyrosine-protein kinase lck": "LCK",
    "tyrosine-protein kinase lyn": "LYN",
    "tyrosine-protein kinase receptor tie-1": "TIE1",
    "tyrosine-protein kinase receptor tyro3": "TYRO3",
    "tyrosine-protein kinase src": "SRC",
    "tyrosine-protein kinase syk": "SYK",
    "tyrosine-protein kinase tec": "TEC",
    "tyrosine-protein kinase transmembrane receptor ror1": "ROR1",
    "tyrosine-protein kinase yes": "YES1",
    "tyrosine-protein kinase zap-70": "ZAP70",
    "tyrosine-protein phosphatase non-receptor type 1": "PTPN1",
    "tyrosine-protein phosphatase non-receptor type 11": "PTPN11",
    "tyrosine-protein phosphatase non-receptor type 2": "PTPN2",
    "tyrosine-protein phosphatase non-receptor type 6": "PTPN6",
    "tyrosine-protein phosphatase non-receptor type substrate 1": "SIRPA",
    "ubiquitin carboxyl-terminal hydrolase 25": "USP25",
    "ubiquitin carboxyl-terminal hydrolase isozyme l1": "UCHL1",
    "ubiquitin+1": "UBB",
    "ubiquitin-conjugating enzyme e2 g2": "UBE2G2",
    "ubiquitin-conjugating enzyme e2 l3": "UBE2L3",
    "ubiquitin-conjugating enzyme e2 n": "UBE2N",
    "ubiquitin-fold modifier 1": "UFM1",
    "ubiquitin-fold modifier-conjugating enzyme 1": "UFC1",
    "ubiquitin-like protein isg15": "ISG15",
    "ump-cmp kinase": "CMPK1",
    "urokinase plasminogen activator surface receptor": "PLAUR",
    "urokinase-type plasminogen activator": "PLAU",
    "vacuolar protein sorting-associated protein vta1 homolog": "VTA1",
    "vascular cell adhesion protein 1": "VCAM1",
    "vascular endothelial growth factor a": "VEGFA",
    "vascular endothelial growth factor c": "VEGFC",
    "vascular endothelial growth factor d": "VEGFD",
    "vascular endothelial growth factor receptor 2": "KDR",
    "vascular endothelial growth factor receptor 3": "FLT4",
    "vasoactive intestinal peptide": "VIP",
    "vesicular integral-membrane protein vip36": "LMAN2",
    "vitamin k-dependent protein c": "PROC",
    "vitamin k-dependent protein s": "PROS1",
    "vitronectin": "VTN",
    "von willebrand factor": "VWF",
    "vps10 domain-containing receptor sorcs2": "SORCS2",
    "wap, kazal, immunoglobulin, kunitz and ntr domain-containing protein 1": "WFIKKN1",
    "wap, kazal, immunoglobulin, kunitz and ntr domain-containing protein 2": "WFIKKN2",
    "wnt inhibitory factor 1": "WIF1",
    "wnt1-inducible-signaling pathway protein 1": "CCN4",
    "wnt1-inducible-signaling pathway protein 3": "CCN6",
    "x-linked interleukin-1 receptor accessory protein-like 2": "IL1RAPL2",
    "x-ray repair cross-complementing protein 6": "XRCC6",
    "xaa-pro aminopeptidase 1": "XPNPEP1",
    "interleukin-8": "CXCL8",
    "interleukin-17a": "IL17A",
    "interleukin-17b": "IL17B",
    "interleukin-17d": "IL17D",
    "interleukin-17f": "IL17F",
    "interleukin-25": "IL25",
    "interleukin-37": "IL37",
    "interleukin-2": "IL2",
    "interleukin-3": "IL3",
    "interleukin-4": "IL4",
    "interleukin-5": "IL5",
    "interleukin-6": "IL6",
    "interleukin-7": "IL7",
    "interleukin-9": "IL9",
    "interleukin-10": "IL10",
    "interleukin-11": "IL11",
    "interleukin-13": "IL13",
    "interleukin-16": "IL16",
    "interleukin-19": "IL19",
    "interleukin-20": "IL20",
    "interleukin-22": "IL22",
    "interleukin-24": "IL24",
    "interleukin-34": "IL34",
    "interferon alpha-2": "IFNA2",
    "interferon alpha-7": "IFNA7",
    "interferon alpha-10": "IFNA10",
    "interferon lambda-1": "IFNL1",
    "interferon lambda-2": "IFNL2",
}

# Systematic nomenclature regexes applied after surface normalization.
_FAMILY_REGEX_RULES: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"^c-c motif chemokine (\d+)-like (\d+)$", re.I), r"CCL\1L\2"),
    (re.compile(r"^c-c motif chemokine (\d+)-like$", re.I), r"CCL\1L1"),
    (re.compile(r"^c-c motif chemokine (\d+)$", re.I), r"CCL\1"),
    (re.compile(r"^c-x-c motif chemokine (\d+)$", re.I), r"CXCL\1"),
    (
        re.compile(
            r"^tumor necrosis factor receptor superfamily member (?!edar$|6$|3$|19l$|27$)([a-z0-9]+)$",
            re.I,
        ),
        r"TNFRSF\1",
    ),
    (
        re.compile(
            r"^tumor necrosis factor ligand superfamily member (?!6$)([a-z0-9]+)$",
            re.I,
        ),
        r"TNFSF\1",
    ),
    (re.compile(r"^ephrin type-a receptor (\d+)$", re.I), r"EPHA\1"),
    (re.compile(r"^ephrin type-b receptor (\d+)$", re.I), r"EPHB\1"),
    (re.compile(r"^ephrin-([ab])(\d+)$", re.I), r"EFN\1\2"),
    (re.compile(r"^matrix metalloproteinase-(\d+)$", re.I), r"MMP\1"),
    (
        re.compile(r"^disintegrin and metalloproteinase domain-containing protein (\d+)$", re.I),
        r"ADAM\1",
    ),
    (
        re.compile(
            r"^a disintegrin and metalloproteinase with thrombospondin motifs (\d+)$",
            re.I,
        ),
        r"ADAMTS\1",
    ),
    (re.compile(r"^metalloproteinase inhibitor (\d+)$", re.I), r"TIMP\1"),
    (re.compile(r"^fibroblast growth factor (\d+)$", re.I), r"FGF\1"),
    (re.compile(r"^fibroblast growth factor receptor (\d+)$", re.I), r"FGFR\1"),
    (re.compile(r"^bone morphogenetic protein (\d+)$", re.I), r"BMP\1"),
    (re.compile(r"^bone morphogenetic protein receptor type-(\d)([ab]?)$", re.I), r"BMPR\1\2"),
    (re.compile(r"^insulin-like growth factor-binding protein (\d+)$", re.I), r"IGFBP\1"),
    (re.compile(r"^caspase-(\d+)$", re.I), r"CASP\1"),
    (re.compile(r"^annexin a(\d+)$", re.I), r"ANXA\1"),
    (re.compile(r"^cadherin-(\d+)$", re.I), r"CDH\1"),
    (re.compile(r"^carbonic anhydrase (\d+)$", re.I), r"CA\1"),
    (re.compile(r"^cathepsin ([bdefghsz])$", re.I), r"CTS\1"),
    (re.compile(r"^c-type lectin domain family (\d+) member ([a-z])$", re.I), r"CLEC\1\2"),
    (re.compile(r"^ficolin-(\d+)$", re.I), r"FCN\1"),
    (re.compile(r"^collectin-(\d+)$", re.I), r"COLEC\1"),
    (re.compile(r"^contactin-(\d+)$", re.I), r"CNTN\1"),
    (re.compile(r"^desmoglein-(\d+)$", re.I), r"DSG\1"),
    (re.compile(r"^desmocollin-(\d+)$", re.I), r"DSC\1"),
    (re.compile(r"^galectin-([23489])$", re.I), r"LGALS\1"),
    (re.compile(r"^glypican-(\d+)$", re.I), r"GPC\1"),
    (re.compile(r"^kallikrein-(\d+)$", re.I), r"KLK\1"),
    (re.compile(r"^matrilin-(\d+)$", re.I), r"MATN\1"),
    (re.compile(r"^neuregulin-(\d+)$", re.I), r"NRG\1"),
    (re.compile(r"^netrin-(\d+)$", re.I), r"NTN\1"),
    (re.compile(r"^neurexin-(\d+)-(?:alpha|beta)$", re.I), r"NRXN\1"),
    (re.compile(r"^plexin-([abc])(\d+)$", re.I), r"PLXN\1\2"),
    (re.compile(r"^protein s100-a(\d+)$", re.I), r"S100A\1"),
    (re.compile(r"^protein fam(\d+[a-z]?)$", re.I), r"FAM\1"),
    (re.compile(r"^r-spondin-(\d+)$", re.I), r"RSPO\1"),
    (re.compile(r"^roundabout homolog (\d+)$", re.I), r"ROBO\1"),
    (re.compile(r"^semaphorin-(\d+[a-z])$", re.I), r"SEMA\1"),
    (re.compile(r"^toll-like receptor (\d+)$", re.I), r"TLR\1"),
    (re.compile(r"^thrombospondin-(\d+)$", re.I), r"THBS\1"),
    (re.compile(r"^trefoil factor (\d+)$", re.I), r"TFF\1"),
    (re.compile(r"^sialic acid-binding ig-like lectin (\d+)$", re.I), r"SIGLEC\1"),
    (re.compile(r"^slam family member ([67])$", re.I), r"SLAMF\1"),
    (re.compile(r"^killer cell immunoglobulin-like receptor ([23]d[ls]\d+)$", re.I), r"KIR\1"),
    (re.compile(r"^leukocyte immunoglobulin-like receptor subfamily b member (\d+)$", re.I), r"LILRB\1"),
    (re.compile(r"^mitogen-activated protein kinase (\d+)$", re.I), r"MAPK\1"),
    (re.compile(r"^dual specificity mitogen-activated protein kinase kinase (\d+)$", re.I), r"MAP2K\1"),
    (re.compile(r"^map kinase-activated protein kinase (\d+)$", re.I), r"MAPKAPK\1"),
    (re.compile(r"^calcium/calmodulin-dependent protein kinase type 1([a-z]?)$", re.I), r"CAMK1\1"),
    (re.compile(r"^serine/threonine-protein kinase pak ([1-6])$", re.I), r"PAK\1"),
    (re.compile(r"^peroxiredoxin-(\d+)$", re.I), r"PRDX\1"),
    (re.compile(r"^peptidyl-prolyl cis-trans isomerase ([a-f])$", re.I), r"PPI\1"),
    (re.compile(r"^trypsin-(\d+)$", re.I), r"PRSS\1"),
    (re.compile(r"^gdnf family receptor alpha-(\d+)$", re.I), r"GFRA\1"),
    (re.compile(r"^angiopoietin-(\d+)$", re.I), r"ANGPT\1"),
    (re.compile(r"^angiopoietin-related protein (\d+)$", re.I), r"ANGPTL\1"),
    (re.compile(r"^dickkopf-related protein (\d+)$", re.I), r"DKK\1"),
    (re.compile(r"^delta-like protein (\d+)$", re.I), r"DLL\1"),
    (re.compile(r"^follistatin-related protein (\d+)$", re.I), r"FSTL\1"),
    (re.compile(r"^growth/differentiation factor (?!8$)(\d+)$", re.I), r"GDF\1"),
    (re.compile(r"^mothers against decapentaplegic homolog (\d+)$", re.I), r"SMAD\1"),
    (re.compile(r"^signal transducer and activator of transcription (\d+).*$", re.I), r"STAT\1"),
    (re.compile(r"^complement (?:component )?c(\d+)$", re.I), r"C\1"),
    (re.compile(r"^nkg2d ligand (\d+)$", re.I), r"ULBP\1"),
    (re.compile(r"^proteasome subunit alpha type-(\d+)$", re.I), r"PSMA\1"),
    (re.compile(r"^proteasome activator complex subunit (\d+)$", re.I), r"PSME\1"),
    (re.compile(r"^receptor tyrosine-protein kinase erbb-(\d+)$", re.I), r"ERBB\1"),
    (re.compile(r"^protein jagged-(\d+)$", re.I), r"JAG\1"),
    (re.compile(r"^slit and ntrk-like protein (\d+)$", re.I), r"SLITRK\1"),
    (re.compile(r"^neurogenic locus notch homolog protein (\d+)$", re.I), r"NOTCH\1"),
    (re.compile(r"^integrin alpha-(iib|v|\d+)$", re.I), r"ITGA\1"),
    (re.compile(r"^integrin beta-(\d+)$", re.I), r"ITGB\1"),
    (re.compile(r"^cd(\d+)(?: antigen| glycoprotein)?$", re.I), r"CD\1"),
)

_GREEK_SUBUNIT_MAP: dict[str, str] = {
    "alpha": "A",
    "beta": "B",
    "gamma": "G",
    "delta": "D",
}


@dataclass(frozen=True)
class SeqIdComponent:
    """One normalized mapping component emitted for a SomaScan SeqId."""

    seq_id: str
    component_key: str
    component_query: str
    source_field: str


def _collapse_ws(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip())


def _is_single_protein_guard(target: str, target_full_name: str) -> bool:
    for text in (target.strip(), target_full_name.strip()):
        for pat in _SINGLE_PROTEIN_GUARDS:
            if pat.match(text):
                return True
    return False


def _split_complex_components(
    target: str, target_full_name: str
) -> tuple[list[str], str] | None:
    """Split true multi-subunit complexes and paralog-slash targets.

    Returns ``(raw_component_strings, source_field)`` for multi-component rows,
    or ``None`` when the row represents a single protein (including slash aliases
    for a single protein).
    """
    t = _collapse_ws(target)
    full = _collapse_ws(target_full_name)

    if _is_single_protein_guard(t, full):
        return None

    t_low = t.casefold()
    full_low = full.casefold()

    # 1. AMPK heterotrimer notation (e.g., "AMPK a1b1g1" / "AMPK a2b2g1").
    m_ampk = re.fullmatch(r"ampk\s+a(\d+)b(\d+)g(\d+)", t_low)
    if m_ampk:
        a, b, g = m_ampk.groups()
        return (
            [
                f"AMP-activated protein kinase catalytic subunit alpha-{a}",
                f"AMP-activated protein kinase subunit beta-{b}",
                f"AMP-activated protein kinase subunit gamma-{g}",
            ],
            "Target",
        )

    # 1b. Integrin compact shorthand (e.g. "Integrin aVb3" when TargetFullName lacks ':').
    m_int = re.fullmatch(r"integrin\s+a([a-z0-9]+)b(\d+)", t_low)
    if m_int and ":" not in full:
        a_sub, b_sub = m_int.groups()
        a_label = "V" if a_sub == "v" else ("IIb" if a_sub == "iib" else a_sub)
        return [f"Integrin alpha-{a_label}", f"Integrin beta-{b_sub}"], "Target"

    # 2. C5b, 6 Complex / Complement C5b-C6 complex.
    if re.fullmatch(r"complement c5b-c6 complex", full_low) or t_low == "c5b, 6 complex":
        return ["Complement C5", "Complement component C6"], "TargetFullName"

    # 3. Casein kinase II 2-alpha('):2-beta heterotetramer.
    m_ck2 = re.fullmatch(
        r"casein kinase ii 2-(alpha'?):2-(beta) heterotetramer", full_low
    )
    if m_ck2:
        a_sub, b_sub = m_ck2.groups()
        return (
            [f"Casein kinase II subunit {a_sub}", f"Casein kinase II subunit {b_sub}"],
            "TargetFullName",
        )

    # 4. Lymphotoxin alpha1:beta2 / alpha2:beta1 heterotrimers.
    if re.fullmatch(r"lymphotoxin alpha\d+:beta\d+", full_low):
        return ["Lymphotoxin-alpha", "Lymphotoxin-beta"], "TargetFullName"

    # 5. Named obligate heteromeric complexes.
    if full_low in _NAMED_COMPLEX_SUBUNITS:
        return list(_NAMED_COMPLEX_SUBUNITS[full_low]), "TargetFullName"

    # 6. Explicit multi-paralog slash targets.
    if full_low in _PARALOG_SLASH_COMPONENTS:
        return list(_PARALOG_SLASH_COMPONENTS[full_low]), "TargetFullName"

    # 7. Colon-delimited complexes / heterodimers / fusions in TargetFullName.
    if ":" in full:
        body = re.sub(
            r"\s+(?:complex|heterodimer|fusion)$", "", full, flags=re.I
        ).strip()
        # Fix SomaScan typo "Integrin alpha-I: beta-1 complex" vs "a1b1" in Target.
        body = re.sub(r"\bIntegrin alpha-I\b", "Integrin alpha-1", body)
        parts = [p.strip() for p in body.split(":") if p.strip()]
        if len(parts) >= 2:
            if parts[1].lower().startswith("beta-") and parts[0].lower().startswith(
                "integrin "
            ):
                parts[1] = f"Integrin {parts[1]}"
            return parts, "TargetFullName"

    # 8. Two-gene slash notation when both tokens look like HGNC symbols.
    if re.fullmatch(r"[A-Z][A-Z0-9]{1,11}/[A-Z][A-Z0-9]{1,11}", t):
        return t.split("/"), "Target"

    return None


def _clean_surface_text(text: str) -> str:
    """Normalize SomaScan surface syntax, modifiers, and cleavage fragments."""
    s = _collapse_ws(text)
    s = s.replace("Muellerian", "Mullerian")
    s = re.sub(r"\[Pyruvate dehydrogenase \(acetyl-transferring\)\]", "Pyruvate dehydrogenase (acetyl-transferring)", s)

    # Strip trailing soluble / secreted / isoform / variant / mutation qualifiers.
    patterns = (
        r",\s*soluble form$",
        r",\s*soluble$",
        r",\s*secreted form$",
        r"\s*\(isoform [^)]+\)$",
        r"\s*\(splice variant [^)]+\)$",
        r",\s*\d+\s*kDa\s*isoform$",
        r",\s*isoform\s+[A-Za-z0-9-]+$",
        r"\s+[Ii]soform\s+(?:alpha|beta|gamma|delta|A|B|C|\d+)$",
        r"\s+variant\s+[IVX0-9]+$",
        r"\s+N-product$",
        r",\s*truncated mutation for\s+.*$",
        r",\s*Mixed Type$",
        r"\s+types I and II$",
        r"\s+Fragment\s+\d+$",
        r",\s*inactivated$",
    )
    for pat in patterns:
        s = re.sub(pat, "", s).strip()

    s = re.sub(r"\s+g-chain dimer$", " gamma chain", s)
    s = re.sub(r"^Alpha-1-antichymotrypsin complex$", "Alpha-1-antichymotrypsin", s)
    s = re.sub(r"^Ck-beta-8-1$", "C-C motif chemokine 23", s, flags=re.I)

    low = s.casefold()
    if low in _CLEAVAGE_FRAGMENT_TO_PARENT:
        s = _CLEAVAGE_FRAGMENT_TO_PARENT[low]

    return _collapse_ws(s)


def _apply_family_rules(cleaned: str) -> str | None:
    """Apply systematic public protein-family nomenclature rules."""
    low = cleaned.casefold()

    # Special multi-group / irregular family cases before uniform regexes.
    if low == "slam family member 5":
        return "CD84"
    if low == "galectin-10":
        return "CLC"
    if low == "growth/differentiation factor 8":
        return "MSTN"
    if low == "cathepsin a":
        return "CTSA"
    if low == "cathepsin l2":
        return "CTSV"
    if low == "c-type lectin domain family 4 member k":
        return "CD207"
    if low == "c-type lectin domain family 4 member m":
        return "CLEC4M"
    if low == "serine/threonine-protein kinase pak 7":
        return "PAK5"
    if low == "tumor necrosis factor receptor superfamily member 19l":
        return "RELT"
    if low == "tumor necrosis factor receptor superfamily member 27":
        return "EDA2R"
    if low == "tumor necrosis factor receptor superfamily member 6":
        return "FAS"
    if low == "tumor necrosis factor receptor superfamily member 3":
        return "LTBR"
    if low == "tumor necrosis factor receptor superfamily member edar":
        return "EDAR"
    if low == "tumor necrosis factor ligand superfamily member 6":
        return "FASLG"
    if low == "tryptase beta-2":
        return "TPSB2"
    if low == "tryptase gamma":
        return "TPSG1"
    if low == "junctional adhesion molecule b":
        return "JAM2"
    if low == "junctional adhesion molecule c":
        return "JAM3"
    if low == "carbonic anhydrase 6":
        return "CA6"

    m_camk2 = re.fullmatch(
        r"calcium/calmodulin-dependent protein kinase type ii subunit (alpha|beta|delta|gamma)",
        low,
    )
    if m_camk2:
        return f"CAMK2{_GREEK_SUBUNIT_MAP[m_camk2.group(1)]}"

    m_integrin_a = re.fullmatch(r"integrin alpha-(iib|v|\d+)", low)
    if m_integrin_a:
        sub = m_integrin_a.group(1).upper()
        return f"ITGA{'2B' if sub == 'IIB' else sub}"

    m_coag = re.fullmatch(r"coagulation factor ([ivx]+)", low)
    if m_coag:
        num = _ROMAN_TO_INT.get(m_coag.group(1).upper())
        if num is not None:
            return f"F{num}"

    for pat, repl in _FAMILY_REGEX_RULES:
        m = pat.fullmatch(cleaned)
        if m:
            out = pat.sub(repl, cleaned)
            # Uppercase letters in symbol outputs (e.g. tnfrsf10b -> TNFRSF10B).
            return out.upper()
    return None


def _canonicalize_component(
    component_text: str,
    *,
    target: str = "",
    target_full_name: str = "",
    is_split: bool = False,
) -> tuple[str, str]:
    """Canonicalize one protein or subunit description into a clean query string."""
    t_clean = _collapse_ws(target)
    f_clean = _collapse_ws(target_full_name)

    # Disambiguate cases where TargetFullName omits isoform/group specificity.
    if not is_split:
        if t_clean == "GIB" and f_clean == "Phospholipase A2":
            return "PLA2G1B", "Target+TargetFullName"
        if t_clean == "GV" and f_clean == "Calcium-dependent phospholipase A2":
            return "PLA2G5", "Target+TargetFullName"
        if t_clean == "c-Jun" and f_clean == "Transcription factor AP-1":
            return "JUN", "Target"

    cleaned = _clean_surface_text(component_text)
    low = cleaned.casefold()
    default_src = "Target" if (is_split and component_text == t_clean) or not f_clean else "TargetFullName"

    if low in _NOMENCLATURE_ALIASES:
        return _NOMENCLATURE_ALIASES[low], default_src

    fam = _apply_family_rules(cleaned)
    if fam is not None:
        return fam, default_src

    # For unmapped single-protein rows, prefer Target when it is an unambiguous
    # uppercase HGNC-style symbol and not an unsafe shorthand or UniProt mnemonic.
    if (
        not is_split
        and t_clean not in _UNSAFE_RAW_TARGETS
        and re.fullmatch(r"[A-Z][A-Z0-9]{1,11}(?:-[A-Z0-9]+)?", t_clean)
    ):
        return t_clean, "Target"

    return cleaned, ("Target" if component_text == t_clean else "TargetFullName")


def normalize_somascan_row(
    seq_id: str, target: str, target_full_name: str
) -> list[SeqIdComponent]:
    """Normalize one SomaScan target row into one or more component queries."""
    sid = seq_id.strip()
    t = _collapse_ws(target)
    f = _collapse_ws(target_full_name)

    split = _split_complex_components(t, f)
    if split is not None:
        raw_parts, split_source = split
        components: list[SeqIdComponent] = []
        for idx, raw_part in enumerate(raw_parts, start=1):
            query, part_src = _canonicalize_component(
                raw_part, target=t, target_full_name=f, is_split=True
            )
            src = split_source if split_source else part_src
            components.append(
                SeqIdComponent(
                    seq_id=sid,
                    component_key=f"{sid}__{idx}",
                    component_query=query,
                    source_field=src,
                )
            )
        return components

    primary_text = f or t
    query, src = _canonicalize_component(
        primary_text, target=t, target_full_name=f, is_split=False
    )
    return [
        SeqIdComponent(
            seq_id=sid,
            component_key=sid,
            component_query=query,
            source_field=src,
        )
    ]


def normalize_somascan_table(
    input_tsv: Path, *, limit: int | None = None
) -> tuple[list[SeqIdComponent], dict[str, int]]:
    """Normalize all rows (or first ``limit`` rows) of a SomaScan target TSV."""
    with input_tsv.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if limit is not None:
        rows = rows[:limit]

    all_components: list[SeqIdComponent] = []
    n_split = 0
    n_single = 0
    for row in rows:
        comps = normalize_somascan_row(
            row["SeqId"], row["Target"], row["TargetFullName"]
        )
        if len(comps) > 1:
            n_split += 1
        else:
            n_single += 1
        all_components.extend(comps)

    stats = {
        "n_input_seqids": len(rows),
        "n_single_seqids": n_single,
        "n_split_seqids": n_split,
        "n_total_components": len(all_components),
        "n_unique_queries": len({c.component_query for c in all_components}),
    }
    return all_components, stats


def write_normalized_inputs(
    input_tsv: Path,
    *,
    components_tsv: Path,
    mapped_tsv: Path,
    limit: int | None = None,
) -> dict[str, int]:
    """Write ``seqid_components.tsv`` and the normalized TwinsUK TSV for mapping."""
    components, stats = normalize_somascan_table(input_tsv, limit=limit)
    components_tsv.parent.mkdir(parents=True, exist_ok=True)
    mapped_tsv.parent.mkdir(parents=True, exist_ok=True)

    with components_tsv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(COMPONENTS_HEADER)
        for comp in components:
            writer.writerow(
                [comp.seq_id, comp.component_key, comp.component_query, comp.source_field]
            )

    with mapped_tsv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(MAPPED_INPUT_HEADER)
        for comp in components:
            writer.writerow(
                [comp.component_key, comp.component_query, comp.seq_id, comp.source_field]
            )

    return stats
