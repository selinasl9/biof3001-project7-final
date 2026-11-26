#!/usr/bin/env python3
"""
mosquito.py

Recognise a mosquito species in free text + pool genera into the groups the figures use.

Used by: hosts/lookup_host_species.py, tables/build_contig_annotation.py,
coi/identify_host_from_coi.py
"""

import re

# Set the genera, the name patterns, the sample-name codes + the genus groups
# Ochlerotatus is now treated as a subgenus of Aedes, so the two are pooled
GENERA = ["Culex", "Aedes", "Ochlerotatus", "Anopheles", "Armigeres", "Mansonia",
          "Culiseta", "Coquillettidia", "Uranotaenia", "Verrallina", "Lutzia",
          "Tripteroides", "Psorophora", "Haemagogus", "Sabethes", "Toxorhynchites",
          "Orthopodomyia", "Wyeomyia", "Deinocerites", "Mimomyia", "Ficalbia"]

NAME = re.compile(r"\b(" + "|".join(GENERA) + r")\b(?:\s+([a-z]{3,}))?")

NOT_SPECIES = {"sp", "spp", "species", "mosquito", "mosquitoes", "pool", "pools",
               "female", "males", "larvae", "adult", "and", "from", "metagenome"}

NAME_CODES = {"cx": "Culex", "cu": "Culex", "ae": "Aedes", "ad": "Aedes",
              "an": "Anopheles", "ar": "Armigeres", "cs": "Culiseta",
              "ma": "Mansonia", "ur": "Uranotaenia"}
CODE_AT_START = re.compile(r"\s*([A-Za-z]{2})[-_. ]")

GENUS_GROUPS = {"Culex": "Culex", "Aedes": "Aedes", "Ochlerotatus": "Aedes"}

NOT_A_SPECIES = {"", "mosquito metagenome", "Culicidae"}


def find_species(text):
    """Return the species in free text, the genus alone if there is no species, else None."""
    match = NAME.search(text or "")
    if not match:
        return None

    genus, species = match.group(1), match.group(2)
    if species and species.lower() in NOT_SPECIES:
        species = None
    return f"{genus} {species}" if species else genus


def genus_from_code(text):
    """Return a genus from a two-letter code at the front of a sample name, or None."""
    match = CODE_AT_START.match(text or "")
    if not match:
        return None
    return NAME_CODES.get(match.group(1).lower())


def genus_and_group(species):
    """Return (genus, group) for a species name, both "unknown" when there is no name."""
    if not species or species in NOT_A_SPECIES:
        return "unknown", "unknown"

    genus = species.split()[0]
    return genus, GENUS_GROUPS.get(genus, "unknown")
