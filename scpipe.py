"""
scpipe.py  --  shared, TESTED helper functions for the immuno-signature scripts.

`is_noise_gene` is the one with real logic worth testing (the analysis scripts
02/03/04 each carry an inline copy of the same rules; going forward, import from
here so there is one tested source of truth).
"""

import re

_NOISE_PREFIXES = ("RPS", "RPL", "MRPS", "MRPL", "MT-")   # ribosomal + mitochondrial
_NOISE_EXACT = {"MALAT1", "NEAT1", "XIST", "TMSB4X", "TMSB10", "FOS", "FOSB",
                "JUN", "JUNB", "JUND", "EGR1", "HSPA1A", "HSPA1B", "HSPA8", "DNAJB1"}
_CLONE = re.compile(r"^(AC|AL|AP|LINC)\d{3,}")            # uncharacterized clone/lncRNA ids


def is_noise_gene(gene: str) -> bool:
    """True if a gene is a technical artifact (drop it from a signature)."""
    g = gene.upper()
    return (g.startswith(_NOISE_PREFIXES)
            or "." in g                  # versioned clone names like AC016596.1
            or bool(_CLONE.match(g))
            or g in _NOISE_EXACT)
