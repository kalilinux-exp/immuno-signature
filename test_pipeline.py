"""
test_pipeline.py  --  check that the code's functions actually behave correctly.

This is how you "check a function" in Python: feed it known inputs, assert the
output is what it MUST be. If any assert fails, Python raises an error pointing
to the broken case. No test framework needed.

Run:  .venv\\Scripts\\python.exe test_pipeline.py
"""

from scpipe import is_noise_gene


def test_drops_known_noise():
    """Ribosomal, mito, stress, and clone/version ids must be dropped."""
    for g in ["RPS14", "RPL13A", "MT-CO3", "MALAT1", "FOS", "AC016596.1", "LINC00115"]:
        assert is_noise_gene(g), f"FAIL: should have dropped {g}"


def test_keeps_real_biology():
    """Real immune genes must be kept."""
    for g in ["TCF7", "CCR7", "IL7R", "CD8A", "CD28", "LAT", "GZMB"]:
        assert not is_noise_gene(g), f"FAIL: should have kept {g}"


def test_tricky_edge_case():
    """ACE2 starts with 'AC' but is a REAL gene, not a clone id -> must be kept."""
    assert not is_noise_gene("ACE2"), "FAIL: ACE2 is real, must not be dropped"


if __name__ == "__main__":
    test_drops_known_noise()
    test_keeps_real_biology()
    test_tricky_edge_case()
    print("All tests passed - OK")
