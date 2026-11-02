"""Complete numeric predictions retain independent evidence semantics."""
from mlvamaps.pipeline import simple_call_rows_from_alleles
from mlvamaps.profile_matching import build_fingerprint
from mlvamaps.models import Locus


def test_zero_depth_exports_keep_selected_count_and_absence():
    alleles = [dict(locus_id="L", called_repeat_count=4, read_depth=0,
                    call_status="LOCUS_DROPOUT")]
    row = simple_call_rows_from_alleles("s", alleles)[0]
    fingerprint, probabilistic = build_fingerprint("s", alleles, [Locus("L")])
    assert row["repeat_count"] == fingerprint[0]["L"] == probabilistic[0]["repeat_count"] == 4
    assert row["present"] == "no" and row["read_depth"] == 0
