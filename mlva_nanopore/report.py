from __future__ import annotations

import html
import math
from pathlib import Path

from .models import Locus


def _safe(value) -> str:
    return html.escape(str(value))


def _called_count(value) -> int | None:
    if value in ("", None):
        return None
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return None


def _amplicon_size(locus: Locus, repeat_count: int | None) -> int | None:
    if repeat_count is None:
        return None
    return (
        len(locus.forward_primer)
        + len(locus.left_flank_sequence)
        + len(locus.repeat_motif) * repeat_count
        + len(locus.right_flank_sequence)
        + len(locus.reverse_primer)
    )


def _best_profile(match_rows: list[dict], profiles: list[dict]) -> dict | None:
    if not match_rows or not profiles:
        return None
    best_id = match_rows[0].get("best_profile_id")
    for profile in profiles:
        if profile.get("profile_id") == best_id:
            return profile
    return None


def _gel_svg(
    sample_id: str,
    loci: list[Locus],
    allele_rows: list[dict],
    best_profile: dict | None,
    asv_rows: list[dict],
) -> str:
    allele_by_locus = {row["locus_id"]: row for row in allele_rows}
    asv_by_locus: dict[str, list[dict]] = {}
    for row in asv_rows:
        asv_by_locus.setdefault(row["locus_id"], []).append(row)

    query_bands = []
    reference_bands = []
    for locus in loci:
        locus_asvs = asv_by_locus.get(locus.locus_id, [])
        if locus_asvs:
            for asv in locus_asvs:
                query_size = _amplicon_size(locus, _called_count(asv.get("repeat_count")))
                if query_size:
                    support = _called_count(asv.get("support_reads")) or 0
                    frequency = float(asv.get("frequency") or 0)
                    query_bands.append((locus.locus_id, query_size, support, frequency))
        else:
            allele = allele_by_locus.get(locus.locus_id, {})
            query_size = _amplicon_size(locus, _called_count(allele.get("called_repeat_count")))
            if query_size:
                support = _called_count(allele.get("read_depth")) or 0
                query_bands.append((locus.locus_id, query_size, support, 1.0 if support else 0.0))
        if best_profile:
            reference_size = _amplicon_size(locus, _called_count(best_profile.get(locus.locus_id)))
            if reference_size:
                reference_bands.append((locus.locus_id, reference_size, 0, 1.0))

    all_sizes = [size for _name, size, _support, _frequency in query_bands + reference_bands]
    if not all_sizes:
        return '<p class="terminal-note">No callable loci were available for gel rendering.</p>'

    marker_sizes = [2000, 1500, 1000, 700, 500, 300, 200, 100, 50]
    all_sizes.extend(marker_sizes)
    max_size = max(all_sizes)
    min_size = max(20, min(all_sizes))
    gel_top = 78
    gel_height = 340

    def y_for_size(size: int) -> float:
        log_max = math.log10(max_size)
        log_min = math.log10(min_size)
        if log_max == log_min:
            return gel_top + gel_height / 2
        return gel_top + ((log_max - math.log10(size)) / (log_max - log_min)) * gel_height

    max_support = max((support for _locus_id, _size, support, _frequency in query_bands), default=1)

    def query_intensity(support: int) -> tuple[float, float]:
        if max_support <= 0:
            return 0.25, 5.0
        scaled = math.sqrt(support / max_support)
        return 0.25 + 0.72 * scaled, 4.0 + 8.0 * scaled

    def query_bands_svg(bands: list[tuple[str, int, int, float]], x: int) -> str:
        rects = []
        for locus_id, size, support, frequency in bands:
            y = y_for_size(size)
            width = 66
            opacity, height = query_intensity(support)
            rects.append(
                f'<g><title>{_safe(locus_id)}: {size} bp; {support} reads; frequency {frequency:.3f}</title>'
                f'<rect class="query-band" x="{x - width / 2:.1f}" y="{y - height / 2:.1f}" '
                f'width="{width}" height="{height:.1f}" rx="2" opacity="{opacity:.3f}" />'
                f'<text class="band-label" x="{x + 43}" y="{y + 3:.1f}">{_safe(locus_id)}</text></g>'
            )
        return "\n".join(rects)

    def reference_bands_svg(bands: list[tuple[str, int, int, float]], x: int) -> str:
        rects = []
        for locus_id, size, _support, _frequency in bands:
            y = y_for_size(size)
            width = 66
            rects.append(
                f'<g><title>{_safe(locus_id)} reference: {size} bp</title>'
                f'<rect class="reference-band" x="{x - width / 2:.1f}" y="{y - 3.5:.1f}" '
                f'width="{width}" height="7" rx="2" />'
                f'<text class="band-label reference-label" x="{x + 43}" y="{y + 3:.1f}">{_safe(locus_id)}</text></g>'
            )
        return "\n".join(rects)

    marker = "\n".join(
        f'<g><rect class="marker-band" x="82" y="{y_for_size(size) - 2.5:.1f}" width="56" height="5" rx="2" />'
        f'<text class="marker-label" x="30" y="{y_for_size(size) + 3:.1f}">{size}</text></g>'
        for size in marker_sizes
        if min_size <= size <= max_size
    )
    reference_name = best_profile.get("profile_id", "best reference") if best_profile else "no reference"
    return f"""
<figure class="gel-panel" aria-label="Generated agarose gel comparison">
  <svg viewBox="0 0 720 500" role="img" aria-labelledby="gel-title gel-desc">
    <title id="gel-title">Generated MLVA agarose gel comparison</title>
    <desc id="gel-desc">Marker, query sample, and best matching reference profile bands estimated from VNTR amplicon sizes.</desc>
    <defs>
      <filter id="glow"><feGaussianBlur stdDeviation="2.4" result="blur"/><feMerge><feMergeNode in="blur"/><feMergeNode in="SourceGraphic"/></feMerge></filter>
      <linearGradient id="gel-bg" x1="0" x2="1" y1="0" y2="1">
        <stop offset="0%" stop-color="#24103f"/>
        <stop offset="55%" stop-color="#0a0928"/>
        <stop offset="100%" stop-color="#060816"/>
      </linearGradient>
      <pattern id="scanlines" width="6" height="6" patternUnits="userSpaceOnUse">
        <rect width="6" height="1" fill="#ffffff" opacity="0.045"/>
      </pattern>
    </defs>
    <rect class="gel-frame" x="18" y="18" width="684" height="458" rx="8"/>
    <rect x="42" y="54" width="636" height="388" rx="6" fill="url(#gel-bg)"/>
    <rect x="42" y="54" width="636" height="388" rx="6" fill="url(#scanlines)"/>
    <line class="well-line" x1="62" y1="72" x2="658" y2="72"/>
    <rect class="well" x="80" y="58" width="60" height="20" rx="3"/>
    <rect class="well" x="288" y="58" width="72" height="20" rx="3"/>
    <rect class="well" x="496" y="58" width="72" height="20" rx="3"/>
    <text class="lane-title" x="110" y="466">LADDER</text>
    <text class="lane-title" x="324" y="466">{_safe(sample_id)}</text>
    <text class="lane-title" x="532" y="466">{_safe(reference_name)}</text>
    <text class="gel-legend" x="324" y="42">query band intensity = fragment read support</text>
    {marker}
    {query_bands_svg(query_bands, 324)}
    {reference_bands_svg(reference_bands, 532)}
  </svg>
  <figcaption>Generated gel image: band position shows estimated fragment length; query band brightness and thickness scale with reads supporting that fragment. Horizontal alignment between query and reference bands indicates matching VNTR fragment lengths.</figcaption>
</figure>
"""


def write_report(
    outdir: str | Path,
    sample_id: str,
    allele_rows: list[dict],
    novelty_rows: list[dict],
    loci: list[Locus] | None = None,
    match_rows: list[dict] | None = None,
    profiles: list[dict] | None = None,
    asv_rows: list[dict] | None = None,
) -> None:
    outdir = Path(outdir)
    loci = loci or []
    match_rows = match_rows or []
    profiles = profiles or []
    asv_rows = asv_rows or []
    passed = sum(1 for row in allele_rows if row.get("call_status") == "PASS")
    low_depth = sum(1 for row in allele_rows if row.get("call_status") == "LOW_DEPTH")
    dropout = sum(1 for row in allele_rows if row.get("call_status") == "LOCUS_DROPOUT")
    novelty = novelty_rows[0] if novelty_rows else {}
    best_match = match_rows[0] if match_rows else {}
    best_profile = _best_profile(match_rows, profiles)
    gel = _gel_svg(sample_id, loci, allele_rows, best_profile, asv_rows)
    rows = "\n".join(
        f"<tr><td>{_safe(row['locus_id'])}</td><td>{_safe(row['called_repeat_count'])}</td>"
        f"<td>{_safe(row['posterior_probability'])}</td><td>{_safe(row['read_depth'])}</td>"
        f"<td>{_safe(row['call_status'])}</td></tr>"
        for row in allele_rows
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>MLVA Nanopore Report - {sample_id}</title>
  <style>
    :root {{
      color-scheme: dark;
      --screen: #07150f;
      --panel: #0b2117;
      --phosphor: #62ff9b;
      --amber: #ffc857;
      --magenta: #ff5fc8;
      --cyan: #5fe8ff;
      --muted: #8dd7aa;
      --line: #245c3a;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background:
        radial-gradient(circle at 50% -15%, rgba(98, 255, 155, 0.14), transparent 34rem),
        linear-gradient(180deg, #030705 0%, #07150f 55%, #030705 100%);
      color: var(--phosphor);
      font-family: "Courier New", Courier, monospace;
      letter-spacing: 0;
    }}
    body::before {{
      content: "";
      position: fixed;
      inset: 0;
      pointer-events: none;
      background: repeating-linear-gradient(180deg, rgba(255,255,255,0.035) 0 1px, transparent 1px 4px);
      mix-blend-mode: screen;
    }}
    main {{ max-width: 1180px; margin: 0 auto; padding: 2rem; }}
    h1, h2 {{ color: var(--amber); text-transform: uppercase; text-shadow: 0 0 10px rgba(255, 200, 87, 0.65); }}
    h1 {{ font-size: clamp(1.7rem, 4vw, 3rem); margin: 0 0 0.35rem; }}
    h2 {{ font-size: 1.1rem; margin-top: 2rem; }}
    .subhead {{ color: var(--muted); margin: 0 0 1.5rem; }}
    .terminal {{
      border: 2px solid var(--line);
      background: linear-gradient(180deg, rgba(11, 33, 23, 0.94), rgba(4, 13, 9, 0.94));
      box-shadow: 0 0 0 1px rgba(98,255,155,0.16), 0 0 28px rgba(98,255,155,0.12);
      border-radius: 8px;
      padding: 1rem;
    }}
    .summary {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 0.75rem; margin: 1rem 0 1.5rem; }}
    .metric {{
      border: 1px solid var(--line);
      background: rgba(98, 255, 155, 0.055);
      padding: 0.75rem;
      border-radius: 6px;
      min-height: 5rem;
    }}
    .metric strong {{ display: block; color: var(--cyan); font-size: 0.8rem; margin-bottom: 0.35rem; }}
    .metric span {{ color: var(--amber); font-size: 1.35rem; }}
    .terminal-note {{ color: var(--muted); }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 0.75rem; }}
    th, td {{ border-bottom: 1px solid var(--line); padding: 0.55rem; text-align: left; }}
    th {{ color: var(--cyan); font-size: 0.82rem; text-transform: uppercase; }}
    td {{ color: #c9ffd8; }}
    .gel-panel {{ margin: 1rem 0 1.5rem; }}
    .gel-panel svg {{ width: 100%; max-height: 620px; display: block; }}
    .gel-panel figcaption {{ color: var(--muted); font-size: 0.9rem; margin-top: 0.5rem; }}
    .gel-frame {{ fill: #080817; stroke: #395089; stroke-width: 2; }}
    .well-line {{ stroke: #7e8cff; stroke-width: 1.2; opacity: 0.55; }}
    .well {{ fill: #03030c; stroke: #7e8cff; opacity: 0.75; }}
    .marker-band {{ fill: #b6d8ff; filter: url(#glow); opacity: 0.92; }}
    .query-band {{ fill: var(--phosphor); filter: url(#glow); }}
    .reference-band {{ fill: var(--magenta); filter: url(#glow); opacity: 0.88; }}
    .lane-title {{ fill: var(--amber); text-anchor: middle; font: 700 15px "Courier New", monospace; }}
    .gel-legend {{ fill: var(--muted); text-anchor: middle; font: 12px "Courier New", monospace; }}
    .band-label {{ fill: #d6ffe2; font: 10px "Courier New", monospace; opacity: 0.85; }}
    .reference-label {{ fill: #ffd5f3; }}
    .marker-label {{ fill: #b6d8ff; font: 11px "Courier New", monospace; text-anchor: end; }}
  </style>
</head>
<body>
  <main>
    <h1>MLVA Nanopore Report: {_safe(sample_id)}</h1>
    <p class="subhead">VNTR allele terminal // agarose gel comparison // probabilistic fingerprint</p>
    <section class="terminal">
      <div class="summary">
        <div class="metric"><strong>PASS loci</strong><span>{passed}</span></div>
        <div class="metric"><strong>LOW_DEPTH loci</strong><span>{low_depth}</span></div>
        <div class="metric"><strong>DROPOUT loci</strong><span>{dropout}</span></div>
        <div class="metric"><strong>Best reference</strong><span>{_safe(best_match.get('best_profile_id', 'NA') or 'NA')}</span></div>
        <div class="metric"><strong>Novelty</strong><span>{_safe(novelty.get('novelty_score', ''))}</span><br>{_safe(novelty.get('interpretation', ''))}</div>
      </div>
      <h2>Generated Gel</h2>
      {gel}
      <h2>Allele Calls</h2>
      <table>
        <thead><tr><th>Locus</th><th>Call</th><th>Posterior</th><th>Depth</th><th>Status</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
  </main>
</body>
</html>
"""
    (outdir / "report.html").write_text(html)
