from __future__ import annotations

import html
import math
from datetime import datetime, timezone
from pathlib import Path

from .models import Locus


def _safe(value) -> str:
    return html.escape(str(value))


def _metric_card(label: str, value, detail: str = "", tone: str = "") -> str:
    detail_html = f'<small>{_safe(detail)}</small>' if detail else ""
    tone_class = f" {tone}" if tone else ""
    return (
        f'<div class="metric{tone_class}"><strong>{_safe(label)}</strong>'
        f'<span>{_safe(value)}</span>{detail_html}</div>'
    )


def _finding(kind: str, title: str, detail: str) -> str:
    return (
        f'<div class="finding {kind}"><strong>{_safe(title)}</strong>'
        f'<span>{_safe(detail)}</span></div>'
    )


def _closest_reference_detail(row: dict) -> str:
    if row.get("whole_genome_exact_match") == "yes":
        return "Exact whole-genome match"
    if row.get("whole_genome_snps") not in ("", None):
        return (
            f"{row.get('whole_genome_snps')} whole-genome SNPs; "
            f"{row.get('whole_genome_align_fraction_query', '')}% query aligned"
        )
    if row.get("combined_marker_distance") not in ("", None):
        return f"Marker distance {row.get('combined_marker_distance')}"
    return "Closest ranked reference"


def _reference_summary_rows(rows: list[dict]) -> str:
    return "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('rank', ''))}</td>"
        f"<td>{_safe(row.get('reference_id', ''))}</td>"
        f"<td>{_safe(row.get('match_status', ''))}</td>"
        f"<td>{_safe(row.get('combined_marker_distance', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_exact_match', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_snps', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_indel_bases', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_align_fraction_ref', ''))} / "
        f"{_safe(row.get('whole_genome_align_fraction_query', ''))}</td>"
        f"<td>{_safe(row.get('collection_date', ''))}</td>"
        f"<td>{_safe(row.get('location', ''))}</td>"
        "</tr>"
        for row in rows[:10]
    )


def _phylogenetic_warning_html(rows: list[dict]) -> str:
    warned = [
        str(row.get("reference_id", ""))
        for row in rows
        if row.get("ranking_warning") == "EXACT_MATCH_OVERRIDES_PLACEMENT"
    ]
    notices = []
    if warned:
        references = ", ".join(_safe(reference_id) for reference_id in warned)
        notices.append(
            '<div class="warning-banner"><strong>Exact-match placement warning:</strong> '
            "direct sequence identity identified an exact database marker match that "
            "EPA-ng's likelihood-weighted placement distance would have ranked below "
            f"another reference. Identity-aware ranking was applied for: {references}.</div>"
        )
    tie_warning = next(
        (
            str(row.get("ranking_warning", ""))
            for row in rows
            if str(row.get("ranking_warning", "")).startswith("DNADIFF_TIE_BREAK_")
        ),
        "",
    )
    if tie_warning:
        detail = (
            "the reference build has no usable reference-assembly mapping"
            if tie_warning == "DNADIFF_TIE_BREAK_UNAVAILABLE"
            else "dnadiff did not return a result for every tied reference"
        )
        notices.append(
            '<div class="warning-banner"><strong>Whole-genome tie-break warning:</strong> '
            f"{detail}; exact marker matches remain tied where whole-genome evidence "
            "is absent.</div>"
        )
    return "".join(notices)


def _called_count(value) -> int | float | None:
    if value in ("", None):
        return None
    try:
        parsed = float(value)
        return int(parsed) if parsed.is_integer() else parsed
    except (TypeError, ValueError):
        return None


def _amplicon_size(locus: Locus, repeat_count: int | float | None) -> int | float | None:
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


def _locus_repeat_unit_bp(locus: Locus) -> int | None:
    repeat_unit_bp = locus.repeat_unit_length_bp or len(locus.repeat_motif)
    return repeat_unit_bp if repeat_unit_bp > 0 else None


def _gel_svg(
    sample_id: str,
    loci: list[Locus],
    allele_rows: list[dict],
    best_profile: dict | None,
    asv_rows: list[dict],
    closest_reference_bands: list[dict] | None = None,
) -> str:
    closest_reference_bands = closest_reference_bands or []
    allele_by_locus = {row["locus_id"]: row for row in allele_rows}
    query_bands = []
    reference_bands = []
    for locus in loci:
        allele = allele_by_locus.get(locus.locus_id, {})
        called_repeat = _called_count(allele.get("called_repeat_count"))
        query_size = _called_count(allele.get("primary_product_size_bp"))
        if query_size is None:
            query_size = _amplicon_size(locus, called_repeat)
        if query_size:
            support = _called_count(
                allele.get("primary_read_depth", allele.get("read_depth"))
            ) or 0
            frequency = float(allele.get("dominant_variant_fraction") or 0)
            query_bands.append(
                (
                    locus.locus_id,
                    query_size,
                    support,
                    frequency if frequency > 0 else (1.0 if support else 0.0),
                )
            )
        if best_profile and not closest_reference_bands:
            reference_repeat = _called_count(best_profile.get(locus.locus_id))
            repeat_unit_bp = _locus_repeat_unit_bp(locus)
            if (
                query_size
                and called_repeat is not None
                and reference_repeat is not None
                and repeat_unit_bp is not None
            ):
                reference_size = query_size + (
                    (reference_repeat - called_repeat) * repeat_unit_bp
                )
            else:
                reference_size = _amplicon_size(locus, reference_repeat)
            if reference_size:
                reference_bands.append((locus.locus_id, reference_size, 0, 1.0))

    if closest_reference_bands:
        reference_bands = [
            (
                row.get("locus_id", ""),
                size,
                0,
                1.0,
            )
            for row in closest_reference_bands
            if (size := _called_count(row.get("product_size_bp"))) is not None
            and size > 0
        ]

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
                f'<g class="band-hit" tabindex="0"><title>{_safe(locus_id)}: {size} bp; {support} reads; frequency {frequency:.3f}</title>'
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
                f'<g class="band-hit" tabindex="0"><title>{_safe(locus_id)} reference: {size} bp</title>'
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
    reference_name = (
        closest_reference_bands[0].get("reference_id", "closest reference")
        if closest_reference_bands
        else best_profile.get("profile_id", "best reference")
        if best_profile
        else "no reference"
    )
    return f"""
<figure class="gel-panel" aria-label="Generated agarose gel comparison">
  <svg viewBox="0 0 720 500" role="img" aria-labelledby="gel-title gel-desc">
    <title id="gel-title">Generated MLVA agarose gel comparison</title>
    <desc id="gel-desc">Marker, query sample, and closest reference bands estimated from VNTR amplicon sizes.</desc>
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
  <figcaption>Generated gel image: band position shows estimated fragment length; query band brightness and thickness scale with reads supporting that fragment. Horizontal alignment between query and reference bands indicates matching VNTR fragment lengths. Hover or focus a band for locus details.</figcaption>
</figure>
"""


def _assembly_gel_svg(
    sample_id: str,
    call_rows: list[dict],
    best_profile: dict | None = None,
    loci: list[Locus] | None = None,
    closest_reference_bands: list[dict] | None = None,
) -> str:
    closest_reference_bands = closest_reference_bands or []
    bands = []
    call_by_locus = {row.get("locus_id", ""): row for row in call_rows}
    for row in call_rows:
        if row.get("present") != "yes":
            continue
        size = _called_count(row.get("product_size_bp"))
        if not size:
            continue
        support = _called_count(row.get("read_depth")) or 0
        repeat_count = row.get("repeat_count", "")
        bands.append((row.get("locus_id", ""), size, support, repeat_count))

    reference_bands = []
    if closest_reference_bands:
        for row in closest_reference_bands:
            size = _called_count(row.get("product_size_bp"))
            if size is None or size <= 0:
                continue
            reference_bands.append(
                (
                    row.get("locus_id", ""),
                    int(round(size)),
                    row.get("repeat_count", ""),
                )
            )
    elif best_profile and loci:
        for locus in loci:
            reference_repeat = _called_count(best_profile.get(locus.locus_id))
            repeat_unit_bp = _locus_repeat_unit_bp(locus)
            if reference_repeat is None or repeat_unit_bp is None:
                continue
            call = call_by_locus.get(locus.locus_id, {})
            query_size = _called_count(call.get("product_size_bp"))
            query_repeat = _called_count(call.get("repeat_count"))
            if query_size and query_repeat is not None:
                reference_size = query_size + ((reference_repeat - query_repeat) * repeat_unit_bp)
            elif locus.expected_product_size_bp and locus.nominal_repeat_units is not None:
                reference_size = locus.expected_product_size_bp + (
                    (reference_repeat - locus.nominal_repeat_units) * repeat_unit_bp
                )
            else:
                continue
            if reference_size > 0:
                reference_bands.append((locus.locus_id, int(round(reference_size)), reference_repeat))

    if not bands and not reference_bands:
        return '<p class="terminal-note">No assembly VNTR products were available for gel rendering.</p>'

    marker_sizes = [2000, 1500, 1000, 700, 500, 300, 200, 100, 50]
    all_sizes = (
        [size for _name, size, _support, _repeat_count in bands]
        + [size for _name, size, _repeat_count in reference_bands]
        + marker_sizes
    )
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

    max_support = max((support for _locus_id, _size, support, _repeat_count in bands), default=0)
    depth_available = max_support > 0

    def band_intensity(support: int) -> tuple[float, float]:
        if not depth_available:
            return 0.74, 8.0
        scaled = math.sqrt(support / max_support)
        return 0.22 + 0.75 * scaled, 4.0 + 10.0 * scaled

    marker = "\n".join(
        f'<g><rect class="marker-band" x="82" y="{y_for_size(size) - 2.5:.1f}" width="56" height="5" rx="2" />'
        f'<text class="marker-label" x="30" y="{y_for_size(size) + 3:.1f}">{size}</text></g>'
        for size in marker_sizes
        if min_size <= size <= max_size
    )
    band_svg = []
    for locus_id, size, support, repeat_count in bands:
        y = y_for_size(size)
        opacity, height = band_intensity(support)
        label = f"{locus_id} ({repeat_count}U)" if repeat_count != "" else locus_id
        support_label = f"{support} reads" if depth_available else "no depth estimate"
        band_svg.append(
            f'<g class="band-hit" tabindex="0"><title>{_safe(label)}: {size} bp; {support_label}</title>'
            f'<rect class="query-band" x="258" y="{y - height / 2:.1f}" width="92" height="{height:.1f}" '
            f'rx="2" opacity="{opacity:.3f}" />'
            f'<text class="band-label" x="364" y="{y + 3:.1f}">{_safe(label)}</text></g>'
        )

    reference_svg = []
    for locus_id, size, repeat_count in reference_bands:
        y = y_for_size(size)
        label = f"{locus_id} ({repeat_count}U)" if repeat_count != "" else locus_id
        reference_svg.append(
            f'<g class="band-hit" tabindex="0"><title>{_safe(label)} reference: {size} bp</title>'
            f'<rect class="reference-band" x="468" y="{y - 3.5:.1f}" width="92" height="7" rx="2" />'
            f'<text class="band-label reference-label" x="574" y="{y + 3:.1f}">{_safe(label)}</text></g>'
        )

    reference_name = (
        closest_reference_bands[0].get("reference_id", "closest reference")
        if closest_reference_bands
        else best_profile.get("profile_id", "closest reference")
        if best_profile
        else "no reference"
    )
    depth_note = (
        "Band brightness and thickness scale with read depth from FASTQ/BAM support."
        if depth_available
        else "No FASTQ/BAM depth support was provided, so present loci are drawn at a uniform default intensity."
    )
    return f"""
<figure class="gel-panel" aria-label="Generated assembly gel electrophoresis image">
  <svg viewBox="0 0 720 500" role="img" aria-labelledby="assembly-gel-title assembly-gel-desc">
    <title id="assembly-gel-title">Generated MLVA assembly gel electrophoresis image</title>
    <desc id="assembly-gel-desc">Marker, assembly VNTR product bands, and closest reference profile bands estimated from primer product sizes. Band intensity reflects read-depth support when available.</desc>
    <defs>
      <filter id="assembly-glow"><feGaussianBlur stdDeviation="2.4" result="blur"/><feMerge><feMergeNode in="blur"/><feMergeNode in="SourceGraphic"/></feMerge></filter>
      <linearGradient id="assembly-gel-bg" x1="0" x2="1" y1="0" y2="1">
        <stop offset="0%" stop-color="#24103f"/>
        <stop offset="55%" stop-color="#0a0928"/>
        <stop offset="100%" stop-color="#060816"/>
      </linearGradient>
      <pattern id="assembly-scanlines" width="6" height="6" patternUnits="userSpaceOnUse">
        <rect width="6" height="1" fill="#ffffff" opacity="0.045"/>
      </pattern>
    </defs>
    <rect class="gel-frame" x="18" y="18" width="684" height="458" rx="8"/>
    <rect x="42" y="54" width="636" height="388" rx="6" fill="url(#assembly-gel-bg)"/>
    <rect x="42" y="54" width="636" height="388" rx="6" fill="url(#assembly-scanlines)"/>
    <line class="well-line" x1="62" y1="72" x2="658" y2="72"/>
    <rect class="well" x="80" y="58" width="60" height="20" rx="3"/>
    <rect class="well" x="268" y="58" width="72" height="20" rx="3"/>
    <rect class="well" x="478" y="58" width="72" height="20" rx="3"/>
    <text class="lane-title" x="110" y="466">LADDER</text>
    <text class="lane-title" x="304" y="466">{_safe(sample_id)}</text>
    <text class="lane-title" x="514" y="466">{_safe(reference_name)}</text>
    <text class="gel-legend" x="360" y="42">VNTR product bands // intensity = depth support // hover bands for labels</text>
    {marker}
    {"".join(band_svg)}
    {"".join(reference_svg)}
  </svg>
  <figcaption>{depth_note} Closest-reference bands are drawn in magenta when a database placement or profile match is available. Hover or focus a band for locus details.</figcaption>
</figure>
"""


_STATUS_COLORS = {
    "PASS": "#62ff9b",
    "LOW_DEPTH": "#ffc857",
    "AMBIGUOUS": "#ff9f5f",
    "MULTIPLE_VARIANTS": "#ff5fc8",
    "OUT_OF_RANGE": "#ff6b6b",
    "LOCUS_DROPOUT": "#64748b",
}


def _repeat_count_svg(rows: list[dict], assembly: bool = False) -> str:
    if not rows:
        return '<p class="terminal-note">No repeat-count calls were available.</p>'
    normalized = []
    for row in sorted(rows, key=lambda item: str(item.get("locus_id", ""))):
        value_key = "repeat_count" if assembly else "called_repeat_count"
        raw_key = "repeat_count_raw" if assembly else value_key
        count = _called_count(row.get(value_key))
        raw = row.get(raw_key, "")
        status = str(row.get("status" if assembly else "call_status", ""))
        normalized.append((str(row.get("locus_id", "")), count, raw, status))
    max_count = max((count for _locus, count, _raw, _status in normalized if count is not None), default=1)
    row_height = 29
    plot_left = 190
    plot_width = 600
    height = 65 + len(normalized) * row_height
    marks = []
    for index, (locus_id, count, raw, status) in enumerate(normalized):
        y = 42 + index * row_height
        color = _STATUS_COLORS.get(status, "#8dd7aa")
        width = 0 if count is None else (count / max(max_count, 1)) * plot_width
        exact = "NA" if count is None else f"{count} repeats"
        if raw not in ("", None) and str(raw) != str(count):
            exact += f" (raw {raw})"
        marks.append(
            f'<text class="chart-label" x="8" y="{y + 14:.1f}">{_safe(locus_id)}</text>'
            f'<rect class="mapping-track" x="{plot_left}" y="{y:.1f}" width="{plot_width}" height="18" rx="3"/>'
            f'<rect x="{plot_left}" y="{y:.1f}" width="{width:.2f}" height="18" rx="3" fill="{color}" opacity="0.82">'
            f'<title>{_safe(locus_id)}: {_safe(exact)}; {_safe(status)}</title></rect>'
            f'<text class="chart-value" x="{plot_left + plot_width + 18}" y="{y + 13:.1f}">{_safe(exact)} · {_safe(status)}</text>'
        )
    return f"""
<figure class="chart-panel" aria-label="Individual locus repeat counts">
  <svg viewBox="0 0 1080 {height}" role="img">
    <title>Individual locus repeat counts</title>
    <desc>Exact repeat count at every panel locus, shown independently of amplicon SNP bands.</desc>
    {"".join(marks)}
  </svg>
  <figcaption>Bar length represents repeat units, with the exact call (and raw assembly estimate when applicable) printed at right. Color indicates call status.</figcaption>
</figure>
"""


def _locus_confidence_svg(allele_rows: list[dict]) -> str:
    if not allele_rows:
        return '<p class="terminal-note">No locus calls were available.</p>'
    rows = sorted(allele_rows, key=lambda row: str(row.get("locus_id", "")))
    row_height = 30
    width = 1000
    plot_left = 190
    plot_width = 620
    height = 74 + (len(rows) * row_height)
    ticks = []
    for value in (0.0, 0.25, 0.5, 0.75, 1.0):
        x = plot_left + (value * plot_width)
        ticks.append(
            f'<line class="chart-grid" x1="{x:.1f}" y1="38" x2="{x:.1f}" y2="{height - 24}"/>'
            f'<text class="chart-axis" x="{x:.1f}" y="28" text-anchor="middle">{value:.2g}</text>'
        )
    marks = []
    for index, row in enumerate(rows):
        y = 54 + (index * row_height)
        posterior = max(0.0, min(1.0, float(row.get("posterior_probability") or 0)))
        depth = max(
            0,
            int(row.get("primary_read_depth", row.get("read_depth")) or 0),
        )
        radius = min(12.0, 4.5 + math.log10(depth + 1) * 2.2)
        status = str(row.get("call_status") or "")
        color = _STATUS_COLORS.get(status, "#8dd7aa")
        x = plot_left + (posterior * plot_width)
        marks.append(
            f'<text class="chart-label" x="8" y="{y + 4:.1f}">{_safe(row.get("locus_id", ""))}</text>'
            f'<line class="confidence-track" x1="{plot_left}" y1="{y:.1f}" x2="{x:.1f}" y2="{y:.1f}"/>'
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius:.1f}" fill="{color}">'
            f'<title>{_safe(row.get("locus_id", ""))}: posterior {posterior:.3f}; depth {depth}; {status}</title></circle>'
            f'<text class="chart-value" x="{plot_left + plot_width + 22}" y="{y + 4:.1f}">'
            f'{_safe(row.get("called_repeat_count", ""))}U · {_safe(status)}</text>'
        )
    return f"""
<figure class="chart-panel" aria-label="Locus call posterior plot">
  <svg viewBox="0 0 {width} {height}" role="img">
    <title>Locus call confidence</title>
    <desc>Posterior probability by locus. Point size reflects primary-cluster read depth and color reflects call status.</desc>
    {"".join(ticks)}
    {"".join(marks)}
  </svg>
  <figcaption>Farther right is more confident. Point size scales with dominant-cluster read depth; color indicates status.</figcaption>
</figure>
"""


def _variant_mixture_svg(
    mixture_rows: list[dict], allele_rows: list[dict]
) -> str:
    if not mixture_rows:
        return '<p class="terminal-note">No retained variants were available for mixture estimation.</p>'
    by_locus: dict[str, list[dict]] = {}
    for row in mixture_rows:
        by_locus.setdefault(str(row.get("locus_id", "")), []).append(row)
    status_by_locus = {
        str(row.get("locus_id", "")): str(row.get("call_status", ""))
        for row in allele_rows
    }
    palette = ["#62ff9b", "#5fe8ff", "#ff5fc8", "#ffc857", "#9d8cff", "#ff8f70"]
    row_height = 36
    plot_left = 180
    plot_width = 650
    loci = sorted(by_locus)
    height = 70 + (len(loci) * row_height)
    rows_svg = []
    for locus_index, locus_id in enumerate(loci):
        y = 48 + (locus_index * row_height)
        estimates = sorted(
            by_locus[locus_id],
            key=lambda row: float(row.get("estimated_fraction") or 0),
            reverse=True,
        )
        meaningful = [
            row
            for row in estimates
            if str(row.get("meaningful", "")).lower() == "yes"
            and float(row.get("estimated_fraction") or 0) > 0
        ]
        candidate_rows = [
            row
            for row in estimates
            if str(row.get("evidence_class", "")).upper() == "CANDIDATE"
        ]
        trace_rows = [
            row
            for row in estimates
            if row not in meaningful and row not in candidate_rows
        ]
        segments = [
            (
                str(row.get("variant_id", "")),
                float(row.get("estimated_fraction") or 0),
                palette[index % len(palette)],
                f'{row.get("repeat_count", "")}U',
            )
            for index, row in enumerate(meaningful)
        ]
        segments.extend(
            (
                f'{row.get("variant_id", "")} candidate',
                float(row.get("estimated_fraction") or 0),
                "#ffc857",
                f'{row.get("repeat_count", "")}U candidate',
            )
            for row in candidate_rows
        )
        trace_fraction = sum(
            float(row.get("estimated_fraction") or 0) for row in trace_rows
        )
        if trace_fraction > 0:
            segments.append(
                (f"trace ({len(trace_rows)} variants)", trace_fraction, "#64748b", "trace")
            )
        cursor = plot_left
        segment_svg = []
        for label, fraction, color, repeat_label in segments:
            segment_width = max(0.0, fraction * plot_width)
            if segment_width <= 0:
                continue
            segment_svg.append(
                f'<rect x="{cursor:.2f}" y="{y:.1f}" width="{segment_width:.2f}" height="20" '
                f'fill="{color}" rx="2"><title>{_safe(locus_id)} · {_safe(label)} · '
                f'{fraction * 100:.2f}% · {_safe(repeat_label)}</title></rect>'
            )
            if segment_width >= 72:
                segment_svg.append(
                    f'<text class="segment-label" x="{cursor + segment_width / 2:.2f}" '
                    f'y="{y + 14:.1f}" text-anchor="middle">{_safe(label.rsplit("_", 1)[-1])} '
                    f'{fraction * 100:.1f}%</text>'
                )
            cursor += segment_width
        dominant = meaningful[0] if meaningful else estimates[0]
        dominant_fraction = float(dominant.get("estimated_fraction") or 0)
        status = status_by_locus.get(locus_id, "")
        status_color = _STATUS_COLORS.get(status, "#8dd7aa")
        rows_svg.append(
            f'<circle cx="12" cy="{y + 10:.1f}" r="5" fill="{status_color}"/>'
            f'<text class="chart-label" x="24" y="{y + 14:.1f}">{_safe(locus_id)}</text>'
            f'<rect class="mixture-track" x="{plot_left}" y="{y:.1f}" width="{plot_width}" height="20" rx="2"/>'
            f'{"".join(segment_svg)}'
            f'<text class="chart-value" x="{plot_left + plot_width + 18}" y="{y + 14:.1f}">'
            f'{_safe(dominant.get("variant_id", ""))} {dominant_fraction * 100:.1f}%</text>'
        )
    return f"""
<figure class="chart-panel" aria-label="EM-estimated variant abundance plot">
  <svg viewBox="0 0 1000 {height}" role="img">
    <title>Variant mixture abundance</title>
    <desc>Stacked estimated fractions of confirmed, candidate, and trace variants at each locus.</desc>
    <text class="chart-axis" x="{plot_left}" y="26">0%</text>
    <text class="chart-axis" x="{plot_left + plot_width}" y="26" text-anchor="end">100%</text>
    {"".join(rows_svg)}
  </svg>
  <figcaption>Abundance estimates from competitive read-mapping groups. Confirmed variants are colored separately, candidates are amber, and trace components are combined in gray.</figcaption>
</figure>
"""


def _mapping_coverage_svg(mapping_rows: list[dict]) -> str:
    if not mapping_rows:
        return '<p class="terminal-note">Representative mapping was disabled or produced no locus summaries.</p>'
    rows = sorted(mapping_rows, key=lambda row: str(row.get("locus_id", "")))
    row_height = 31
    plot_left = 180
    plot_width = 620
    height = 66 + (len(rows) * row_height)
    bars = []
    for index, row in enumerate(rows):
        y = 42 + (index * row_height)
        coverage = max(0.0, min(100.0, float(row.get("coverage_percent") or 0)))
        bar_width = plot_width * coverage / 100.0
        snps = int(row.get("snp_count") or 0)
        bars.append(
            f'<text class="chart-label" x="8" y="{y + 14:.1f}">{_safe(row.get("locus_id", ""))}</text>'
            f'<rect class="mapping-track" x="{plot_left}" y="{y:.1f}" width="{plot_width}" height="20" rx="3"/>'
            f'<rect class="mapping-bar" x="{plot_left}" y="{y:.1f}" width="{bar_width:.2f}" height="20" rx="3">'
            f'<title>{coverage:.1f}% covered; mean depth {row.get("mean_depth", 0)}; {snps} SNPs</title></rect>'
            f'<text class="chart-value" x="{plot_left + plot_width + 18}" y="{y + 14:.1f}">'
            f'{coverage:.1f}% · {snps} SNP</text>'
        )
    return f"""
<figure class="chart-panel" aria-label="Representative mapping coverage plot">
  <svg viewBox="0 0 1000 {height}" role="img">
    <title>Representative mapping coverage</title>
    <desc>Percent of each sample-derived representative covered by quality-filtered aligned bases.</desc>
    {"".join(bars)}
  </svg>
  <figcaption>Coverage is relative to the sample-derived dominant amplicon. Hover bars for depth and SNP counts.</figcaption>
</figure>
"""


def write_report(
    outdir: str | Path,
    sample_id: str,
    allele_rows: list[dict],
    loci: list[Locus] | None = None,
    match_rows: list[dict] | None = None,
    profiles: list[dict] | None = None,
    asv_rows: list[dict] | None = None,
    mapping_rows: list[dict] | None = None,
    snp_rows: list[dict] | None = None,
    mixture_rows: list[dict] | None = None,
    phylogenetic_rows: list[dict] | None = None,
    closest_reference_bands: list[dict] | None = None,
    presence_rows: list[dict] | None = None,
    local_assembly_rows: list[dict] | None = None,
    taxon_screen_summary: dict | None = None,
) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    loci = loci or []
    match_rows = match_rows or []
    profiles = profiles or []
    asv_rows = asv_rows or []
    mapping_rows = mapping_rows or []
    snp_rows = snp_rows or []
    mixture_rows = mixture_rows or []
    phylogenetic_rows = phylogenetic_rows or []
    closest_reference_bands = closest_reference_bands or []
    presence_rows = presence_rows or []
    local_assembly_rows = local_assembly_rows or []
    taxon_screen_summary = taxon_screen_summary or {}
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    passed = sum(1 for row in allele_rows if row.get("call_status") == "PASS")
    low_depth = sum(1 for row in allele_rows if row.get("call_status") == "LOW_DEPTH")
    dropout = sum(1 for row in allele_rows if row.get("call_status") == "LOCUS_DROPOUT")
    multiple = sum(
        1 for row in allele_rows if row.get("call_status") == "MULTIPLE_VARIANTS"
    )
    best_match = match_rows[0] if match_rows else {}
    best_profile = _best_profile(match_rows, profiles)
    phylogenetic_best = phylogenetic_rows[0] if phylogenetic_rows else {}
    total_loci = len(allele_rows)
    flagged = low_depth + dropout + multiple + sum(
        row.get("call_status") in {"AMBIGUOUS", "OUT_OF_RANGE"}
        for row in allele_rows
    )
    summary_cards = [
        _metric_card("Confident calls", f"{passed}/{total_loci}", "PASS loci"),
        _metric_card(
            "Loci needing review",
            flagged,
            "Low depth, mixed, ambiguous, out-of-range, or missing",
            "warn" if flagged else "good",
        ),
    ]
    if taxon_screen_summary:
        screen_input = int(taxon_screen_summary.get("seqs_in", 0))
        screen_retained = int(taxon_screen_summary.get("seqs_out", 0))
        retained_fraction = (
            screen_retained / screen_input if screen_input else 0.0
        )
        summary_cards.append(
            _metric_card(
                "Target-taxon screen",
                f"{screen_retained:,}/{screen_input:,}",
                f"{retained_fraction:.2%} of reads retained by Deacon",
                "good" if screen_retained else "warn",
            )
        )
    if presence_rows:
        detected = sum(
            row.get("presence_status") != "NO_EVIDENCE"
            for row in presence_rows
        )
        genotyped = sum(
            row.get("presence_status")
            in {"PRESENT_GENOTYPED", "PRESENT_PROVISIONAL"}
            for row in presence_rows
        )
        summary_cards.append(
            _metric_card(
                "Recruited loci",
                f"{detected}/{len(presence_rows)}",
                f"{genotyped} with repeat-informative evidence",
            )
        )
    if local_assembly_rows:
        poa_passed = sum(
            row.get("pcr_status") == "PASS" for row in local_assembly_rows
        )
        summary_cards.append(
            _metric_card(
                "POA assembly calls",
                f"{poa_passed}/{len(local_assembly_rows)}",
                "Dominant locus consensuses resolved by assembly PCR",
                "good" if poa_passed == len(local_assembly_rows) else "warn",
            )
        )
    if mapping_rows:
        summary_cards.extend(
            [
                _metric_card(
                    "Assigned locus reads",
                    f"{sum(int(row.get('mapped_reads') or 0) for row in mapping_rows)}/"
                    f"{sum(int(row.get('total_reads') or 0) for row in mapping_rows)}",
                ),
                _metric_card("Amplicon SNP observations", len(snp_rows)),
            ]
        )
    if best_match.get("best_profile_id"):
        summary_cards.append(
            _metric_card(
                "Closest MLVA profile",
                best_match["best_profile_id"],
                f"Distance {best_match.get('distance', '')}; confidence {best_match.get('confidence', '')}",
            )
        )
    if phylogenetic_best.get("reference_id"):
        summary_cards.append(
            _metric_card(
                "Closest reference genome",
                phylogenetic_best["reference_id"],
                _closest_reference_detail(phylogenetic_best),
                "good" if phylogenetic_best.get("whole_genome_exact_match") == "yes" else "",
            )
        )
    findings = []
    if taxon_screen_summary and not int(taxon_screen_summary.get("seqs_out", 0)):
        findings.append(
            _finding(
                "warn",
                "Target taxon not detected",
                "The Deacon screen retained no reads; downstream locus calls "
                "represent target-taxon dropout rather than an unfiltered sample.",
            )
        )
    for status, title in (
        ("LOCUS_DROPOUT", "Missing loci"),
        ("LOW_DEPTH", "Low-depth loci"),
        ("MULTIPLE_VARIANTS", "Mixed loci"),
        ("AMBIGUOUS", "Ambiguous loci"),
        ("OUT_OF_RANGE", "Out-of-range loci"),
    ):
        affected = [str(row.get("locus_id", "")) for row in allele_rows if row.get("call_status") == status]
        if affected:
            findings.append(_finding("warn", title, ", ".join(affected)))
    poa_fallbacks = [
        str(row.get("locus_id", ""))
        for row in local_assembly_rows
        if row.get("pcr_status") != "PASS"
    ]
    if poa_fallbacks:
        findings.append(
            _finding(
                "warn",
                "Local assembly fallbacks",
                ", ".join(poa_fallbacks),
            )
        )
    if not findings:
        findings.append(_finding("good", "Panel quality", "No locus-level review flags were detected."))
    if phylogenetic_best.get("reference_id"):
        findings.append(
            _finding(
                "info",
                "Reference interpretation",
                f"{phylogenetic_best['reference_id']}: {_closest_reference_detail(phylogenetic_best)}",
            )
        )
    summary_html = "".join(summary_cards)
    findings_html = "".join(findings)
    gel = _gel_svg(
        sample_id,
        loci,
        allele_rows,
        best_profile,
        asv_rows,
        closest_reference_bands,
    )
    repeat_count_plot = _repeat_count_svg(allele_rows)
    confidence_plot = _locus_confidence_svg(allele_rows)
    mixture_plot = _variant_mixture_svg(mixture_rows, allele_rows)
    mapping_plot = _mapping_coverage_svg(mapping_rows)
    rows = "\n".join(
        f"<tr><td>{_safe(row['locus_id'])}</td><td>{_safe(row['called_repeat_count'])}</td>"
        f"<td>{_safe(row.get('primary_product_size_bp', ''))}</td>"
        f"<td>{_safe(row.get('primary_repeat_count_raw', ''))}</td>"
        f"<td>{_safe(row.get('primary_measurement_source', ''))}</td>"
        f"<td>{_safe(row['posterior_probability'])}</td>"
        f"<td>{_safe(row.get('primary_read_depth', ''))}/{_safe(row['read_depth'])}</td>"
        f"<td>{_safe(row.get('num_meaningful_variants', ''))}</td>"
        f"<td>{_safe(row.get('dominant_variant_fraction', ''))}</td>"
        f"<td>{_safe(row['call_status'])}</td></tr>"
        for row in allele_rows
    )
    local_assembly_table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('locus_id', ''))}</td>"
        f"<td>{_safe(row.get('dominant_variant_id', ''))}</td>"
        f"<td>{_safe(row.get('input_reads', ''))}</td>"
        f"<td>{_safe(row.get('unique_sequences', ''))}</td>"
        f"<td>{_safe(row.get('observed_min_product_bp', ''))} / "
        f"{_safe(row.get('observed_modal_product_bp', ''))} / "
        f"{_safe(row.get('observed_max_product_bp', ''))}</td>"
        f"<td>{_safe(row.get('poa_consensus_bp', ''))}</td>"
        f"<td>{_safe(row.get('pcr_product_size_bp', ''))}</td>"
        f"<td>{_safe(row.get('raw_repeat_count', ''))}</td>"
        f"<td>{_safe(row.get('called_repeat_count', ''))}</td>"
        f"<td>{_safe(row.get('measurement_source', ''))}</td>"
        f"<td><span class=\"status-pill {'status-good' if row.get('pcr_status') == 'PASS' else 'status-warn'}\">"
        f"{_safe(row.get('pcr_status', ''))}</span></td>"
        "</tr>"
        for row in local_assembly_rows
    )
    local_assembly_section = ""
    if local_assembly_rows:
        local_assembly_section = f"""
      <section class="report-section">
        <h2>FASTQ Local Assembly Concordance</h2>
        <p class="section-intro">SPOARS assembles the dominant mapping-derived product group, then the standard assembly PCR caller measures the consensus. Compare the PCR product and final repeat columns directly with the corresponding assembly report. Min / mode / max shows the uncorrected read-product length distribution.</p>
        <div class="table-scroll"><table>
          <thead><tr><th>Locus</th><th>Dominant mapped group</th><th>Reads</th><th>Unique products</th><th>Raw bp min / mode / max</th><th>POA bp</th><th>Assembly PCR bp</th><th>Raw repeats</th><th>Final repeats</th><th>Call source</th><th>POA status</th></tr></thead>
          <tbody>{local_assembly_table_rows}</tbody>
        </table></div>
      </section>
"""
    presence_table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('locus_id', ''))}</td>"
        f"<td>{_safe(row.get('presence_status', ''))}</td>"
        f"<td>{_safe(row.get('mapped_reads', ''))}</td>"
        f"<td>{_safe(row.get('full_product_reads', ''))}</td>"
        f"<td>{_safe(row.get('genotype_informative_reads', ''))}</td>"
        f"<td>{_safe(row.get('candidate_alleles', ''))}</td>"
        f"<td>{_safe(row.get('reference_source', ''))}</td>"
        "</tr>"
        for row in presence_rows
    )
    presence_detail_section = ""
    if presence_rows:
        presence_detail_section = f"""
      <details>
        <summary>Locus recruitment and presence evidence</summary>
        <div class="table-scroll"><table>
          <thead><tr><th>Locus</th><th>Presence</th><th>Mapped reads</th><th>Full products</th><th>Repeat-informative</th><th>Candidate alleles</th><th>Reference source</th></tr></thead>
          <tbody>{presence_table_rows}</tbody>
        </table></div>
      </details>
"""
    mixture_table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('locus_id', ''))}</td>"
        f"<td>{_safe(row.get('variant_id', ''))}</td>"
        f"<td>{_safe(row.get('repeat_count', ''))}</td>"
        f"<td>{_safe(row.get('observed_reads', ''))}</td>"
        f"<td>{float(row.get('observed_fraction') or 0) * 100:.3f}%</td>"
        f"<td>{_safe(row.get('estimated_reads', ''))}</td>"
        f"<td>{float(row.get('estimated_fraction') or 0) * 100:.3f}%</td>"
        f"<td>{_safe(row.get('abundance_class', ''))}</td>"
        f"<td>{_safe(row.get('evidence_class', ''))}</td>"
        f"<td>{_safe(row.get('meaningful', ''))}</td>"
        "</tr>"
        for row in mixture_rows
    )
    if not mixture_table_rows:
        mixture_table_rows = (
            '<tr><td colspan="10">No retained variants were available for mixture estimation.</td></tr>'
        )
    mapping_table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('locus_id', ''))}</td>"
        f"<td>{_safe(row.get('reference_variant_id', ''))}</td>"
        f"<td>{_safe(row.get('mapped_reads', ''))}/{_safe(row.get('total_reads', ''))}</td>"
        f"<td>{_safe(row.get('mapping_rate', ''))}</td>"
        f"<td>{_safe(row.get('mean_depth', ''))}</td>"
        f"<td>{_safe(row.get('coverage_percent', ''))}</td>"
        f"<td>{_safe(row.get('snp_count', ''))}</td>"
        "</tr>"
        for row in mapping_rows
    )
    if not mapping_table_rows:
        mapping_table_rows = (
            '<tr><td colspan="7">Locus representative mapping was disabled or no retained '
            "representatives were available.</td></tr>"
        )
    snp_table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('locus_id', ''))}</td>"
        f"<td>{_safe(row.get('reference_variant_id', ''))}</td>"
        f"<td>{_safe(row.get('position', ''))}</td>"
        f"<td>{_safe(row.get('reference_base', ''))}&gt;{_safe(row.get('alternate_base', ''))}</td>"
        f"<td>{_safe(row.get('alternate_depth', ''))}/{_safe(row.get('depth', ''))}</td>"
        f"<td>{_safe(row.get('alternate_frequency', ''))}</td>"
        f"<td>{_safe(row.get('mean_alternate_base_quality', ''))}</td>"
        "</tr>"
        for row in snp_rows
    )
    if not snp_table_rows:
        snp_table_rows = (
            '<tr><td colspan="7">No SNPs passed the configured depth, base-quality, '
            "support, and allele-frequency thresholds.</td></tr>"
        )
    profile_table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('best_profile_id', ''))}</td>"
        f"<td>{_safe(row.get('strain_id', ''))}</td>"
        f"<td>{_safe(row.get('distance', ''))}</td>"
        f"<td>{_safe(row.get('matched_loci', ''))}</td>"
        f"<td>{_safe(row.get('mismatched_loci', ''))}</td>"
        f"<td>{_safe(row.get('confidence', ''))}</td>"
        "</tr>"
        for row in match_rows[:10]
    )
    profile_section = ""
    if profiles:
        if not profile_table_rows:
            profile_table_rows = '<tr><td colspan="6">No comparable profile rows were available.</td></tr>'
        profile_section = f"""
      <section class="report-section">
        <h2>Closest MLVA Profiles</h2>
        <p class="section-intro">Direct repeat-count comparison against the supplied profile collection.</p>
        <div class="table-scroll"><table>
          <thead><tr><th>Profile</th><th>Strain</th><th>Distance</th><th>Matched loci</th><th>Mismatched loci</th><th>Confidence</th></tr></thead>
          <tbody>{profile_table_rows}</tbody>
        </table></div>
      </section>
"""
    reference_summary_rows = _reference_summary_rows(phylogenetic_rows)
    phylogenetic_table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('rank', ''))}</td>"
        f"<td>{_safe(row.get('reference_id', ''))}</td>"
        f"<td>{_safe(row.get('combined_marker_distance', ''))}</td>"
        f"<td>{_safe(row.get('total_normalized_snp_distance', ''))}</td>"
        f"<td>{_safe(row.get('total_placement_normalized_snp_distance', ''))}</td>"
        f"<td>{_safe(row.get('total_normalized_direct_snp_distance', ''))}</td>"
        f"<td>{_safe(row.get('total_normalized_repeat_distance', ''))}</td>"
        f"<td>{_safe(row.get('compared_loci', ''))}</td>"
        f"<td>{_safe(row.get('exact_marker_loci', ''))}</td>"
        f"<td>{_safe(row.get('match_status', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_exact_match', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_snps', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_indel_bases', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_align_fraction_ref', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_align_fraction_query', ''))}</td>"
        f"<td>{_safe(row.get('tie_break_status', ''))}</td>"
        f"<td>{_safe(row.get('distance_gap_to_next', ''))}</td>"
        f"<td>{_safe(row.get('collection_date', ''))}</td>"
        f"<td>{_safe(row.get('location', ''))}</td>"
        "</tr>"
        for row in phylogenetic_rows[:10]
    )
    phylogenetic_section = ""
    if phylogenetic_rows:
        phylogenetic_warning = _phylogenetic_warning_html(phylogenetic_rows)
        phylogenetic_section = f"""
      <section class="report-section">
        <h2>Closest Reference Genomes</h2>
        {phylogenetic_warning}
        <p class="section-intro">Marker similarity identifies the candidate group. Exact assembly ties are resolved by canonical genome identity and MUMmer whole-genome SNP comparison.</p>
        <div class="table-scroll"><table class="primary-table">
          <thead><tr><th>Rank</th><th>Reference</th><th>Marker status</th><th>Marker distance</th><th>Exact genome</th><th>WG SNPs</th><th>Indel bases</th><th>Aligned % ref/query</th><th>Date</th><th>Location</th></tr></thead>
          <tbody>{reference_summary_rows}</tbody>
        </table></div>
        <details>
          <summary>Technical marker-distance components</summary>
          <p class="terminal-note">EPA/tree, direct aligned-sequence, repeat, and tie-break components used to construct the ranking.</p>
          <div class="table-scroll"><table>
            <thead><tr><th>Rank</th><th>Reference</th><th>Combined distance</th><th>Hybrid SNP</th><th>EPA/tree SNP</th><th>Direct SNP</th><th>Normalized repeat</th><th>Compared loci</th><th>Exact marker loci</th><th>Match status</th><th>Exact genome</th><th>WG SNPs</th><th>Indel bases</th><th>Ref AF</th><th>Query AF</th><th>Tie break</th><th>Gap to next</th><th>Date</th><th>Location</th></tr></thead>
            <tbody>{phylogenetic_table_rows}</tbody>
          </table></div>
        </details>
      </section>
"""
    mixture_overview_section = ""
    mixture_detail_section = ""
    if mixture_rows:
        mixture_overview_section = f"""
      <h2>Variant Mixture Abundance</h2>
      <p class="section-intro">Estimated within-locus variant fractions distinguish dominant alleles from meaningful mixtures and trace evidence.</p>
      <div class="chart-scroll">{mixture_plot}</div>
"""
        mixture_detail_section = f"""
      <details>
        <summary>Variant mixture details</summary>
        <div class="table-scroll"><table>
          <thead><tr><th>Locus</th><th>Variant</th><th>Repeat</th><th>Observed reads</th><th>Observed %</th><th>EM reads</th><th>EM %</th><th>Abundance</th><th>Evidence tier</th><th>Meaningful</th></tr></thead>
          <tbody>{mixture_table_rows}</tbody>
        </table></div>
      </details>
"""
    mapping_overview_section = ""
    mapping_detail_section = ""
    if mapping_rows:
        mapping_overview_section = f"""
      <h2>Representative Mapping Coverage</h2>
      <p class="section-intro">Coverage and SNP observations are relative to each sample-derived dominant amplicon, not chromosomal coordinates.</p>
      <div class="chart-scroll">{mapping_plot}</div>
"""
        mapping_detail_section = f"""
      <details>
        <summary>Representative mapping details</summary>
        <div class="table-scroll"><table>
          <thead><tr><th>Locus</th><th>POA reference</th><th>Mapped</th><th>Rate</th><th>Mean depth</th><th>Covered %</th><th>SNPs</th></tr></thead>
          <tbody>{mapping_table_rows}</tbody>
        </table></div>
      </details>
      <details>
        <summary>SNP evidence details ({len(snp_rows)} rows)</summary>
        <div class="table-scroll"><table>
          <thead><tr><th>Locus</th><th>POA reference</th><th>Position</th><th>Change</th><th>Alt/depth</th><th>Frequency</th><th>Mean alt Q</th></tr></thead>
          <tbody>{snp_table_rows}</tbody>
        </table></div>
      </details>
"""
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta http-equiv="Cache-Control" content="no-cache, no-store, must-revalidate">
  <meta http-equiv="Pragma" content="no-cache">
  <meta http-equiv="Expires" content="0">
  <title>MLVAMaps Report - {sample_id}</title>
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
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
      letter-spacing: 0;
    }}
    main {{ max-width: 1320px; margin: 0 auto; padding: 2rem; }}
    h1, h2 {{ color: var(--amber); }}
    h1 {{ font-size: clamp(1.7rem, 4vw, 3rem); margin: 0 0 0.35rem; }}
    h2 {{ font-size: 1.1rem; margin-top: 2rem; }}
    .subhead {{ color: var(--muted); margin: 0 0 1.5rem; }}
    .generated-at {{ color: var(--muted); font-size: 0.78rem; margin: -1rem 0 1.5rem; }}
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
    .metric small {{ display: block; color: var(--muted); line-height: 1.35; margin-top: 0.35rem; }}
    .metric.good {{ border-color: rgba(98,255,155,0.65); }}
    .metric.warn {{ border-color: rgba(255,200,87,0.75); }}
    .findings {{ display: grid; gap: 0.55rem; margin: 0 0 1.5rem; }}
    .finding {{ display: grid; grid-template-columns: minmax(150px, 220px) 1fr; gap: 0.8rem; padding: 0.75rem 0.9rem; border-left: 4px solid var(--cyan); background: rgba(95,232,255,0.06); border-radius: 4px; }}
    .finding strong {{ color: var(--cyan); }}
    .finding span {{ color: #d8fbe2; }}
    .finding.warn {{ border-left-color: var(--amber); background: rgba(255,200,87,0.07); }}
    .finding.warn strong {{ color: var(--amber); }}
    .finding.good {{ border-left-color: var(--phosphor); }}
    .finding.good strong {{ color: var(--phosphor); }}
    .report-section {{ margin-top: 2.2rem; padding-top: 0.2rem; }}
    .section-intro {{ color: var(--muted); max-width: 72rem; line-height: 1.5; }}
    .terminal-note {{ color: var(--muted); }}
    .warning-banner {{ margin: 0.75rem 0; padding: 0.8rem 1rem; color: #1b1200; background: var(--amber); border: 2px solid #ff8c42; border-radius: 5px; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 0.75rem; }}
    th, td {{ border-bottom: 1px solid var(--line); padding: 0.55rem; text-align: left; }}
    th {{ color: var(--cyan); font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.04em; }}
    td {{ color: #d8fbe2; }}
    .primary-table tbody tr:first-child {{ background: rgba(98,255,155,0.08); }}
    .status-pill {{ display: inline-block; border: 1px solid; border-radius: 999px; padding: 0.16rem 0.48rem; white-space: nowrap; font-size: 0.78rem; }}
    .status-good {{ color: var(--phosphor); border-color: var(--phosphor); background: rgba(98,255,155,0.08); }}
    .status-warn {{ color: var(--amber); border-color: var(--amber); background: rgba(255,200,87,0.08); }}
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
    .band-hit .band-label {{ fill: #d6ffe2; font: 10px "Courier New", monospace; opacity: 0; pointer-events: none; transition: opacity 120ms ease; }}
    .band-hit:hover .band-label, .band-hit:focus .band-label {{ opacity: 0.95; }}
    .band-hit:focus {{ outline: none; }}
    .reference-label {{ fill: #ffd5f3; }}
    .marker-label {{ fill: #b6d8ff; font: 11px "Courier New", monospace; text-anchor: end; }}
    .chart-grid {{ stroke: rgba(141, 215, 170, 0.18); stroke-width: 1; }}
    .chart-axis {{ fill: var(--muted); font: 12px "Courier New", monospace; }}
    .chart-label {{ fill: #d6ffe2; font: 12px "Courier New", monospace; }}
    .chart-value {{ fill: var(--muted); font: 11px "Courier New", monospace; }}
    .segment-label {{ fill: #03130b; font: 700 10px "Courier New", monospace; pointer-events: none; }}
    .confidence-track {{ stroke: rgba(95, 232, 255, 0.42); stroke-width: 3; }}
    .mixture-track, .mapping-track {{ fill: rgba(141, 215, 170, 0.1); stroke: var(--line); stroke-width: 1; }}
    .mapping-bar {{ fill: var(--cyan); opacity: 0.82; }}
    .chart-panel {{ margin: 0.75rem 0 1.5rem; border: 1px solid var(--line); border-radius: 7px; padding: 0.6rem; background: rgba(3, 12, 8, 0.58); }}
    .chart-panel svg {{ width: 100%; display: block; min-width: 720px; }}
    .chart-panel figcaption {{ color: var(--muted); font-size: 0.88rem; margin-top: 0.45rem; }}
    .chart-scroll {{ overflow-x: auto; }}
    details {{ border: 1px solid var(--line); border-radius: 6px; margin: 0.7rem 0; padding: 0.65rem 0.8rem; background: rgba(98, 255, 155, 0.035); }}
    summary {{ color: var(--cyan); cursor: pointer; font-weight: 700; }}
    .table-scroll {{ overflow-x: auto; }}
    @media (max-width: 700px) {{
      main {{ padding: 1rem; }}
      .finding {{ grid-template-columns: 1fr; gap: 0.2rem; }}
    }}
  </style>
</head>
<body>
  <main>
    <h1>MLVAMaps Report: {_safe(sample_id)}</h1>
    <p class="subhead">VNTR calls, sample quality, mixture evidence, and closest-reference interpretation.</p>
    <p class="generated-at">Generated {_safe(generated_at)} UTC</p>
    <section class="terminal">
      <h2>Sample Overview</h2>
      <div class="summary">
        {summary_html}
      </div>
      <div class="findings">{findings_html}</div>
      <h2>Individual Locus Repeat Counts</h2>
      <div class="chart-scroll">{repeat_count_plot}</div>
      {local_assembly_section}
      <h2>Locus Confidence</h2>
      <div class="chart-scroll">{confidence_plot}</div>
      {mixture_overview_section}
      {mapping_overview_section}
      <h2>Generated Gel</h2>
      {gel}
      {profile_section}
      {phylogenetic_section}
      <h2>Detailed Evidence</h2>
      <p class="terminal-note">The plots above are the primary interpretation view. Expand a section below when exact values are needed.</p>
      <details>
        <summary>Allele call details</summary>
        <div class="table-scroll"><table>
          <thead><tr><th>Locus</th><th>Primary call</th><th>Product bp</th><th>Raw repeat</th><th>Measurement source</th><th>Confidence</th><th>Primary/total depth</th><th>Meaningful variants</th><th>Dominant fraction</th><th>Status</th></tr></thead>
          <tbody>{rows}</tbody>
        </table></div>
      </details>
      {presence_detail_section}
      {mixture_detail_section}
      {mapping_detail_section}
    </section>
  </main>
</body>
</html>
"""
    (outdir / "report.html").write_text(html)


def write_assembly_report(
    outdir: str | Path,
    sample_id: str,
    call_rows: list[dict],
    product_rows: list[dict],
    match_rows: list[dict] | None = None,
    profiles: list[dict] | None = None,
    loci: list[Locus] | None = None,
    phylogenetic_rows: list[dict] | None = None,
    closest_reference_bands: list[dict] | None = None,
) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    match_rows = match_rows or []
    profiles = profiles or []
    phylogenetic_rows = phylogenetic_rows or []
    closest_reference_bands = closest_reference_bands or []
    present = sum(1 for row in call_rows if row.get("present") == "yes")
    not_found = sum(1 for row in call_rows if row.get("present") != "yes")
    with_depth = sum(1 for row in call_rows if _called_count(row.get("read_depth")))
    best_match = match_rows[0] if match_rows else {}
    best_profile = _best_profile(match_rows, profiles)
    phylogenetic_best = phylogenetic_rows[0] if phylogenetic_rows else {}
    gel = _assembly_gel_svg(
        sample_id,
        call_rows,
        best_profile,
        loci,
        closest_reference_bands,
    )
    repeat_count_plot = _repeat_count_svg(call_rows, assembly=True)
    summary_cards = [
        _metric_card(
            "Panel recovered",
            f"{present}/{len(call_rows)}",
            "Loci with a selected primer product",
            "good" if not not_found else "warn",
        ),
        _metric_card("Primer products", len(product_rows)),
    ]
    if any(row.get("read_depth") not in ("", None) for row in call_rows):
        summary_cards.append(
            _metric_card("Read-supported loci", f"{with_depth}/{present}")
        )
    if best_match.get("best_profile_id"):
        summary_cards.append(
            _metric_card(
                "Closest MLVA profile",
                best_match["best_profile_id"],
                f"Distance {best_match.get('distance', '')}; confidence {best_match.get('confidence', '')}",
            )
        )
    if phylogenetic_best.get("reference_id"):
        summary_cards.append(
            _metric_card(
                "Closest reference genome",
                phylogenetic_best["reference_id"],
                _closest_reference_detail(phylogenetic_best),
                "good" if phylogenetic_best.get("whole_genome_exact_match") == "yes" else "",
            )
        )
    findings = []
    missing_loci = [
        str(row.get("locus_id", ""))
        for row in call_rows
        if row.get("present") != "yes"
    ]
    if missing_loci:
        findings.append(_finding("warn", "Loci not recovered", ", ".join(missing_loci)))
    review_rows = [
        row
        for row in call_rows
        if row.get("present") == "yes"
        and row.get("status") not in {"PASS", "PRESENT", "INCLUDED", ""}
    ]
    if review_rows:
        findings.append(
            _finding(
                "warn",
                "Calls needing review",
                ", ".join(
                    f"{row.get('locus_id')} ({row.get('status')})" for row in review_rows
                ),
            )
        )
    if not findings:
        findings.append(_finding("good", "Panel quality", "All configured loci were recovered without call flags."))
    if phylogenetic_best.get("reference_id"):
        findings.append(
            _finding(
                "info",
                "Reference interpretation",
                f"{phylogenetic_best['reference_id']}: {_closest_reference_detail(phylogenetic_best)}",
            )
        )
    summary_html = "".join(summary_cards)
    findings_html = "".join(findings)
    table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('locus_id', ''))}</td>"
        f"<td>{_safe(row.get('present', ''))}</td>"
        f"<td>{_safe(row.get('repeat_count', ''))}</td>"
        f"<td>{_safe(row.get('repeat_count_raw', ''))}</td>"
        f"<td>{_safe(row.get('product_size_bp', ''))}</td>"
        f"<td>{_safe(row.get('read_depth', ''))}</td>"
        f"<td>{_safe(row.get('mean_coverage', ''))}</td>"
        f"<td>{_safe(row.get('allele_confidence', ''))}</td>"
        f"<td>{_safe(row.get('second_best_repeat_count', ''))}</td>"
        f"<td>{_safe(row.get('second_best_probability', ''))}</td>"
        f"<td>{_safe(row.get('inference_method', ''))}</td>"
        f"<td>{_safe(row.get('status', ''))}</td>"
        f"<td>{_safe(row.get('evidence', ''))}</td>"
        "</tr>"
        for row in call_rows
    )
    product_table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('locus_id', ''))}</td>"
        f"<td>{_safe(row.get('contig', ''))}</td>"
        f"<td>{_safe(row.get('contig_start', ''))}-{_safe(row.get('contig_end', ''))}</td>"
        f"<td>{_safe(row.get('orientation', ''))}</td>"
        f"<td>{_safe(row.get('product_size_bp', ''))}</td>"
        f"<td>{_safe(row.get('forward_mismatches', ''))}</td>"
        f"<td>{_safe(row.get('reverse_mismatches', ''))}</td>"
        "</tr>"
        for row in product_rows
    )
    if not product_table_rows:
        product_table_rows = '<tr><td colspan="7">No primer products found in the assembly.</td></tr>'
    match_table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('best_profile_id', ''))}</td>"
        f"<td>{_safe(row.get('strain_id', ''))}</td>"
        f"<td>{_safe(row.get('distance', ''))}</td>"
        f"<td>{_safe(row.get('matched_loci', ''))}</td>"
        f"<td>{_safe(row.get('mismatched_loci', ''))}</td>"
        f"<td>{_safe(row.get('confidence', ''))}</td>"
        "</tr>"
        for row in match_rows[:10]
    )
    profile_section = ""
    if profiles:
        if not match_table_rows:
            match_table_rows = '<tr><td colspan="6">No comparable profile rows were available.</td></tr>'
        profile_section = f"""
      <h2>Closest MLVA Profiles</h2>
      <table>
        <thead><tr><th>Profile</th><th>Strain</th><th>Distance</th><th>Matched loci</th><th>Mismatched loci</th><th>Confidence</th></tr></thead>
        <tbody>{match_table_rows}</tbody>
      </table>
"""
    reference_summary_rows = _reference_summary_rows(phylogenetic_rows)
    phylogenetic_table_rows = "\n".join(
        "<tr>"
        f"<td>{_safe(row.get('rank', ''))}</td>"
        f"<td>{_safe(row.get('reference_id', ''))}</td>"
        f"<td>{_safe(row.get('combined_marker_distance', ''))}</td>"
        f"<td>{_safe(row.get('total_normalized_snp_distance', ''))}</td>"
        f"<td>{_safe(row.get('total_placement_normalized_snp_distance', ''))}</td>"
        f"<td>{_safe(row.get('total_normalized_direct_snp_distance', ''))}</td>"
        f"<td>{_safe(row.get('total_normalized_repeat_distance', ''))}</td>"
        f"<td>{_safe(row.get('compared_loci', ''))}</td>"
        f"<td>{_safe(row.get('exact_marker_loci', ''))}</td>"
        f"<td>{_safe(row.get('match_status', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_exact_match', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_snps', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_indel_bases', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_align_fraction_ref', ''))}</td>"
        f"<td>{_safe(row.get('whole_genome_align_fraction_query', ''))}</td>"
        f"<td>{_safe(row.get('tie_break_status', ''))}</td>"
        f"<td>{_safe(row.get('distance_gap_to_next', ''))}</td>"
        f"<td>{_safe(row.get('collection_date', ''))}</td>"
        f"<td>{_safe(row.get('location', ''))}</td>"
        "</tr>"
        for row in phylogenetic_rows[:10]
    )
    phylogenetic_section = ""
    if phylogenetic_rows:
        phylogenetic_warning = _phylogenetic_warning_html(phylogenetic_rows)
        phylogenetic_section = f"""
      <section class="report-section">
        <h2>Closest Reference Genomes</h2>
        {phylogenetic_warning}
        <p class="section-intro">Marker similarity identifies the candidate group. Exact assembly ties are resolved by canonical genome identity and MUMmer whole-genome SNP comparison.</p>
        <div class="table-scroll"><table class="primary-table">
          <thead><tr><th>Rank</th><th>Reference</th><th>Marker status</th><th>Marker distance</th><th>Exact genome</th><th>WG SNPs</th><th>Indel bases</th><th>Aligned % ref/query</th><th>Date</th><th>Location</th></tr></thead>
          <tbody>{reference_summary_rows}</tbody>
        </table></div>
        <details>
          <summary>Technical marker-distance components</summary>
          <p class="terminal-note">EPA/tree, direct aligned-sequence, repeat, and tie-break components used to construct the ranking.</p>
          <div class="table-scroll"><table>
            <thead><tr><th>Rank</th><th>Reference</th><th>Combined distance</th><th>Hybrid SNP</th><th>EPA/tree SNP</th><th>Direct SNP</th><th>Normalized repeat</th><th>Compared loci</th><th>Exact marker loci</th><th>Match status</th><th>Exact genome</th><th>WG SNPs</th><th>Indel bases</th><th>Ref AF</th><th>Query AF</th><th>Tie break</th><th>Gap to next</th><th>Date</th><th>Location</th></tr></thead>
            <tbody>{phylogenetic_table_rows}</tbody>
          </table></div>
        </details>
      </section>
"""

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>MLVAMaps Assembly Report - {_safe(sample_id)}</title>
  <style>
    :root {{
      color-scheme: dark;
      --screen: #07150f;
      --panel: #0b2117;
      --phosphor: #62ff9b;
      --amber: #ffc857;
      --cyan: #5fe8ff;
      --magenta: #ff5fc8;
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
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
      letter-spacing: 0;
    }}
    main {{ max-width: 1320px; margin: 0 auto; padding: 2rem; }}
    h1, h2 {{ color: var(--amber); }}
    h1 {{ font-size: clamp(1.7rem, 4vw, 3rem); margin: 0 0 0.35rem; }}
    h2 {{ font-size: 1.1rem; margin-top: 2rem; }}
    .subhead, .terminal-note {{ color: var(--muted); }}
    .warning-banner {{ margin: 0.75rem 0; padding: 0.8rem 1rem; color: #1b1200; background: var(--amber); border: 2px solid #ff8c42; border-radius: 5px; }}
    .terminal {{
      border: 2px solid var(--line);
      background: linear-gradient(180deg, rgba(11, 33, 23, 0.94), rgba(4, 13, 9, 0.94));
      box-shadow: 0 0 0 1px rgba(98,255,155,0.16), 0 0 28px rgba(98,255,155,0.12);
      border-radius: 8px;
      padding: 1rem;
      overflow-x: auto;
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
    .metric small {{ display: block; color: var(--muted); line-height: 1.35; margin-top: 0.35rem; }}
    .metric.good {{ border-color: rgba(98,255,155,0.65); }}
    .metric.warn {{ border-color: rgba(255,200,87,0.75); }}
    .findings {{ display: grid; gap: 0.55rem; margin: 0 0 1.5rem; }}
    .finding {{ display: grid; grid-template-columns: minmax(150px, 220px) 1fr; gap: 0.8rem; padding: 0.75rem 0.9rem; border-left: 4px solid var(--cyan); background: rgba(95,232,255,0.06); border-radius: 4px; }}
    .finding strong {{ color: var(--cyan); }}
    .finding span {{ color: #d8fbe2; }}
    .finding.warn {{ border-left-color: var(--amber); background: rgba(255,200,87,0.07); }}
    .finding.warn strong {{ color: var(--amber); }}
    .finding.good {{ border-left-color: var(--phosphor); }}
    .finding.good strong {{ color: var(--phosphor); }}
    .report-section {{ margin-top: 2.2rem; padding-top: 0.2rem; }}
    .section-intro {{ color: var(--muted); max-width: 72rem; line-height: 1.5; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 0.75rem; min-width: 900px; }}
    th, td {{ border-bottom: 1px solid var(--line); padding: 0.55rem; text-align: left; vertical-align: top; }}
    th {{ color: var(--cyan); font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.04em; }}
    td {{ color: #d8fbe2; }}
    .primary-table tbody tr:first-child {{ background: rgba(98,255,155,0.08); }}
    .gel-panel {{ margin: 1rem 0 1.5rem; }}
    .gel-panel svg {{ width: 100%; max-height: 620px; display: block; }}
    .gel-panel figcaption {{ color: var(--muted); font-size: 0.9rem; margin-top: 0.5rem; }}
    .gel-frame {{ fill: #080817; stroke: #395089; stroke-width: 2; }}
    .well-line {{ stroke: #7e8cff; stroke-width: 1.2; opacity: 0.55; }}
    .well {{ fill: #03030c; stroke: #7e8cff; opacity: 0.75; }}
    .marker-band {{ fill: #b6d8ff; filter: url(#assembly-glow); opacity: 0.92; }}
    .query-band {{ fill: var(--phosphor); filter: url(#assembly-glow); }}
    .reference-band {{ fill: var(--magenta); filter: url(#assembly-glow); opacity: 0.88; }}
    .lane-title {{ fill: var(--amber); text-anchor: middle; font: 700 15px "Courier New", monospace; }}
    .gel-legend {{ fill: var(--muted); text-anchor: middle; font: 12px "Courier New", monospace; }}
    .band-hit .band-label {{ fill: #d6ffe2; font: 10px "Courier New", monospace; opacity: 0; pointer-events: none; transition: opacity 120ms ease; }}
    .band-hit:hover .band-label, .band-hit:focus .band-label {{ opacity: 0.95; }}
    .band-hit:focus {{ outline: none; }}
    .reference-label {{ fill: #ffd5f3; }}
    .marker-label {{ fill: #b6d8ff; font: 11px "Courier New", monospace; text-anchor: end; }}
    .chart-scroll {{ overflow-x: auto; }}
    .chart-panel {{ margin: 0.75rem 0 1.5rem; border: 1px solid var(--line); border-radius: 7px; padding: 0.6rem; background: rgba(3, 12, 8, 0.58); }}
    .chart-panel svg {{ width: 100%; display: block; min-width: 720px; }}
    .chart-panel figcaption {{ color: var(--muted); font-size: 0.88rem; margin-top: 0.45rem; }}
    .chart-label {{ fill: #d6ffe2; font: 12px "Courier New", monospace; }}
    .chart-value {{ fill: var(--muted); font: 11px "Courier New", monospace; }}
    .mapping-track {{ fill: rgba(141, 215, 170, 0.1); stroke: var(--line); stroke-width: 1; }}
    details {{ border: 1px solid var(--line); border-radius: 6px; margin: 0.7rem 0; padding: 0.65rem 0.8rem; background: rgba(98,255,155,0.035); }}
    summary {{ color: var(--cyan); cursor: pointer; font-weight: 700; }}
    .table-scroll {{ overflow-x: auto; }}
    @media (max-width: 700px) {{
      main {{ padding: 1rem; }}
      .finding {{ grid-template-columns: 1fr; gap: 0.2rem; }}
    }}
  </style>
</head>
<body>
  <main>
    <h1>MLVAMaps Assembly Report: {_safe(sample_id)}</h1>
    <p class="subhead">Assembly-derived VNTR calls, panel quality, and closest-reference interpretation.</p>
    <section class="terminal">
      <h2>Sample Overview</h2>
      <div class="summary">
        {summary_html}
      </div>
      <div class="findings">{findings_html}</div>
      <h2>Individual Locus Repeat Counts</h2>
      <div class="chart-scroll">{repeat_count_plot}</div>
      <h2>Generated Gel</h2>
      {gel}
      {profile_section}
      {phylogenetic_section}
      <section class="report-section">
        <h2>Detailed Evidence</h2>
        <p class="section-intro">Expand these tables when exact product coordinates or alternative-call probabilities are needed.</p>
        <details>
          <summary>Locus call details</summary>
          <div class="table-scroll"><table>
            <thead><tr><th>Locus</th><th>Present</th><th>Repeat count</th><th>Raw count</th><th>Product bp</th><th>Reads</th><th>Coverage</th><th>Confidence</th><th>Second allele</th><th>Second probability</th><th>Inference</th><th>Status</th><th>Evidence</th></tr></thead>
            <tbody>{table_rows}</tbody>
          </table></div>
        </details>
        <details>
          <summary>Assembly Amplicons: coordinates and primer mismatches</summary>
          <div class="table-scroll"><table>
            <thead><tr><th>Locus</th><th>Contig</th><th>Coordinates</th><th>Orientation</th><th>Product bp</th><th>Forward mismatches</th><th>Reverse mismatches</th></tr></thead>
            <tbody>{product_table_rows}</tbody>
          </table></div>
        </details>
      </section>
    </section>
  </main>
</body>
</html>
"""
    (outdir / "report.html").write_text(html)
