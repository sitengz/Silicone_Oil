#!/usr/bin/env python3
"""Build explicit PDMS chain-count configs and reporting data.

This is an offline case-design tool. The C++ simulation generator reads only
the resulting ``chain_count = length count`` rows and knows nothing about PDI.
"""

from __future__ import annotations

import csv
import io
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
MONO_LENGTHS = (4, 8, 16, 32, 64, 128)
MEANS = (16, 32, 64)
TARGETS = (1.05, 1.10, 1.15, 1.20, 1.25, 1.30)
COLORS = ("#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE")


def case_name(mean: int, target: float) -> str:
    return f"N{mean}_PDI{target:.2f}".rstrip("0").rstrip(".")


def probabilities(lower: int, upper: int, shape: float, rate: float) -> list[float]:
    logs = [(shape - 1.0) * math.log(n) - rate * n
            for n in range(lower, upper + 1)]
    largest = max(logs)
    weights = [math.exp(value - largest) for value in logs]
    total = sum(weights)
    return [weight / total for weight in weights]


def moments(lower: int, probabilities_by_length: list[float]) -> tuple[float, float]:
    mean = sum((lower + i) * p for i, p in enumerate(probabilities_by_length))
    second = sum((lower + i) ** 2 * p for i, p in enumerate(probabilities_by_length))
    return mean, second / (mean * mean)


def fit_sz(mean: int, target: float, lower: int, upper: int) -> tuple[float, float, list[float]]:
    """Fit a truncated, discrete number distribution proportional to N^(k-1)e^(-rate*N)."""
    def at_shape(shape: float) -> tuple[float, float, list[float]]:
        at_zero = probabilities(lower, upper, shape, 0.0)
        if moments(lower, at_zero)[0] < mean - 1e-10:
            raise ValueError("This shape cannot reach the requested mean with positive scale")
        lo, hi = 0.0, 1.0
        while moments(lower, probabilities(lower, upper, shape, hi))[0] > mean:
            hi *= 2.0
        for _ in range(85):
            mid = (lo + hi) / 2.0
            if moments(lower, probabilities(lower, upper, shape, mid))[0] > mean:
                lo = mid
            else:
                hi = mid
        rate = (lo + hi) / 2.0
        p = probabilities(lower, upper, shape, rate)
        return moments(lower, p)[1], rate, p

    lo, hi = 1e-8, 10.0
    for _ in range(85):
        mid = (lo + hi) / 2.0
        if moments(lower, probabilities(lower, upper, mid, 0.0))[0] < mean:
            lo = mid
        else:
            hi = mid
    minimum_shape = (lo + hi) / 2.0
    maximum_pdi, _, _ = at_shape(minimum_shape)
    if target > maximum_pdi + 1e-9:
        raise ValueError(f"PDI {target} exceeds the SZ limit {maximum_pdi:.6f} for N={lower}..{upper}")
    lo, hi = minimum_shape, 16.0
    while at_shape(hi)[0] > target:
        hi *= 2.0
    for _ in range(85):
        mid = (lo + hi) / 2.0
        if at_shape(mid)[0] > target:
            lo = mid
        else:
            hi = mid
    shape = (lo + hi) / 2.0
    _, rate, p = at_shape(shape)
    return shape, rate, p


def statistics(counts: dict[int, int]) -> tuple[int, int, float, float, float]:
    chains = sum(counts.values())
    repeats = sum(length * count for length, count in counts.items())
    second = sum(length * length * count for length, count in counts.items())
    mn = repeats / chains
    mw = second / repeats
    return chains, repeats, mn, mw, mw / mn


def is_unimodal(values: list[int]) -> bool:
    index = 0
    while index + 1 < len(values) and values[index + 1] >= values[index]:
        index += 1
    return all(values[i + 1] <= values[i] for i in range(index, len(values) - 1))


def smooth_histogram_at_shape(mean: int, lower: int, upper: int, chains: int,
                              shape: float) -> tuple[dict[int, int], float, list[float]] | None:
    """Round an SZ curve to one continuous, unimodal integer-count profile."""
    zero_mean = moments(lower, probabilities(lower, upper, shape, 0.0))[0]
    if zero_mean < mean - 1e-10:
        return None
    low, high = 0.0, 1.0
    while moments(lower, probabilities(lower, upper, shape, high))[0] > mean:
        high *= 2.0
    for _ in range(55):
        rate = (low + high) / 2.0
        if moments(lower, probabilities(lower, upper, shape, rate))[0] > mean:
            low = rate
        else:
            high = rate
    rate = (low + high) / 2.0
    p = probabilities(lower, upper, shape, rate)
    quotas = [chains * value for value in p]
    values = [round(quota) for quota in quotas]
    if not is_unimodal(values):
        return None

    target_repeats = chains * mean
    actual_chains = sum(values)
    actual_repeats = sum((lower + i) * count for i, count in enumerate(values))
    while actual_chains != chains:
        direction = 1 if actual_chains < chains else -1
        remaining = chains - actual_chains
        desired_length = (target_repeats - actual_repeats) / remaining
        options = []
        for i, count in enumerate(values):
            if count + direction < 0:
                continue
            values[i] += direction
            if is_unimodal(values):
                curve_error = (values[i] - quotas[i]) ** 2 - (count - quotas[i]) ** 2
                options.append(((lower + i - desired_length) ** 2 +
                                0.05 * curve_error, i))
            values[i] -= direction
        if not options:
            return None
        _, index = min(options)
        values[index] += direction
        actual_chains += direction
        actual_repeats += direction * (lower + index)

    while actual_repeats != target_repeats:
        direction = 1 if actual_repeats < target_repeats else -1
        options = []
        for i, count in enumerate(values):
            neighbor = i + direction
            if count < 1 or neighbor < 0 or neighbor >= len(values):
                continue
            values[i] -= 1
            values[neighbor] += 1
            if is_unimodal(values):
                curve_error = sum((values[j] - quotas[j]) ** 2
                                  for j in (i, neighbor))
                options.append((curve_error, i, neighbor))
            values[i] += 1
            values[neighbor] -= 1
        if not options:
            return None
        _, index, neighbor = min(options)
        values[index] -= 1
        values[neighbor] += 1
        actual_repeats += direction

    occupied = [lower + i for i, count in enumerate(values) if count]
    if occupied != list(range(occupied[0], occupied[-1] + 1)):
        return None
    result = {lower + i: count for i, count in enumerate(values) if count}
    return result, rate, p


def refine_pdi(counts: dict[int, int], lower: int, upper: int,
               target: float, chains: int, mean: int) -> dict[int, int]:
    """Move two chains together/apart, preserving chains, Mn, and unimodality."""
    values = [counts.get(n, 0) for n in range(lower, upper + 1)]
    denominator = chains * mean * mean
    desired_second = target * denominator
    second = sum((lower + i) ** 2 * count for i, count in enumerate(values))
    for _ in range(5):
        if abs(second - desired_second) / denominator < 0.0001:
            break
        increase = second < desired_second
        best = None
        for a in range(len(values)):
            for b in range(a + 3, len(values)):
                for c in range(a + 1, (a + b) // 2 + 1):
                    d = a + b - c
                    sources, destinations = ((c, d), (a, b)) if increase else ((a, b), (c, d))
                    changes: dict[int, int] = {}
                    for index in sources:
                        changes[index] = changes.get(index, 0) - 1
                    for index in destinations:
                        changes[index] = changes.get(index, 0) + 1
                    if any(values[index] + delta < 0 for index, delta in changes.items()):
                        continue
                    delta_second = sum((lower + index) ** 2 * delta
                                       for index, delta in changes.items())
                    error = abs(second + delta_second - desired_second)
                    if error >= abs(second - desired_second) or (best and error >= best[0]):
                        continue
                    for index, delta in changes.items():
                        values[index] += delta
                    valid = is_unimodal(values)
                    for index, delta in changes.items():
                        values[index] -= delta
                    if valid:
                        best = error, changes, delta_second
        if best is None:
            break
        _, changes, delta_second = best
        for index, delta in changes.items():
            values[index] += delta
        second += delta_second
    return {lower + i: value for i, value in enumerate(values) if value}


def fit_integer_histogram(mean: int, target: float, lower: int, upper: int,
                          chains: int) -> tuple[dict[int, int], float, float]:
    shape, _, _ = fit_sz(mean, target, lower, upper)
    best = None
    # Small shape offsets compensate for integer rounding while the continuous
    # SZ distribution remains the template. The two-chain adjustment then
    # fine-tunes PDI without changing Mn or introducing holes/oscillations.
    for offset in (0, -0.002, 0.002, -0.005, 0.005, -0.01, 0.01,
                   -0.02, 0.02, -0.04, 0.04):
        candidate = smooth_histogram_at_shape(mean, lower, upper, chains,
                                               shape * (1 + offset))
        if candidate is None:
            continue
        counts, rate, _ = candidate
        counts = refine_pdi(counts, lower, upper, target, chains, mean)
        error = abs(statistics(counts)[4] - target)
        score = error, abs(offset)
        if best is None or score < best[0]:
            best = score, counts, shape * (1 + offset), rate
        if error < 0.0001:
            break
    if best is None or best[0][0] >= 0.0001:
        raise ValueError(f"Cannot fit a smooth integer histogram for Mn={mean}, PDI={target}")
    _, counts, shape, rate = best
    return counts, shape, rate


def write_config(name: str, counts: dict[int, int]) -> None:
    folder = HERE / name
    folder.mkdir(exist_ok=True)
    rows = [
        "# PDMS chain-length/count histogram. The generator reads these rows directly.",
        "# Format: chain_count = number_of_DMS_repeats number_of_chains",
    ]
    rows.extend(f"chain_count = {length} {count}" for length, count in counts.items())
    rows.extend((
        "mps_percent = 0",
        "sequence = random",
        "density = 0.1",
        "target_density = 0.8",
        "min_separation = 4.5",
        "seed = 20260727",
        "velocity_seed = 492845",
        f"output = data.{name}",
    ))
    (folder / "model.conf").write_text("\n".join(rows) + "\n")


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    handle = io.StringIO()
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    content = handle.getvalue()
    if not path.exists() or path.read_text() != content:
        path.write_text(content)


def markdown_table(rows: list[dict[str, object]]) -> list[str]:
    keys = list(rows[0])
    result = ["| " + " | ".join(keys) + " |",
              "| " + " | ".join("---" for _ in keys) + " |"]
    result.extend("| " + " | ".join(str(row[key]) for key in keys) + " |" for row in rows)
    return result


def nice_y_axis(maximum: float) -> tuple[float, int]:
    """Choose a 1/2/5 tick interval with no more than six y labels."""
    if maximum <= 0:
        raise ValueError("The plotted chain counts must be positive")
    raw_step = maximum / 5
    magnitude = 10 ** math.floor(math.log10(raw_step))
    step = next(value * magnitude for value in (1, 2, 5, 10)
                if value * magnitude >= raw_step)
    intervals = math.ceil(maximum / step)
    return intervals * step, intervals


def write_svg(mean: int, series: list[tuple[float, dict[int, int], float, float]]) -> None:
    """Overlay integer counts, continuous SZ curves, and Mn/Mw references."""
    width, height = 1120, 840  # 4:3 figure canvas
    left, right, top, bottom = 90, 40, 30, 85
    x0 = 4
    x1 = 34 if mean == 16 else 2 * mean
    curves = []
    maximum = 0.0
    for target, counts, shape, rate in series:
        upper = 34 if mean == 16 and target == 1.30 else 2 * mean
        normalization = sum(n ** (shape - 1) * math.exp(-rate * n)
                            for n in range(x0, upper + 1))
        chains = sum(counts.values())
        samples = [(x0 + i * (upper - x0) / 240,
                    chains * (x0 + i * (upper - x0) / 240) ** (shape - 1)
                    * math.exp(-rate * (x0 + i * (upper - x0) / 240)) / normalization)
                   for i in range(241)]
        curves.append(samples)
        maximum = max(maximum, max(counts.values()), max(y for _, y in samples))
    y1, y_intervals = nice_y_axis(maximum)
    y_step = y1 / y_intervals
    plot_w, plot_h = width - left - right, height - top - bottom

    def sx(value: float) -> float:
        return left + (value - x0) * plot_w / (x1 - x0)

    def sy(value: float) -> float:
        return top + plot_h * (1.0 - value / y1)

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    for index in range(y_intervals + 1):
        tick = round(index * y_step)
        y = sy(tick)
        lines.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="#dedede"/>')
        lines.append(f'<text x="{left - 12}" y="{y + 5:.2f}" text-anchor="end" font-family="Arial" font-size="15" fill="#333">{tick}</text>')
    step = 5 if mean == 16 else 10 if mean == 32 else 20
    for tick in range(0, x1 + 1, step):
        if tick < x0:
            continue
        x = sx(tick)
        lines.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}" stroke="#eeeeee"/>')
        lines.append(f'<text x="{x:.2f}" y="{top + plot_h + 28}" text-anchor="middle" font-family="Arial" font-size="15" fill="#333">{tick}</text>')
    lines.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#333"/>')
    # Draw references first so curves and observed counts remain visible.
    mn_x = sx(mean)
    lines.append(f'<line x1="{mn_x:.2f}" y1="{top}" x2="{mn_x:.2f}" y2="{top + plot_h}" stroke="#222" stroke-width="1.6" stroke-dasharray="7 5"/>')
    for i, (_, counts, _, _) in enumerate(series):
        mw = statistics(counts)[3]
        x = sx(mw)
        lines.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}" stroke="{COLORS[i]}" stroke-width="1.4" stroke-dasharray="6 5" opacity="0.85"/>')
    for i, (target, counts, _, _) in enumerate(series):
        color = COLORS[i]
        points = " ".join(f"{sx(n):.2f},{sy(y):.2f}" for n, y in curves[i])
        lines.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.2" stroke-linejoin="round"/>')
        for n, count in counts.items():
            lines.append(f'<circle cx="{sx(n):.2f}" cy="{sy(count):.2f}" r="2.5" fill="{color}" stroke="white" stroke-width="0.6"/>')
        legend_y = top + 28 + i * 32
        legend_x = left + plot_w - 155
        lines.append(f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 42}" y2="{legend_y}" stroke="{color}" stroke-width="3"/>')
        lines.append(f'<circle cx="{legend_x + 21}" cy="{legend_y}" r="3.5" fill="{color}" stroke="white" stroke-width="0.8"/>')
        lines.append(f'<text x="{legend_x + 54}" y="{legend_y + 5}" font-family="Arial" font-size="17" fill="#222">PDI {target:.2f}</text>')
    lines.extend((
        f'<text x="{left + plot_w / 2}" y="{height - 25}" text-anchor="middle" font-family="DejaVu Math TeX Gyre, serif" font-style="italic" font-size="18" fill="#222">n</text>',
        f'<text x="25" y="{top + plot_h / 2}" text-anchor="middle" transform="rotate(-90 25 {top + plot_h / 2})" font-family="DejaVu Math TeX Gyre, serif" font-style="italic" font-size="18" fill="#222">M(n)</text>',
        '</svg>',
    ))
    (HERE / f"SZ_N{mean}.svg").write_text("\n".join(lines) + "\n")


def main() -> None:
    mono_rows = []
    for length in MONO_LENGTHS:
        chains = 100_000 // length
        counts = {length: chains}
        name = case_name(length, 1.0)
        write_config(name, counts)
        _, repeats, mn, mw, pdi = statistics(counts)
        mono_rows.append({"Case": name, "N": length, "Chains": chains,
                          "Beads": repeats, "Mn": f"{mn:.0f}",
                          "Mw": f"{mw:.0f}", "PDI": f"{pdi:.2f}"})

    write_csv(HERE / "table_PDI1.csv", mono_rows)
    markdown = ["# PDMS chain-length series", "",
                "One DMS repeat is one CG bead. Terminals are excluded from Mn and PDI.",
                "All model.conf files contain explicit chain_count rows; the C++ generator",
                "does not calculate or fit a dispersity.", "", "## PDI = 1.00", ""]
    markdown.extend(markdown_table(mono_rows))

    for mean in MEANS:
        chains = 99_999 // mean
        table = []
        plot_series = []
        for target in TARGETS:
            upper = 34 if mean == 16 and target == 1.30 else 2 * mean
            counts, shape, rate = fit_integer_histogram(mean, target, 4, upper, chains)
            name = case_name(mean, target)
            write_config(name, counts)
            actual_chains, repeats, mn, mw, pdi = statistics(counts)
            if actual_chains != chains or repeats != mean * chains:
                raise ValueError(f"Incorrect chain or repeat count for {name}")
            table.append({
                "Case": name,
                "Target PDI": f"{target:.2f}",
                "Realized PDI": f"{pdi:.6f}",
                "Mn": f"{mn:.0f}",
                "Mw": f"{mw:.6f}",
                "Chains": actual_chains,
                "Beads": repeats,
                "Occupied N": f"{min(counts)}-{max(counts)}",
                "SZ k": f"{shape:.6f}",
                "SZ theta": f"{1.0 / rate:.6f}",
                "SZ rate": f"{rate:.8f}",
            })
            plot_series.append((target, counts, shape, rate))
        write_csv(HERE / f"table_N{mean}.csv", table)
        markdown.extend(("", f"## Mn = {mean}", ""))
        markdown.extend(markdown_table(table))
        write_svg(mean, plot_series)
    (HERE / "pdms_series_tables.md").write_text("\n".join(markdown) + "\n")


if __name__ == "__main__":
    main()
