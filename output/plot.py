from pathlib import Path
import bisect
import struct
import zlib


BASE_DIR = Path(__file__).resolve().parent
DATA_PATH = BASE_DIR / "result.dat"
PNG_PATH = BASE_DIR / "contour.png"
SVG_PATH = BASE_DIR / "contour.svg"


def read_result(path):
    points = []

    with path.open() as f:
        for line in f:
            values = line.split()

            if len(values) < 3:
                continue

            x, y, phi = map(float, values[:3])
            points.append((x, y, phi))

    return points


def to_grid(points):
    xs = sorted({point[0] for point in points})
    ys = sorted({point[1] for point in points})
    values = {(x, y): phi for x, y, phi in points}
    phi_grid = [[values[(x, y)] for x in xs] for y in ys]

    return xs, ys, phi_grid


def tick_values(vmin, vmax, intervals):
    return [vmin + (vmax - vmin) * i / intervals for i in range(intervals + 1)]


def format_tick(value):
    if abs(value - round(value)) < 1.0e-10:
        return str(round(value))

    return f"{value:.2f}".rstrip("0").rstrip(".")


def plot_with_matplotlib(xs, ys, phi_grid):
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 3.2), constrained_layout=True)

    contour = ax.contourf(xs, ys, phi_grid, levels=30, cmap="viridis")
    ax.contour(xs, ys, phi_grid, levels=15, colors="black", linewidths=0.4)

    fig.colorbar(contour, ax=ax, label="Φ")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Laplace Equation FEM Result")
    ax.set_xticks(tick_values(min(xs), max(xs), 5))
    ax.set_yticks(tick_values(min(ys), max(ys), 4))
    ax.set_aspect("equal", adjustable="box")

    fig.savefig(PNG_PATH, dpi=300)
    print(f"saved: {PNG_PATH}")


def viridis_rgb(t):
    colors = [
        (68, 1, 84),
        (59, 82, 139),
        (33, 145, 140),
        (94, 201, 98),
        (253, 231, 37),
    ]

    t = min(1.0, max(0.0, t))
    scaled = t * (len(colors) - 1)
    i = min(int(scaled), len(colors) - 2)
    local = scaled - i

    r0, g0, b0 = colors[i]
    r1, g1, b1 = colors[i + 1]

    r = round(r0 + (r1 - r0) * local)
    g = round(g0 + (g1 - g0) * local)
    b = round(b0 + (b1 - b0) * local)

    return r, g, b


def viridis_like(t):
    r, g, b = viridis_rgb(t)

    return f"#{r:02x}{g:02x}{b:02x}"


def project(x, y, xmin, xmax, ymin, ymax, width, height, margin):
    px = margin + (x - xmin) / (xmax - xmin) * (width - 2 * margin)
    py = height - margin - (y - ymin) / (ymax - ymin) * (height - 2 * margin)

    return px, py


def interpolate(p0, p1, v0, v1, level):
    if abs(v1 - v0) < 1.0e-14:
        t = 0.5
    else:
        t = (level - v0) / (v1 - v0)

    return (
        p0[0] + (p1[0] - p0[0]) * t,
        p0[1] + (p1[1] - p0[1]) * t,
    )


def contour_segments(xs, ys, phi_grid, levels):
    segments = []

    for iy in range(len(ys) - 1):
        for ix in range(len(xs) - 1):
            p = [
                (xs[ix], ys[iy]),
                (xs[ix + 1], ys[iy]),
                (xs[ix + 1], ys[iy + 1]),
                (xs[ix], ys[iy + 1]),
            ]
            v = [
                phi_grid[iy][ix],
                phi_grid[iy][ix + 1],
                phi_grid[iy + 1][ix + 1],
                phi_grid[iy + 1][ix],
            ]

            for level in levels:
                crossings = []

                for i0, i1 in ((0, 1), (1, 2), (2, 3), (3, 0)):
                    v0 = v[i0]
                    v1 = v[i1]

                    if (v0 - level) * (v1 - level) < 0.0:
                        crossings.append(interpolate(p[i0], p[i1], v0, v1, level))

                if len(crossings) == 2:
                    segments.append((crossings[0], crossings[1]))
                elif len(crossings) == 4:
                    segments.append((crossings[0], crossings[1]))
                    segments.append((crossings[2], crossings[3]))

    return segments


def plot_with_svg(xs, ys, phi_grid):
    width = 1060
    height = 340
    margin = 75

    xmin = min(xs)
    xmax = max(xs)
    ymin = min(ys)
    ymax = max(ys)
    phis = [value for row in phi_grid for value in row]
    pmin = min(phis)
    pmax = max(phis)

    def color(value):
        return viridis_like((value - pmin) / (pmax - pmin))

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="22">Laplace Equation FEM Result</text>',
    ]

    for iy in range(len(ys) - 1):
        for ix in range(len(xs) - 1):
            x0, x1 = xs[ix], xs[ix + 1]
            y0, y1 = ys[iy], ys[iy + 1]
            value = (
                phi_grid[iy][ix]
                + phi_grid[iy][ix + 1]
                + phi_grid[iy + 1][ix]
                + phi_grid[iy + 1][ix + 1]
            ) / 4.0

            polygon = [
                project(x0, y0, xmin, xmax, ymin, ymax, width, height, margin),
                project(x1, y0, xmin, xmax, ymin, ymax, width, height, margin),
                project(x1, y1, xmin, xmax, ymin, ymax, width, height, margin),
                project(x0, y1, xmin, xmax, ymin, ymax, width, height, margin),
            ]
            points = " ".join(f"{px:.2f},{py:.2f}" for px, py in polygon)
            lines.append(f'<polygon points="{points}" fill="{color(value)}" stroke="none"/>')

    levels = [pmin + (pmax - pmin) * i / 10.0 for i in range(1, 10)]

    for p0, p1 in contour_segments(xs, ys, phi_grid, levels):
        x0, y0 = project(p0[0], p0[1], xmin, xmax, ymin, ymax, width, height, margin)
        x1, y1 = project(p1[0], p1[1], xmin, xmax, ymin, ymax, width, height, margin)
        lines.append(
            f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" '
            'stroke="black" stroke-width="1" opacity="0.65"/>'
        )

    x_axis_y = height - margin
    y_axis_x = margin
    lines.append(f'<line x1="{margin}" y1="{x_axis_y}" x2="{width - margin}" y2="{x_axis_y}" stroke="black"/>')
    lines.append(f'<line x1="{y_axis_x}" y1="{margin}" x2="{y_axis_x}" y2="{height - margin}" stroke="black"/>')

    for tick in tick_values(xmin, xmax, 5):
        x, _ = project(tick, ymin, xmin, xmax, ymin, ymax, width, height, margin)
        lines.append(f'<line x1="{x:.2f}" y1="{margin}" x2="{x:.2f}" y2="{x_axis_y}" stroke="black" opacity="0.18"/>')
        lines.append(f'<line x1="{x:.2f}" y1="{x_axis_y}" x2="{x:.2f}" y2="{x_axis_y + 6}" stroke="black"/>')
        lines.append(
            f'<text x="{x:.2f}" y="{x_axis_y + 22}" text-anchor="middle" '
            f'font-family="sans-serif" font-size="13">{format_tick(tick)}</text>'
        )

    for tick in tick_values(ymin, ymax, 4):
        _, y = project(xmin, tick, xmin, xmax, ymin, ymax, width, height, margin)
        lines.append(f'<line x1="{margin}" y1="{y:.2f}" x2="{width - margin}" y2="{y:.2f}" stroke="black" opacity="0.18"/>')
        lines.append(f'<line x1="{y_axis_x - 6}" y1="{y:.2f}" x2="{y_axis_x}" y2="{y:.2f}" stroke="black"/>')
        lines.append(
            f'<text x="{y_axis_x - 10}" y="{y + 4:.2f}" text-anchor="end" '
            f'font-family="sans-serif" font-size="13">{format_tick(tick)}</text>'
        )

    lines.append(f'<text x="{width / 2}" y="{height - 12}" text-anchor="middle" font-family="sans-serif" font-size="16">x</text>')
    lines.append(f'<text x="16" y="{height / 2}" text-anchor="middle" font-family="sans-serif" font-size="16" transform="rotate(-90 16 {height / 2})">y</text>')

    bar_x = width - 55
    bar_y = margin
    bar_w = 16
    bar_h = height - 2 * margin

    for i in range(80):
        t0 = i / 80.0
        y = bar_y + bar_h * (1.0 - (i + 1) / 80.0)
        lines.append(
            f'<rect x="{bar_x}" y="{y:.2f}" width="{bar_w}" height="{bar_h / 80.0 + 0.5:.2f}" '
            f'fill="{viridis_like(t0)}"/>'
        )

    for tick in tick_values(pmin, pmax, 4):
        t = (tick - pmin) / (pmax - pmin)
        y = bar_y + bar_h * (1.0 - t)
        lines.append(f'<line x1="{bar_x + bar_w}" y1="{y:.2f}" x2="{bar_x + bar_w + 5}" y2="{y:.2f}" stroke="black"/>')
        lines.append(
            f'<text x="{bar_x + bar_w + 9}" y="{y + 4:.2f}" font-family="sans-serif" '
            f'font-size="13">{format_tick(tick)}</text>'
        )

    lines.append(f'<text x="{bar_x - 4}" y="{bar_y - 10}" text-anchor="middle" font-family="sans-serif" font-size="14">Φ</text>')

    lines.append("</svg>")

    SVG_PATH.write_text("\n".join(lines))
    print(f"saved: {SVG_PATH}")


def write_png(path, width, height, pixels):
    def chunk(kind, data):
        return (
            struct.pack(">I", len(data))
            + kind
            + data
            + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)
        )

    raw = b"".join(b"\x00" + bytes(row) for row in pixels)

    with path.open("wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n")
        f.write(chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)))
        f.write(chunk(b"IDAT", zlib.compress(raw, 9)))
        f.write(chunk(b"IEND", b""))


def set_pixel(pixels, x, y, color):
    height = len(pixels)
    width = len(pixels[0]) // 3

    if 0 <= x < width and 0 <= y < height:
        i = x * 3
        pixels[y][i : i + 3] = bytes(color)


FONT = {
    "0": ["111", "101", "101", "101", "111"],
    "1": ["010", "110", "010", "010", "111"],
    "2": ["111", "001", "111", "100", "111"],
    "3": ["111", "001", "111", "001", "111"],
    "4": ["101", "101", "111", "001", "001"],
    "5": ["111", "100", "111", "001", "111"],
    "6": ["111", "100", "111", "101", "111"],
    "7": ["111", "001", "010", "010", "010"],
    "8": ["111", "101", "111", "101", "111"],
    "9": ["111", "101", "111", "001", "111"],
    ".": ["0", "0", "0", "0", "1"],
    "-": ["000", "000", "111", "000", "000"],
    "P": ["110", "101", "110", "100", "100"],
    "h": ["100", "100", "110", "101", "101"],
    "i": ["1", "0", "1", "1", "1"],
    "p": ["000", "110", "101", "110", "100"],
    "Φ": ["01110", "10101", "11111", "10101", "01110"],
    "x": ["101", "101", "010", "101", "101"],
    "y": ["101", "101", "111", "001", "111"],
}


def text_size(text, scale):
    width = 0

    for index, char in enumerate(text):
        glyph = FONT.get(char, ["000", "000", "000", "000", "000"])
        width += len(glyph[0]) * scale

        if index != len(text) - 1:
            width += scale

    return width, 5 * scale


def draw_text(pixels, x, y, text, color=(0, 0, 0), scale=2, anchor="left"):
    width, _ = text_size(text, scale)

    if anchor == "center":
        x -= width // 2
    elif anchor == "right":
        x -= width

    cursor = round(x)
    y = round(y)

    for char in text:
        glyph = FONT.get(char, ["000", "000", "000", "000", "000"])

        for gy, row in enumerate(glyph):
            for gx, value in enumerate(row):
                if value != "1":
                    continue

                for sy in range(scale):
                    for sx in range(scale):
                        set_pixel(pixels, cursor + gx * scale + sx, y + gy * scale + sy, color)

        cursor += (len(glyph[0]) + 1) * scale


def draw_line(pixels, x0, y0, x1, y1, color, thickness=1):
    x0 = round(x0)
    y0 = round(y0)
    x1 = round(x1)
    y1 = round(y1)

    dx = abs(x1 - x0)
    dy = -abs(y1 - y0)
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    err = dx + dy

    while True:
        for oy in range(-(thickness // 2), thickness // 2 + 1):
            for ox in range(-(thickness // 2), thickness // 2 + 1):
                set_pixel(pixels, x0 + ox, y0 + oy, color)

        if x0 == x1 and y0 == y1:
            break

        e2 = 2 * err

        if e2 >= dy:
            err += dy
            x0 += sx

        if e2 <= dx:
            err += dx
            y0 += sy


def grid_value(xs, ys, phi_grid, x, y):
    ix = bisect.bisect_right(xs, x) - 1
    iy = bisect.bisect_right(ys, y) - 1

    ix = max(0, min(ix, len(xs) - 2))
    iy = max(0, min(iy, len(ys) - 2))

    x0 = xs[ix]
    x1 = xs[ix + 1]
    y0 = ys[iy]
    y1 = ys[iy + 1]

    tx = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
    ty = 0.0 if y1 == y0 else (y - y0) / (y1 - y0)

    v00 = phi_grid[iy][ix]
    v10 = phi_grid[iy][ix + 1]
    v01 = phi_grid[iy + 1][ix]
    v11 = phi_grid[iy + 1][ix + 1]

    return (
        v00 * (1.0 - tx) * (1.0 - ty)
        + v10 * tx * (1.0 - ty)
        + v01 * (1.0 - tx) * ty
        + v11 * tx * ty
    )


def plot_with_png(xs, ys, phi_grid):
    width = 1060
    height = 340
    margin = 75

    pixels = [bytearray([255, 255, 255] * width) for _ in range(height)]

    xmin = min(xs)
    xmax = max(xs)
    ymin = min(ys)
    ymax = max(ys)
    phis = [value for row in phi_grid for value in row]
    pmin = min(phis)
    pmax = max(phis)

    plot_left = margin
    plot_right = width - margin
    plot_top = margin
    plot_bottom = height - margin

    for py in range(plot_top, plot_bottom + 1):
        y = ymin + (plot_bottom - py) / (plot_bottom - plot_top) * (ymax - ymin)

        for px in range(plot_left, plot_right + 1):
            x = xmin + (px - plot_left) / (plot_right - plot_left) * (xmax - xmin)
            value = grid_value(xs, ys, phi_grid, x, y)
            t = (value - pmin) / (pmax - pmin)
            set_pixel(pixels, px, py, viridis_rgb(t))

    for tick in tick_values(xmin, xmax, 5):
        x, _ = project(tick, ymin, xmin, xmax, ymin, ymax, width, height, margin)
        draw_line(pixels, x, plot_top, x, plot_bottom, (0, 0, 0), thickness=1)
        draw_line(pixels, x, plot_bottom, x, plot_bottom + 7, (0, 0, 0), thickness=2)
        draw_text(pixels, x, plot_bottom + 13, format_tick(tick), anchor="center")

    for tick in tick_values(ymin, ymax, 4):
        _, y = project(xmin, tick, xmin, xmax, ymin, ymax, width, height, margin)
        draw_line(pixels, plot_left, y, plot_right, y, (0, 0, 0), thickness=1)
        draw_line(pixels, plot_left - 7, y, plot_left, y, (0, 0, 0), thickness=2)
        draw_text(pixels, plot_left - 12, y - 5, format_tick(tick), anchor="right")

    draw_text(pixels, (plot_left + plot_right) // 2, height - 18, "x", anchor="center")
    draw_text(pixels, 24, (plot_top + plot_bottom) // 2 - 5, "y", anchor="center")

    levels = [pmin + (pmax - pmin) * i / 10.0 for i in range(1, 10)]

    for p0, p1 in contour_segments(xs, ys, phi_grid, levels):
        x0, y0 = project(p0[0], p0[1], xmin, xmax, ymin, ymax, width, height, margin)
        x1, y1 = project(p1[0], p1[1], xmin, xmax, ymin, ymax, width, height, margin)
        draw_line(pixels, x0, y0, x1, y1, (0, 0, 0), thickness=2)

    draw_line(pixels, plot_left, plot_bottom, plot_right, plot_bottom, (0, 0, 0), thickness=2)
    draw_line(pixels, plot_left, plot_top, plot_left, plot_bottom, (0, 0, 0), thickness=2)

    bar_x0 = width - 55
    bar_x1 = width - 39

    for py in range(plot_top, plot_bottom + 1):
        t = (plot_bottom - py) / (plot_bottom - plot_top)

        for px in range(bar_x0, bar_x1 + 1):
            set_pixel(pixels, px, py, viridis_rgb(t))

    draw_line(pixels, bar_x0, plot_top, bar_x1, plot_top, (0, 0, 0))
    draw_line(pixels, bar_x1, plot_top, bar_x1, plot_bottom, (0, 0, 0))
    draw_line(pixels, bar_x1, plot_bottom, bar_x0, plot_bottom, (0, 0, 0))
    draw_line(pixels, bar_x0, plot_bottom, bar_x0, plot_top, (0, 0, 0))

    for tick in tick_values(pmin, pmax, 4):
        t = (tick - pmin) / (pmax - pmin)
        y = plot_bottom - t * (plot_bottom - plot_top)
        draw_line(pixels, bar_x1, y, bar_x1 + 6, y, (0, 0, 0), thickness=1)
        draw_text(pixels, bar_x1 + 10, y - 5, format_tick(tick), anchor="left")

    draw_text(pixels, (bar_x0 + bar_x1) // 2, plot_top - 19, "Φ", anchor="center")

    write_png(PNG_PATH, width, height, pixels)
    print(f"saved: {PNG_PATH}")


def main():
    points = read_result(DATA_PATH)
    xs, ys, phi_grid = to_grid(points)

    try:
        plot_with_matplotlib(xs, ys, phi_grid)
    except ModuleNotFoundError:
        print("matplotlib is not installed. Creating PNG and SVG with the Python standard library.")
        plot_with_png(xs, ys, phi_grid)
        plot_with_svg(xs, ys, phi_grid)


if __name__ == "__main__":
    main()
