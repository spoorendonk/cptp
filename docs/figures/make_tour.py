"""Draw the docs site's problem figure, logo and favicon from a real cptp solve.

Run it from the repository root, with the extension built (`uv sync`):

    uv run python docs/figures/make_tour.py

It solves INSTANCE to optimality through the cptp Python API and writes, next
to itself, `tour-light.svg` and `tour-dark.svg` (the landing page's figure,
one per Furo colour mode), `logo-light.svg` and `logo-dark.svg` (the sidebar
logo: the tour cut out of the figure) and `favicon.svg`. The docs build uses
the committed SVGs and never runs this script.
"""

import math
from dataclasses import dataclass
from pathlib import Path

import cptp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
# A Roberti Set 3 instance: 75 customers spread over the plane, of which the
# optimum visits a clear subset. The smaller committed instances with
# coordinates visit a single customer.
INSTANCE = ROOT / "benchmarks" / "instances" / "roberti" / "E-n76-k7_a.vrp"
# One thread, so the tour is the same on every run.
OPTIONS = [("output_flag", "false"), ("threads", "1")]

FIGURE_WIDTH = 640.0
MARGIN = 16.0


@dataclass(frozen=True)
class Palette:
    """Colours of one Furo mode; `accent` is conf.py's color-brand-primary."""

    accent: str
    faded: str
    depot: str


LIGHT = Palette(accent="#0f766e", faded="#6b7280", depot="#111827")
DARK = Palette(accent="#2dd4bf", faded="#9ca3af", depot="#f3f4f6")
# A browser tab may be light or dark: a teal between the two modes'.
FAVICON_ACCENT = "#14b8a6"

Point = tuple[float, float]


def read_coords(path: Path) -> list[Point]:
    """The NODE_COORD_SECTION of a TSPLIB file, indexed as cptp numbers nodes."""
    coords: dict[int, Point] = {}
    in_section = False
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "NODE_COORD_SECTION":
            in_section = True
        elif in_section and fields[0].lstrip("-").isdigit():
            coords[int(fields[0]) - 1] = (float(fields[1]), float(fields[2]))
        elif in_section:
            break
    return [coords[i] for i in range(len(coords))]


@dataclass(frozen=True)
class Frame:
    """Maps instance coordinates onto an SVG canvas, y pointing up."""

    x0: float
    y1: float
    scale: float
    width: float
    height: float
    left: float
    top: float

    @staticmethod
    def fit(points: list[Point], size: float, margin: float, square: bool = False) -> "Frame":
        """The points' bounding box, its longer side `size` with the margins.

        The canvas hugs the box, or with `square` is `size` by `size` with
        the box centred in it.
        """
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        span_x = max(xs) - min(xs)
        span_y = max(ys) - min(ys)
        scale = (size - 2 * margin) / max(span_x, span_y)
        width = span_x * scale + 2 * margin
        height = span_y * scale + 2 * margin
        left = top = margin
        if square:
            left += (size - width) / 2
            top += (size - height) / 2
            width = height = size
        return Frame(min(xs), max(ys), scale, width, height, left, top)

    def __call__(self, p: Point) -> Point:
        return (
            self.left + (p[0] - self.x0) * self.scale,
            self.top + (self.y1 - p[1]) * self.scale,
        )


def svg(frame: Frame, title: str, body: list[str]) -> str:
    """A standalone SVG document with a transparent background."""
    w, h = f"{frame.width:.0f}", f"{frame.height:.0f}"
    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
            f'viewBox="0 0 {w} {h}" role="img">',
            f"<title>{title}</title>",
            *body,
            "</svg>",
            "",
        ]
    )


def tour_path(frame: Frame, coords: list[Point], tour: list[int], color: str, width: float) -> str:
    points = " ".join(f"{x:.1f},{y:.1f}" for x, y in (frame(coords[i]) for i in tour))
    return (
        f'<polyline points="{points}" fill="none" stroke="{color}" '
        f'stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round"/>'
    )


def depot_square(frame: Frame, p: Point, side: float, color: str) -> str:
    x, y = frame(p)
    return (
        f'<rect x="{x - side / 2:.1f}" y="{y - side / 2:.1f}" '
        f'width="{side:.1f}" height="{side:.1f}" fill="{color}"/>'
    )


def circle(frame: Frame, p: Point, r: float, color: str, opacity: float = 1.0) -> str:
    x, y = frame(p)
    alpha = "" if opacity == 1.0 else f' fill-opacity="{opacity}"'
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="{color}"{alpha}/>'


def draw_figure(coords: list[Point], profits: list[float], tour: list[int], pal: Palette) -> str:
    """Customers sized by profit, the optimal tour drawn, the rest faded."""
    frame = Frame.fit(coords, FIGURE_WIDTH, MARGIN)
    top = max(*profits[1:], 0.0) or 1.0

    def radius(i: int) -> float:
        return 2.5 + 7.5 * math.sqrt(max(profits[i], 0.0) / top)

    visited = set(tour)
    body = [
        circle(frame, coords[i], radius(i), pal.faded, 0.4)
        for i in range(1, len(coords))
        if i not in visited
    ]
    body.append(tour_path(frame, coords, tour, pal.accent, 2.5))
    body += [circle(frame, coords[i], radius(i), pal.accent) for i in tour if i != 0]
    body.append(depot_square(frame, coords[0], 12.0, pal.depot))
    return svg(frame, "The optimal tour of a CPTP instance", body)


def draw_logo(
    coords: list[Point], tour: list[int], accent: str, depot: str, size: float, icon: bool
) -> str:
    """The tour alone, cropped to its own nodes: the depot and its loop.

    An `icon` is square, with strokes heavy enough to read at 16 pixels.
    """
    stops = [coords[i] for i in tour]
    frame = Frame.fit(stops, size, size / (8 if icon else 10), square=icon)
    line, dot, side = (size / 10, size / 14, size / 5) if icon else (size / 16, size / 24, size / 8)
    body = [
        tour_path(frame, coords, tour, accent, line),
        *(circle(frame, coords[i], dot, accent) for i in tour if i != 0),
        depot_square(frame, coords[0], side, depot),
    ]
    return svg(frame, "cptp", body)


def main() -> None:
    problem = cptp.load(str(INSTANCE))
    model = cptp.Model()
    model.set_problem(problem)
    result = model.solve(OPTIONS)
    if not result.is_optimal():
        raise SystemExit(f"{INSTANCE.name}: no optimum ({result.status})")
    tour = [int(i) for i in result.tour]
    coords = read_coords(INSTANCE)
    profits = [float(p) for p in problem.profits]
    if len(coords) != problem.num_nodes:
        raise SystemExit(f"{INSTANCE.name}: {len(coords)} coordinates, {problem.num_nodes} nodes")

    for mode, pal in (("light", LIGHT), ("dark", DARK)):
        (HERE / f"tour-{mode}.svg").write_text(draw_figure(coords, profits, tour, pal))
        (HERE / f"logo-{mode}.svg").write_text(
            draw_logo(coords, tour, pal.accent, pal.depot, 128.0, icon=False)
        )
    (HERE / "favicon.svg").write_text(
        draw_logo(coords, tour, FAVICON_ACCENT, FAVICON_ACCENT, 64.0, icon=True)
    )
    print(
        f"{INSTANCE.name}: objective {result.objective:.3f}, "
        f"{len(tour) - 2} of {problem.num_nodes - 1} customers visited"
    )


if __name__ == "__main__":
    main()
