---
# The link-preview and search description: the page's first paragraph, which
# sphinxext-opengraph would use, is the README's badge row.
og:description: >-
  cptp is a branch-and-cut solver for the Capacitated Profitable Tour Problem,
  with its open s-t path variant and several capacity rows: a C++ library, a
  command-line solver and Python bindings, on HiGHS.
myst:
  html_meta:
    description: >-
      cptp is a branch-and-cut solver for the Capacitated Profitable Tour
      Problem, with its open s-t path variant and several capacity rows: a C++
      library, a command-line solver and Python bindings, on HiGHS.
---

```{include} ../README.md
:relative-docs: docs/
:relative-images:
:end-before: <!-- docs/index.md puts the docs site's problem figure here. -->
```

% Drawn by figures/make_tour.py from a cptp solve; one copy per colour mode.

```{figure} figures/tour-light.svg
:figclass: only-light
:align: center
:width: 420px
:alt: Customers scattered around a depot, sized by profit; the optimal tour visits a group of them to one side of the depot, and the rest are faded.

The optimal tour of `E-n76-k7_a` (Roberti Set 3) from the depot (the square):
of 75 customers, sized by profit, it visits the 11 worth their travel cost
within the vehicle's capacity.
```

```{figure} figures/tour-dark.svg
:figclass: only-dark
:align: center
:width: 420px
:alt: Customers scattered around a depot, sized by profit; the optimal tour visits a group of them to one side of the depot, and the rest are faded.

The optimal tour of `E-n76-k7_a` (Roberti Set 3) from the depot (the square):
of 75 customers, sized by profit, it visits the 11 worth their travel cost
within the vehicle's capacity.
```

```{include} ../README.md
:relative-docs: docs/
:relative-images:
:start-after: <!-- docs/index.md puts the docs site's problem figure here. -->
```

```{toctree}
:hidden:

options
development
cpp
autoapi/cptp/index
```
