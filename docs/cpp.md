# C++ reference

The user-facing part of the C++ library, in namespace `cptp`. Build a
`Problem` with `io::load` or through `Model`'s setters, and call
`Model::solve`. Its options are `SolverOptions`, a
`std::vector<std::pair<std::string, std::string>>` of the `(key, value)` pairs
listed in [Options](options.md).

## Model

```{doxygenclass} cptp::Model
```

```{doxygenstruct} cptp::Resource
```

## Problem

```{doxygenclass} cptp::Problem
```

```{doxygenstruct} cptp::Edge
```

## Results

```{doxygenstruct} cptp::SolveResult
```

```{doxygenstruct} cptp::SeparatorStats
```

## Loading instances

```{doxygennamespace} cptp::io
:content-only:
```
