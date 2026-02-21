"""Test the Python RCSPP solver bindings."""

import pytest
import numpy as np


def test_model_basic():
    """Test basic model creation and solve."""
    import rcspp_bac

    model = rcspp_bac.Model()

    edges = np.array([[0, 1], [1, 0], [0, 2], [2, 0], [1, 2], [2, 1]], dtype=np.int32)
    costs = np.array([5.0, 5.0, 3.0, 3.0, 4.0, 4.0])

    model.set_graph(3, edges, costs)
    model.set_depot(0)
    model.set_profits(np.array([0.0, 10.0, 8.0]))
    model.add_capacity_resource(np.array([0.0, 2.0, 3.0]), 5.0)

    result = model.solve([("time_limit", "30")])

    assert result.has_solution()
    assert result.objective <= 0.0  # objective = travel_cost - profit, negative is good
    assert len(result.tour) > 0


def test_model_negative_costs():
    """Test model with negative edge costs."""
    import rcspp_bac

    model = rcspp_bac.Model()

    edges = np.array([[0, 1], [1, 0], [0, 2], [2, 0], [1, 2], [2, 1]], dtype=np.int32)
    costs = np.array([-2.0, 5.0, 3.0, 3.0, 4.0, -1.0])

    model.set_graph(3, edges, costs)
    model.set_depot(0)
    model.set_profits(np.array([0.0, 5.0, 5.0]))
    model.add_capacity_resource(np.array([0.0, 1.0, 1.0]), 10.0)

    result = model.solve([("time_limit", "30")])

    assert result.has_solution()
    assert result.objective < 0.0  # negative costs make the objective more negative
