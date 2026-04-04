"""Smoke test: verify all public package modules import without error."""


def test_imports():
    import sc_cell_state_benchmark
    import sc_cell_state_benchmark.config
    import sc_cell_state_benchmark.data
    import sc_cell_state_benchmark.scoring
    import sc_cell_state_benchmark.evaluation
    import sc_cell_state_benchmark.plotting
