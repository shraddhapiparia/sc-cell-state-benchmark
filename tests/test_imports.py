"""Test imports."""

def test_imports():
    """Test that all modules import successfully."""
    try:
        import sc_cell_state_benchmark
        import sc_cell_state_benchmark.config
        import sc_cell_state_benchmark.data
        import sc_cell_state_benchmark.scoring
        import sc_cell_state_benchmark.controls
        import sc_cell_state_benchmark.evaluation
        import sc_cell_state_benchmark.plotting
        print("All imports successful.")
    except ImportError as e:
        raise AssertionError(f"Import failed: {e}")
