# Tests Directory

This directory contains validation tests for the pySCENIC pipeline.

## Test Files

- `test_pipeline.py` - Basic pipeline functionality tests
- `test_complete_pipeline.py` - End-to-end pipeline tests
- `test_all_cell_lines.py` - Multi-cell line validation
- `test_pipeline_simple.py` - Simplified validation tests

## Test Results

Test results are stored in the `test_results*/` subdirectories for reference.

## Running Tests

```bash
# Run basic tests
python tests/test_pipeline.py

# Run complete pipeline test
python tests/test_complete_pipeline.py
```

**Note**: These tests were used during development and validation. For production use, follow the main pipeline documentation in the README.
