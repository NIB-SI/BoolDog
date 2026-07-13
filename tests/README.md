# BoolDog tests

To run all tests:

    cd tests
    python -m unittest discover .

## Coverage

Requires the `test` extra (`pip install booldog[test]`), or just `pip install coverage`.

    cd tests
    coverage run --source=booldog -m unittest discover .
    coverage report -m      # terminal summary
    coverage html            # detailed report in htmlcov/
