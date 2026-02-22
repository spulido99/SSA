.PHONY: test lint typecheck format install install-dev clean

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

install-all:
	pip install -e ".[all]"

test:
	pytest tests/ -v

test-fast:
	pytest tests/ -v -m "not slow and not integration and not gpu"

test-cov:
	pytest tests/ -v --cov=ssame --cov-report=term-missing

lint:
	ruff check ssame/ tests/

format:
	ruff format ssame/ tests/
	ruff check --fix ssame/ tests/

typecheck:
	mypy ssame/

clean:
	rm -rf build/ dist/ *.egg-info .pytest_cache .mypy_cache .ruff_cache htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
