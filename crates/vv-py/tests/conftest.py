from pathlib import Path

import pytest
import vizviz

# tests/ -> vv-py -> crates -> repository root
FIXTURES = Path(__file__).resolve().parents[3] / "fixtures" / "small"


@pytest.fixture(scope="session")
def hemoglobin():
    return vizviz.load(FIXTURES / "4HHB.cif")


@pytest.fixture(scope="session")
def hemoglobin_pdb():
    return vizviz.load(FIXTURES / "4HHB.pdb")
