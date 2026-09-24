import pytest
from main import state


@pytest.fixture
def atom_dict():
    return state.atom_dict