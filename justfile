install:
    uv sync

partitionnement file="molecular_orbital":
    uv run -m partitionnement.{{file}} 

precommit:
    uv run task lint
    uv run task ty
