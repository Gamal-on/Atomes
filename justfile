install:
    uv sync

partitionnement file="molecular_orbital":
    uv run -m partitionnement.{{file}} 

lint:
    uv run task lint

ty :
    uv run task ty

precommit:
    uv run task lint
    uv run task ty
