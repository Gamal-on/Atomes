install:
    uv sync

partitionnement file="molecular_orbital":
    uv run -m partitionnement.{{file}} 
