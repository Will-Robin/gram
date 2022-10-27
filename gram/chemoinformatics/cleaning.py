from .core import cleanup as clean


def cleanup(reactions: list[str]) -> list[str]:
    return clean.cleanup(reactions)
