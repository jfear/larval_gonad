def centerify(text, width=-1):
    """Center multiline text."""
    lines = text.split(" ")
    width = max(map(len, lines)) if width == -1 else width
    return "\n".join(line.center(width) for line in lines)
