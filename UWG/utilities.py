"""Collection of useful methods."""

def zeros(h, w):
    """create a matrix of zeros.

    Args:
        h: Height of the matrix.
        w: Width of the matrix.
    """
    return [[0 for x in range(w)] for y in range(h)]
