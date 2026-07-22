

import numpy as np

DEFAULT_BOX_SIZES = (1, 2, 4, 8, 16, 32, 50)


def box_count(binary_grid, box_size):
    n_rows, n_cols = binary_grid.shape
    n_boxes_rows = int(np.ceil(n_rows / box_size))
    n_boxes_cols = int(np.ceil(n_cols / box_size))

    padded = np.zeros((n_boxes_rows * box_size, n_boxes_cols * box_size), dtype=bool)
    padded[:n_rows, :n_cols] = np.asarray(binary_grid, dtype=bool)

    boxes = padded.reshape(n_boxes_rows, box_size, n_boxes_cols, box_size)
    return int(boxes.any(axis=(1, 3)).sum())


def box_counting_curve(binary_grid, box_sizes=DEFAULT_BOX_SIZES):
    binary_grid = np.asarray(binary_grid)
    sizes, counts = [], []
    for box_size in box_sizes:
        n = box_count(binary_grid, box_size)
        if n > 0:
            sizes.append(box_size)
            counts.append(n)
    return np.array(sizes, dtype=float), np.array(counts, dtype=float)


def box_counting_dimension(binary_grid, box_sizes=DEFAULT_BOX_SIZES):
    sizes, counts = box_counting_curve(binary_grid, box_sizes)
    if len(sizes) < 2:
        return np.nan

    log_inv_b = np.log(1.0 / sizes)
    log_n = np.log(counts)
    slope, _ = np.polyfit(log_inv_b, log_n, 1)
    return float(slope)
