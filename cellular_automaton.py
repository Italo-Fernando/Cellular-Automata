import numpy as np
from scipy.signal import convolve2d

STATE_S = 0
STATE_E = 1
STATE_I = 2
STATE_R = 3

GRID_SIZE = 100
MAX_DAYS = 500


VARIANT_PARAMETERS = {
    "original": {
        "label": "Original (Wuhan)",
        "p": 0.04,
        "sigma": 0.15,
        "gamma": 0.11,
        "r0": 2.79,
        "latent_period_days": 6.57,
        "infectious_period_days": 9.0,
    },
    "delta": {
        "label": "Delta",
        "p": 0.06,
        "sigma": 0.23,
        "gamma": 0.10,
        "r0": 5.08,
        "latent_period_days": 4.41,
        "infectious_period_days": 10.0,
    },
    "omicron": {
        "label": "Ômicron",
        "p": 0.23,
        "sigma": 0.29,
        "gamma": 0.19,
        "r0": 9.5,
        "latent_period_days": 3.42,
        "infectious_period_days": 5.16,
    },
}

# Vizinhança de Moore
_MOORE_KERNEL = np.array([
    [1, 1, 1],
    [1, 0, 1],
    [1, 1, 1],
])


class SEIRCellularAutomaton:

    def __init__(self, p, sigma, gamma, grid_size=GRID_SIZE, max_days=MAX_DAYS, rng=None):
        self.p = p
        self.sigma = sigma
        self.gamma = gamma
        self.grid_size = grid_size
        self.max_days = max_days
        self.rng = rng if rng is not None else np.random.default_rng()

        self.grid = np.full((grid_size, grid_size), STATE_S, dtype=np.int8)
        center = grid_size // 2
        self.grid[center, center] = STATE_I
        self.day = 0

    def _infectious_neighbor_counts(self):
        infectious_mask = (self.grid == STATE_I).astype(np.int8)
        return convolve2d(infectious_mask, _MOORE_KERNEL, mode="same", boundary="fill", fillvalue=0)

    def step(self):
        """Aplica as regras de transição simultaneamente a todas as células (1 passo = 1 dia)."""
        neighbor_counts = self._infectious_neighbor_counts()
        rand = self.rng.random(self.grid.shape)
        new_grid = self.grid.copy()

        infection_prob = 1.0 - (1.0 - self.p) ** neighbor_counts
        s_mask = self.grid == STATE_S
        new_grid[s_mask & (rand < infection_prob)] = STATE_E

        e_mask = self.grid == STATE_E
        new_grid[e_mask & (rand < self.sigma)] = STATE_I

        i_mask = self.grid == STATE_I
        new_grid[i_mask & (rand < self.gamma)] = STATE_R

        self.grid = new_grid
        self.day += 1

    def is_active(self):
        return bool(np.any((self.grid == STATE_E) | (self.grid == STATE_I)))

    def run(self):
        yield self.day, self.grid.copy()
        while self.is_active() and self.day < self.max_days:
            self.step()
            yield self.day, self.grid.copy()
